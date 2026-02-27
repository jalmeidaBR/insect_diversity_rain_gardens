#####################################################################################################################
#######TESTES MODELOS

####CONJUNTO DE DADOS FINAL
models_data_final <- read.csv(
  "models_data_chapter_II.csv",
  sep = ",",
  fileEncoding = "latin1"
)
####ADICIONANDO OS VALORES DE BETA VEGETACAO
models_data_final <- models_data_final %>%
  left_join(
    explanatory_data %>%
      select(ID_RG, bray_combined, jaccard_native, jaccard_exotic),
    by = "ID_RG"
  )


library(bbmle)
library(MuMIn)
library(visreg)
library(DHARMa)
library(lme4)
library(effects)
library(MASS)
library(ggplot2)
library(visreg)
library(survival)
library(car)
library(spaMM)
library(GGally)
library(performance)
library(corrplot)
library(Hmisc)
library(ggeffects)
library(gridExtra)
library(dplyr)
library(merTools)
library(glmmTMB)


##################MODELOS DE PAISAGEM PARA ZE - INSETOS NOS JARDINS DE CHUVA######################

colnames(models_data_final)
exp_data_dredge <- models_data_final %>% dplyr::select(area_m2, IDHM_DIST, manag_level_RG, prop_veg_200, prop_veg_200.1, prop_veg_200.2, prop_edif_200, volume_mea_200, MESH_arb_200, MESH_herb_200,MESH_veg_200,
                                                       area_m2_z, volume_mea_200_z, MESH_arb_200_z, MESH_herb_200_z, MESH_veg_200_z)          

#Crie um data.frame com as variáveis preditoras
data_pred_na <- na.omit(exp_data_dredge)

#P4 Correlation analysis, gera uma matriz de correlação. Cuidado com modelos incluindo variáveis altamente correlacionadas (>|0.7|). Pode ignorar se tiver um X, porque é uma correlação não significativa 
cor.table <- cor(data_pred_na)
p.table <- rcorr(as.matrix(data_pred_na,type="pearson"))$P
corrplot::corrplot(cor(data_pred_na), method = 'color',addgrid.col = 'black',type = 'upper', diag = F, order = 'original',tl.col = 'black',
                   tl.srt = 90, tl.cex = 0.55, p.mat = p.table, sig.level = 0.05, insig = 'pch', pch.cex = 1.25,
                   addCoef.col='black', number.cex = 0.75, number.font = 2, number.digits = 2)
#par(mfrow = c(1, 1)) #limpar a divisão de plots

# P5 Fazendo uma matriz de correlação para o dredge

# Criar a função para testar correlação
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.6) {
  if (j >= i) return(NA)  # Evita calcular diagonal e valores duplicados
  
  # Seleciona as colunas corretamente
  x <- pull(data, i)  # Usa pull() para extrair coluna como vetor numérico
  y <- pull(data, j)
  
  # Teste de correlação
  ct <- cor.test(x, y, conf.level = conf.level)
  
  # Retorna TRUE se não for altamente correlacionado
  return(ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff)
}

# Vetoriza a função para uso no outer()
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Criar a matriz lógica de correlação
smat <- outer(1:ncol(data_pred_na), 1:ncol(data_pred_na), vCorrelated, data = data_pred_na)


nm <- colnames(data_pred_na)

dimnames(smat) <- list(nm, nm)

####MODELOS DISTRIBUIÇÃO POISSON - RIQUEZA ZE
#P5 Criar um modelo cheio, contendo todas as variaveis preditoras para escala 50m com morfotipos
options(na.action = "na.omit")

full_model <- glmer(alpha_geral_ZE ~ area_m2+ IDHM_DIST+ manag_level_RG+ prop_veg_200+ prop_veg_200.1+ prop_veg_200.2+ prop_edif_200+ volume_mea_200+ MESH_arb_200+ MESH_herb_200+MESH_veg_200 + esforco_PT_ZE + (1|ID_RG),
                    family = poisson,
                    data = models_data_final)
resid <- resid(full_model)
plot(resid)
summary(full_model)

###########
#teste de pressupostos 
library(lme4)
library(performance)   # Diagnóstico e qualidade do modelo
library(DHARMa)        # Simulação de resíduos
library(car)           # VIF
library(MuMIn)         # QAIC e dredge (se precisar)
library(see)           # Gráficos bonitos de diagnóstico (opcional)
library(glmmTMB)       # Caso queira testar NegBin



#-------------------------------------------
# 5️⃣ RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).


###ARRUMANDO O PROBLEMA NA ESCALA DAS VARIÁVEIS EXPLICATIVAS - REESCALANDO

vars_to_scale <- c("area_m2", "volume_mea_200", "MESH_arb_200", 
                   "MESH_herb_200", "MESH_veg_200")

models_data_final[ , paste0(vars_to_scale, "_z")] <- 
  scale(models_data_final[ , vars_to_scale])
###########################################################################
####modelo cheio com variaveis escalonadas
full_model <- glmer(alpha_geral_ZE ~ 
                      area_m2_z + IDHM_DIST + manag_level_RG +
                      prop_veg_200 + prop_veg_200.1 + prop_veg_200.2 +
                      prop_edif_200 + volume_mea_200_z +
                      MESH_arb_200_z + MESH_herb_200_z + MESH_veg_200_z +
                      esforco_PT_ZE + (1|ID_RG),
                    family = poisson,
                    data = models_data_final)


#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,12), fixed="esforco_PT_ZE" ,subset = smat)

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_alphainsetosZE_paisagem_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 1)[[1]])
summary(BEST_MODEL)

#########SELECIONOU O MODELO NULO

##########################################################################################################
#MODELOS ALPHA INSETOS ZSE
##aqui funcionou melhor sem usar OLRE
full_model <- glm(alpha_geral_ZSE ~ area_m2_z+ IDHM_DIST+ prop_veg_200+ prop_veg_200.1+ prop_veg_200.2+ prop_edif_200+ volume_mea_200_z+ MESH_arb_200_z+ MESH_herb_200_z+MESH_veg_200_z + esforco_PT_ZSE ,
                  family = poisson,
                  data = models_data_final) 
resid <- resid(full_model)
plot(resid)
summary(full_model)

###########
#teste de pressupostos 
#-------------------------------------------
# RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).


#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,11), fixed="esforco_PT_ZSE" ,subset = smat) #diminui de 12 pra 11 pq tirei manejo

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models,delta<2) #2 modelos
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_alphainsetosZSE_paisagem_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 2)[[1]])
summary(BEST_MODEL)

#modelo 1
model_1 <- glm(formula = alpha_geral_ZSE ~ area_m2 + IDHM_DIST + prop_veg_200.1 + 
                 volume_mea_200 + esforco_PT_ZSE, family = poisson, 
               data = models_data_final)
summary(model_1)
#modelo 2
model_2 <- glm(formula = alpha_geral_ZSE ~ IDHM_DIST + prop_veg_200.1 + 
                 volume_mea_200 + esforco_PT_ZSE, family = poisson, 
               data = models_data_final)
summary(model_2)



#####MODELOS ALPHA POLINIZADORES ZE
full_model_pol<- glmmTMB(alpha_polinizadores_ZE ~ area_m2_z+ IDHM_DIST+ manag_level_RG+ prop_veg_200+ prop_veg_200.1+ prop_veg_200.2+ prop_edif_200+ volume_mea_200_z+ MESH_arb_200+ MESH_herb_200_z + MESH_veg_200_z + esforco_PT_ZE, family = compois(link = "log"), data = models_data_final)

res <- simulateResiduals(full_model_pol_g)
plot(res)  # Diagnóstico gráfico de resíduos simulados

###########
#teste de pressupostos 
#-------------------------------------------
# RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model_pol)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).


#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model_pol, m.lim = c(0,12), fixed="cond(esforco_PT_ZE)" ,subset = smat) 

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models, delta < 2 & !is.na(delta)) 
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_alphapol_ZE_paisagem_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 2)[[1]])
summary(BEST_MODEL)

##########################################################################################################
#####MODELOS POLINIZADORES ZSE

full_model_pol_ZSE<- glmmTMB(alpha_polinizadores_ZSE ~ area_m2_z+ IDHM_DIST+ prop_veg_200+ prop_veg_200.1+ prop_veg_200.2+ prop_edif_200+ volume_mea_200_z+ MESH_arb_200+ MESH_herb_200_z + MESH_veg_200_z + esforco_PT_ZSE, family = compois(link="log"), data = models_data_final)

#teste de pressupostos 
#-------------------------------------------
# RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model_pol_ZSE)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).


#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model_pol_ZSE, m.lim = c(0,11), fixed="cond(esforco_PT_ZSE)" ,subset = smat) 

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models, delta < 2 & !is.na(delta)) 
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_alphapol_ZSE_paisagem_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 2)[[1]])
summary(BEST_MODEL)

##############################################################################################################MODELOS COM COMPOSICAO VEGETAL

colnames(models_data_final)
exp_data_dredge <- models_data_final %>% dplyr::select(area_m2_z, manag_level_RG, native_RG,exotic_RG,      riqueza_jardim, riqueza_invasive_RG, garbage_RG, combined_simpson_inside, simpson_inside_native,   simpson_inside_exotic)          

#Crie um data.frame com as variáveis preditoras
data_pred_na <- na.omit(exp_data_dredge)

#P4 Correlation analysis, gera uma matriz de correlação. Cuidado com modelos incluindo variáveis altamente correlacionadas (>|0.7|). Pode ignorar se tiver um X, porque é uma correlação não significativa 
cor.table <- cor(data_pred_na)
p.table <- rcorr(as.matrix(data_pred_na,type="pearson"))$P
corrplot::corrplot(cor(data_pred_na), method = 'color',addgrid.col = 'black',type = 'upper', diag = F, order = 'original',tl.col = 'black',
                   tl.srt = 90, tl.cex = 0.55, p.mat = p.table, sig.level = 0.05, insig = 'pch', pch.cex = 1.25,
                   addCoef.col='black', number.cex = 0.75, number.font = 2, number.digits = 2)
#par(mfrow = c(1, 1)) #limpar a divisão de plots

# P5 Fazendo uma matriz de correlação para o dredge

# Criar a função para testar correlação
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.6) {
  if (j >= i) return(NA)  # Evita calcular diagonal e valores duplicados
  
  # Seleciona as colunas corretamente
  x <- pull(data, i)  # Usa pull() para extrair coluna como vetor numérico
  y <- pull(data, j)
  
  # Teste de correlação
  ct <- cor.test(x, y, conf.level = conf.level)
  
  # Retorna TRUE se não for altamente correlacionado
  return(ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff)
}

# Vetoriza a função para uso no outer()
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Criar a matriz lógica de correlação
smat_JC <- outer(1:ncol(data_pred_na), 1:ncol(data_pred_na), vCorrelated, data = data_pred_na)


nm <- colnames(data_pred_na)

dimnames(smat_JC) <- list(nm, nm)

####MODELOS DISTRIBUIÇÃO POISSON - RIQUEZA ZE insetos
#P5 Criar um modelo cheio, contendo todas as variaveis preditoras para escala 50m com morfotipos
options(na.action = "na.omit")

full_model_veg <- glmer(alpha_geral_ZE ~ area_m2_z+ manag_level_RG+ native_RG+exotic_RG+ riqueza_jardim+ riqueza_invasive_RG+ garbage_RG + combined_simpson_inside+ simpson_inside_native+ simpson_inside_exotic + esforco_PT_ZE + (1|ID_RG),
                        family = poisson,
                        data = models_data_final) 
resid <- resid(full_model_veg)
plot(resid)
summary(full_model_veg)

#-------------------------------------------
# RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model_veg)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).


#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model_veg, m.lim = c(0,11), fixed="esforco_PT_ZE" ,subset = smat_JC) 

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models, delta < 2 & !is.na(delta)) 
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_alphainsetos_ZE_veg_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 2)[[1]])
summary(BEST_MODEL)


####MODELOS DISTRIBUIÇÃO POISSON - RIQUEZA ZE insetos
#P5 Criar um modelo cheio, contendo todas as variaveis preditoras para escala 50m com morfotipos
options(na.action = "na.omit")

full_model_veg_pol <- glmmTMB(alpha_polinizadores_ZE ~ area_m2_z+ manag_level_RG+ native_RG+exotic_RG+ riqueza_jardim+ riqueza_invasive_RG+ garbage_RG+ combined_simpson_inside+ simpson_inside_native+ simpson_inside_exotic + esforco_PT_ZE ,
                              family = compois(link="log"),
                              data = models_data_final)

AIC(full_model_veg_pol, full_model_veg_pol_2)

resid <- resid(full_model_veg_pol)
plot(resid)
summary(full_model_veg_pol)

#-------------------------------------------
# RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model_veg_pol)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).


#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model_veg_pol, m.lim = c(0,11), fixed="cond(esforco_PT_ZE)" ,subset = smat_JC) 

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models, delta < 2 & !is.na(delta)) 
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_alphapol_ZSE_veg_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 2)[[1]])
summary(BEST_MODEL)

###########
##GAMMA PAISAGEM 

#INSETOS
####modelo cheio com variaveis escalonadas
full_model <- glmer(gama_geral ~ 
                      area_m2_z + IDHM_DIST + manag_level_RG +
                      prop_veg_200 + prop_veg_200.1 + prop_veg_200.2 +
                      prop_edif_200 + volume_mea_200_z +
                      MESH_arb_200_z + MESH_herb_200_z + MESH_veg_200_z +
                      esforco_geral_PT + (1|ID_RG),
                    family = poisson,
                    data = models_data_final)

#RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).



#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,12), fixed="esforco_geral_PT" ,subset = smat)

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_gama_insetos_paisagem_.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 2)[[1]])
summary(BEST_MODEL)

##########################################################################################################
##POLINIZADORES
####modelo cheio com variaveis escalonadas
full_model <- glm(gama_polinizadores ~ 
                    area_m2_z + IDHM_DIST + manag_level_RG +
                    prop_veg_200 + prop_veg_200.1 + prop_veg_200.2 +
                    prop_edif_200 + volume_mea_200_z +
                    MESH_arb_200_z + MESH_herb_200_z + MESH_veg_200_z +
                    esforco_geral_PT,
                  family = poisson,
                  data = models_data_final)

#RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).



#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,12), fixed="esforco_geral_PT" ,subset = smat)

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_gama_pol_paisagem.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 1)[[1]])
summary(BEST_MODEL)

##########################################################################################################
##MODELOS BETA

##INSETOS
##Correcao dos valores diferentes de 0 e 1
models_data_final$beta_geral_cor<- (models_data_final$beta_geral * (30 - 1) + 0.5) / 30


####modelo cheio com variaveis escalonadas
full_model <- glmmTMB(beta_geral_cor ~ 
                        area_m2_z + IDHM_DIST + manag_level_RG +
                        prop_veg_200 + prop_veg_200.1 + prop_veg_200.2 +
                        prop_edif_200 + volume_mea_200_z +
                        MESH_arb_200_z + MESH_herb_200_z + MESH_veg_200_z +
                        esforco_geral_PT + (1|ID_RG),
                      family = beta_family(link = "logit"),
                      data = models_data_final)

#RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).



#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,12), fixed="cond(esforco_geral_PT)" ,subset = smat)

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models, delta < 2 & !is.na(delta)) 
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_beta_insetos_paisagem.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 7)[[1]])
summary(BEST_MODEL)

##########################################################################################################
##POLINIZADORES
##Correcao dos valores diferentes de 0 e 1
models_data_final$beta_polinizadores_cor<- (models_data_final$beta_polinizadores * (30 - 1) + 0.5) / 30

# Para remover os NA para os polinizadores
models_data_final_pol <- models_data_final[!is.na(models_data_final$beta_polinizadores),]

full_model_beta <- glmmTMB(beta_polinizadores_cor ~ 
                             area_m2_z + IDHM_DIST + manag_level_RG +
                             prop_veg_200 + prop_veg_200.1 + prop_veg_200.2 +
                             prop_edif_200 + volume_mea_200_z +
                             MESH_arb_200_z + MESH_herb_200_z + MESH_veg_200_z +
                             esforco_geral_PT + (1|ID_RG),
                           family = beta_family("cloglog"),
                           data = models_data_final_pol, 
                           control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

performance::check_collinearity(full_model_beta)


#RESÍDUOS (DHARMa)
#-------------------------------------------
res <- simulateResiduals(full_model_beta)
plot(res)  # Diagnóstico gráfico de resíduos simulados
#QQPlot - Se os pontos seguem a linha diagonal vermelha, o modelo ajusta bem aos dados.

testQuantiles(res)      # Testa formalmente os desvios de quantis
testResiduals(res)      # Versão geral de todos os testes DHARMa
plotResiduals(res)      # Resíduos vs. variáveis específicas

# Testes formais:
testDispersion(res)     # confirma sobredispersão via simulação - p > 0.05 → sem sobredispersão ou subdispersão.
testZeroInflation(res)  # confirma inflação de zeros - p-valor seria < 0.05, sugere o uso de um modelo zero-inflated.
testUniformity(res)     # Interpretação: resíduos estão uniformemente distribuídos (p > 0.05).



#P7 Fazer a seleção pelo método dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model_beta, m.lim = c(0,12), fixed="cond(esforco_geral_PT)" ,subset = smat)

#Selecionar os modelos com delta AIC menor que 2
candidate_models_sub <- subset(candidate_models, delta < 2 & !is.na(delta)) 
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Criar um novo dataframe com os melhores modelos (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

write.csv(best_models_df, "dredge_beta_pol_gamma_paisagem.csv", row.names = FALSE)

#####SALVAR TABELA COM AS VARIÁVEIS

### Gráfico mostrando a tabela de seleção pelo dredge
x11()
par(mar = c(3,5,6,4))
plot(candidate_models_sub, labAsExpr = T)

#P8 SUMMARY do melhor modelo. Trocar o primeiro 1 para 2 ou 3 dependendo de quantos melhores modelos tiverem, para ver todos
BEST_MODEL <- (get.models(candidate_models, 4)[[1]])
summary(BEST_MODEL)


##############################################
######GRAFICOS MODELOS
setwd("DIRECTORY")
models_data_2 <- read.csv("models_data_chapter_II_zvar.csv", sep=",")
library(lme4)
library(ggeffects)
library(patchwork)
library(tidyr)
library(gamlss)
library(vegan)
library(MuMIn)
library(pscl)
library(performance)
library(ggplot2)
#install.packages("patchwork")
#install.packages("gamlss")

### Model alpha pollinators ZE

model_alpha_pol_ZE_1 <- glmmTMB(alpha_polinizadores_ZE ~ prop_veg_200.1 + MESH_arb_200+ esforco_PT_ZE, 
                                family = compois(link = "log"), data = models_data_2)
summary(model_alpha_pol_ZE_1)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(model_alpha_pol_ZE_1))
null_ll  <- as.numeric(logLik(update(model_alpha_pol_ZE_1, . ~ 1)))
k        <- length(fixef(model_alpha_pol_ZE_1)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

model_alpha_pol_ZE_2 <- glmmTMB(alpha_polinizadores_ZE ~ prop_veg_200.1 +esforco_PT_ZE, 
                                family = compois(link = "log"), data = models_data_2)
summary(model_alpha_pol_ZE_2)

################################################
full_ll  <- as.numeric(logLik(model_alpha_pol_ZE_2))
null_ll  <- as.numeric(logLik(update(model_alpha_pol_ZE_2, . ~ 1)))
k        <- length(fixef(model_alpha_pol_ZE_2)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

### Graphic
# Predictions
pred_alpha_pol_ZE <- ggpredict(model_alpha_pol_ZE_1, terms = c("prop_veg_200.1[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_alpha_pol_ZE <- min(pred_alpha_pol_ZE$conf.low, na.rm = TRUE)
y_max_alpha_pol_ZE <- max(pred_alpha_pol_ZE$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),
    #axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_alpha_pol_ZE <- ggplot(pred_alpha_pol_ZE, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Proportion of tree cover", y = "Alpha diversity (effect zone)") +
  my_theme +
  scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_alpha_pol_ZE, y_max_alpha_pol_ZE)

plot(plot_alpha_pol_ZE)

ggsave("DIRECTORY/plot_alpha_pol_ZE.png", 
       plot = plot_alpha_pol_ZE, 
       width = 12, height = 9.5, dpi = 300)


################MODELS ALPHA INSECTS ZSE
alpha_insetos_ZSE_1<-glm(alpha_geral_ZSE ~ area_m2+ IDHM_DIST+ prop_veg_200.1+volume_mea_200+ esforco_PT_ZSE ,
                         family = poisson,
                         data = models_data_2)
summary(alpha_insetos_ZSE_1)

####calcular R²
R2 <- 1 - alpha_insetos_ZSE_1$deviance / alpha_insetos_ZSE_1$null.deviance

n <- nobs(alpha_insetos_ZSE_1)
p <- length(coef(alpha_insetos_ZSE_1)) - 1   # número de preditores

R2_adj <- 1 - (1 - R2) * ((n - 1) / (n - p - 1))
R2_adj
###################

alpha_insetos_ZSE_2<-glm(alpha_geral_ZSE ~ IDHM_DIST+ prop_veg_200.1+volume_mea_200_z+ esforco_PT_ZSE ,
                         family = poisson,
                         data = models_data_2)
summary(alpha_insetos_ZSE_2)

####calcular R²
R2 <- 1 - alpha_insetos_ZSE_2$deviance / alpha_insetos_ZSE_2$null.deviance

n <- nobs(alpha_insetos_ZSE_2)
p <- length(coef(alpha_insetos_ZSE_2)) - 1   # número de preditores

R2_adj <- 1 - (1 - R2) * ((n - 1) / (n - p - 1))
R2_adj
###################


### Graphic
# Predictions
pred_area <- ggpredict(alpha_insetos_ZSE_1, terms = c("area_m2[all]"))

pred_idh  <- ggpredict(alpha_insetos_ZSE_1, terms = c("IDHM_DIST[all]"))

pred_prop_1  <- ggpredict(alpha_insetos_ZSE_1, terms = c("prop_veg_200.1[all]"))

pred_vol  <- ggpredict(alpha_insetos_ZSE_1, terms = c("volume_mea_200[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_alpha_insect_ZSE <- min(pred_area$conf.low, pred_idh $conf.low, pred_prop_1$conf.low, pred_vol$conf.low,  na.rm = TRUE)
y_max_alpha_insect_ZSE <- max(pred_area$conf.high, pred_idh $conf.high, pred_prop_1$conf.high, pred_vol$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_area <- ggplot(pred_area, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Rain garden area (m²)", y = "Insect alpha diversity (no-effect zone)") +
  my_theme +
  ylim(y_min_alpha_insect_ZSE, y_max_alpha_insect_ZSE)

plot(plot_area)

# Second term
plot_IDH <- ggplot(pred_idh, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Income per capita", y = "Insect alpha diversity (no-effect zone)") +
  my_theme +
  ylim(y_min_alpha_insect_ZSE, y_max_alpha_insect_ZSE)

plot(plot_IDH)

#Third term
plot_prop_arb <- ggplot(pred_prop_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Proportion of tree cover", y = "Insect alpha diversity (no-effect zone)") +
  my_theme +
  scale_x_continuous(n.breaks = 6) +
  ylim(y_min_alpha_insect_ZSE, y_max_alpha_insect_ZSE)

plot(plot_prop_arb)

#4 term
plot_vol <- ggplot(pred_vol, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Average building volume (m³)", y = "Insect alpha diversity (no-effect zone)") +
  my_theme +
  ylim(y_min_alpha_insect_ZSE, y_max_alpha_insect_ZSE)

plot(plot_vol)

# Combined graphic
plot_alpha_insects_ZSE<- (plot_area + plot_IDH + plot_prop_arb + plot_vol) + 
  plot_layout(ncol = 2, nrow = 2, axis_titles = "collect_y") 

print(plot_alpha_insects_ZSE)


ggsave("DIRECTORY/plot_alpha_insects_ZSE_v2.png", 
       plot = plot_alpha_insects_ZSE, 
       width = 13, height = 11, dpi = 300)




##GAMMA PAISAGEM 
##POLINIZADORES
####modelo cheio com variaveis escalonadas
model_gamma_pol <- glm(gama_polinizadores ~ prop_veg_200.1 + esforco_geral_PT,
                       family = poisson,
                       data = models_data_2)

summary(model_gamma_pol)

####calcular R²
R2 <- 1 - model_gamma_pol$deviance / model_gamma_pol$null.deviance

n <- nobs(model_gamma_pol)
p <- length(coef(model_gamma_pol)) - 1   # número de preditores

R2_adj <- 1 - (1 - R2) * ((n - 1) / (n - p - 1))
R2_adj
###################

### Graphic
# Predictions
pred_gamma_pol <- ggpredict(model_gamma_pol, terms = c("prop_veg_200.1[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_gamma_pol <- min(pred_gamma_pol$conf.low, na.rm = TRUE)
y_max_gamma_pol <- max(pred_gamma_pol$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_gamma_pol <- ggplot(pred_gamma_pol, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Proportion of tree cover", y = "Gamma diversity") +
  my_theme +
  scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_gamma_pol, y_max_gamma_pol)

plot(plot_gamma_pol)

ggsave("DIRECTORY/plot_gamma_pol.png", 
       plot = plot_gamma_pol, 
       width = 12, height = 9.5, dpi = 300)


#GAMMA INSETOS
gamma_insects_1 <- glmer(gama_geral ~ 
                           prop_veg_200.1 + volume_mea_200 +
                           esforco_geral_PT + (1|ID_RG),
                         family = poisson,
                         data = models_data_2)
summary(gamma_insects_1)


####calcular R²
R2m <- r.squaredGLMM(gamma_insects_1)[1]  # marginal R²
R2m
n <- nobs(gamma_insects_1)
p <- length(fixef(gamma_insects_1)) - 1

R2_adj <- 1 - (1 - R2m) * ((n - 1) / (n - p - 1))
R2_adj
###################

gamma_insects_2 <- glmer(gama_geral ~ area_m2_z +
                           prop_veg_200.1 + volume_mea_200_z +
                           esforco_geral_PT + (1|ID_RG),
                         family = poisson,
                         data = models_data_2)
summary(gamma_insects_2)


####calcular R²
R2m <- r.squaredGLMM(gamma_insects_2)[1]  # marginal R²
R2m
n <- nobs(gamma_insects_2)
p <- length(fixef(gamma_insects_2)) - 1

R2_adj <- 1 - (1 - R2m) * ((n - 1) / (n - p - 1))
R2_adj
###################

### Graphic
# Predictions
pred_prop_arb<- ggpredict(gamma_insects_1, terms = c("prop_veg_200.1[all]"))

pred_volume  <- ggpredict(gamma_insects_1, terms = c("volume_mea_200[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_gamma_insect <- min(pred_prop_arb$conf.low, pred_volume $conf.low, na.rm = TRUE)
y_max_gamma_insect <- max(pred_prop_arb$conf.high, pred_volume $conf.high,  na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_arb <- ggplot(pred_prop_arb, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Proportion of tree cover", y = "Insect gamma diversity") +
  my_theme +
  scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_gamma_insect, y_max_gamma_insect)

plot(plot_arb)

# Second term
plot_volume <- ggplot(pred_volume, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Average building volume", y = "Insect gamma diversity") +
  my_theme +
  ylim(y_min_gamma_insect, y_max_gamma_insect)

plot(plot_volume)

# Combined graphic
plot_insect_gamma<- (plot_arb + plot_volume) + 
  plot_layout(ncol = 2, axis_titles = "collect_y" ) 

print(plot_insect_gamma)


ggsave("DIRECTORY/plot_insect_gamma.png", 
       plot = plot_insect_gamma, 
       width = 13, height = 11, dpi = 300)


##MODELOS BETA

##INSETOS
##Correcao dos valores diferentes de 0 e 1
models_data_2$beta_geral_cor<- (models_data_2$beta_geral * (30 - 1) + 0.5) / 30


####modelo beta insetos
beta_insetos_1 <- glmmTMB(beta_geral_cor ~ IDHM_DIST + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_1)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_1))
null_ll  <- as.numeric(logLik(update(beta_insetos_1, . ~ 1)))
k        <- length(fixef(beta_insetos_1)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_2 <- glmmTMB(beta_geral_cor ~ IDHM_DIST + MESH_herb_200_z + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_2)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_2))
null_ll  <- as.numeric(logLik(update(beta_insetos_2, . ~ 1)))
k        <- length(fixef(beta_insetos_2)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_3 <- glmmTMB(beta_geral_cor ~ IDHM_DIST + MESH_arb_200_z + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_3)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_3))
null_ll  <- as.numeric(logLik(update(beta_insetos_3, . ~ 1)))
k        <- length(fixef(beta_insetos_3)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_4 <- glmmTMB(beta_geral_cor ~ IDHM_DIST +  prop_veg_200.2 + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_4)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_4))
null_ll  <- as.numeric(logLik(update(beta_insetos_4, . ~ 1)))
k        <- length(fixef(beta_insetos_4)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_5 <- glmmTMB(beta_geral_cor ~ IDHM_DIST +  prop_veg_200.1 + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_5)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_5))
null_ll  <- as.numeric(logLik(update(beta_insetos_5, . ~ 1)))
k        <- length(fixef(beta_insetos_5)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_6 <- glmmTMB(beta_geral_cor ~ IDHM_DIST +  MESH_arb_200 + MESH_veg_200 + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_6)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_6))
null_ll  <- as.numeric(logLik(update(beta_insetos_6, . ~ 1)))
k        <- length(fixef(beta_insetos_6)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_7 <- glmmTMB(beta_geral_cor ~ IDHM_DIST +  MESH_arb_200_z + prop_veg_200.2 + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_7)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_7))
null_ll  <- as.numeric(logLik(update(beta_insetos_7, . ~ 1)))
k        <- length(fixef(beta_insetos_7)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_insetos_8 <- glmmTMB(beta_geral_cor ~ IDHM_DIST +  MESH_arb_200_z + MESH_herb_200_z + esforco_geral_PT + (1|ID_RG),
                          family = beta_family(link = "logit"),
                          data = models_data_2)
summary(beta_insetos_8)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_insetos_8))
null_ll  <- as.numeric(logLik(update(beta_insetos_8, . ~ 1)))
k        <- length(fixef(beta_insetos_8)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

### Graphic
# Predictions
pred_beta_insect_idh <- ggpredict(beta_insetos_1, terms = c("IDHM_DIST[all]"))
pred_beta_insect_mesh_arb <- ggpredict(beta_insetos_6, terms = c("MESH_arb_200[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_beta_insect <- min(pred_beta_insect_idh$conf.low, pred_beta_insect_mesh_arb$conf.low, na.rm = TRUE)
y_max_beta_insect <- max(pred_beta_insect_idh$conf.high, pred_beta_insect_mesh_arb$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_beta_insect_idh <- ggplot(pred_beta_insect_idh, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Income per capita", y = "Insect beta diversity") +
  my_theme +
  #scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_beta_insect, y_max_beta_insect)

plot(plot_beta_insect_idh)

# Second term
plot_beta_insect_mesh<- ggplot(pred_beta_insect_mesh_arb, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Effective mesh size of arboreal cover", y = "Insect beta diversity") +
  my_theme +
  scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_beta_insect, y_max_beta_insect)

plot(plot_beta_insect_mesh)

# Combined graphic
plot_insect_beta<- (plot_beta_insect_idh + plot_beta_insect_mesh) + 
  plot_layout(ncol = 2, axis_titles = "collect_y" ) 

print(plot_insect_beta)

ggsave("DIRECTORY/plot_insect_beta.png", 
       plot = plot_insect_beta, 
       width = 12.55, height = 9.5, dpi = 300)


##POLINIZADORES
##Correcao dos valores diferentes de 0 e 1
models_data_2$beta_polinizadores_cor<- (models_data_2$beta_polinizadores * (30 - 1) + 0.5) / 30

# Para remover os NA para os polinizadores
models_data_final_pol <- models_data_2[!is.na(models_data_2$beta_polinizadores),]

beta_pol_1 <- glmmTMB(beta_polinizadores_cor ~ MESH_herb_200+
                        esforco_geral_PT + (1|ID_RG),
                      family = beta_family("cloglog"),
                      data = models_data_final_pol, 
                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

summary(beta_pol_1)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_pol_1))
null_ll  <- as.numeric(logLik(update(beta_pol_1, . ~ 1)))
k        <- length(fixef(beta_pol_1)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_pol_2 <- glmmTMB(beta_polinizadores_cor ~ prop_veg_200.2 +
                        esforco_geral_PT + (1|ID_RG),
                      family = beta_family("cloglog"),
                      data = models_data_final_pol, 
                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

summary(beta_pol_2)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_pol_2))
null_ll  <- as.numeric(logLik(update(beta_pol_2, . ~ 1)))
k        <- length(fixef(beta_pol_2)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_pol_3 <- glmmTMB(beta_polinizadores_cor ~ MESH_herb_200_z + MESH_veg_200_z +
                        esforco_geral_PT + (1|ID_RG),
                      family = beta_family("cloglog"),
                      data = models_data_final_pol, 
                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

summary(beta_pol_3)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_pol_3))
null_ll  <- as.numeric(logLik(update(beta_pol_3, . ~ 1)))
k        <- length(fixef(beta_pol_3)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

beta_pol_4 <- glmmTMB(beta_polinizadores_cor ~ MESH_veg_200_z + prop_veg_200.2 +
                        esforco_geral_PT + (1|ID_RG),
                      family = beta_family("cloglog"),
                      data = models_data_final_pol, 
                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

summary(beta_pol_4)

#######################R2 McFadden ajustado#########################
full_ll  <- as.numeric(logLik(beta_pol_4))
null_ll  <- as.numeric(logLik(update(beta_pol_4, . ~ 1)))
k        <- length(fixef(beta_pol_4)$cond)

R2_adj_McFadden <- 1 - ((full_ll - k) / null_ll)
R2_adj_McFadden
#########################################################

### Graphic
# Predictions
pred_beta_pol_herb <- ggpredict(beta_pol_2, terms = c("prop_veg_200.2[all]"))
pred_beta_pol_mesh_herb <- ggpredict(beta_pol_1, terms = c("MESH_herb_200[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_beta_pol <- min(pred_beta_pol_herb$conf.low, pred_beta_pol_mesh_herb$conf.low, na.rm = TRUE)
y_max_beta_pol <- max(pred_beta_pol_herb$conf.high, pred_beta_pol_mesh_herb$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_beta_pol_prop_herb <- ggplot(pred_beta_pol_herb, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Proportion of herb–shrub cover", y = "Beta diversity") +
  my_theme +
  #scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_beta_pol, y_max_beta_pol)

plot(plot_beta_pol_prop_herb)

# Second term
plot_beta_pol_mesh<- ggplot(pred_beta_pol_mesh_herb, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Effective mesh size of herb-shrub cover", y = "Beta diversity") +
  my_theme +
  #scale_x_continuous(n.breaks = 6) +   # <<<<<< AQUI
  ylim(y_min_beta_pol, y_max_beta_pol)

plot(plot_beta_pol_mesh)

# Combined graphic
plot_pol_beta<- (plot_beta_pol_prop_herb + plot_beta_pol_mesh) + 
  plot_layout(ncol = 2, axis_titles = "collect_y" ) 

print(plot_pol_beta)

ggsave("DIRECTORY/plot_pol_beta.png", 
       plot = plot_pol_beta, 
       width = 13, height = 9.5, dpi = 300)


####COMBINAR GRAFICOS POLINIZADORES
PLOT_POLINIZADORES <- ( plot_alpha_pol_ZE +plot_gamma_pol + plot_beta_pol_prop_herb + plot_beta_pol_mesh)+ 
  plot_layout(ncol = 2, nrow = 2, heights = c(1, 1)) 
print(PLOT_POLINIZADORES)

ggsave("DIRECTORY/plot_pol_todos_v2.png", 
       plot = PLOT_POLINIZADORES, 
       width = 15, height = 11, dpi = 300)




#############MODELOS VARIAVEIS JARDINS
model_veg_1 <- glmer(alpha_geral_ZE ~ garbage_RG  + esforco_PT_ZE + (1|ID_RG),
                     family = poisson,
                     data = models_data_2) 
summary(model_veg_1)
####calcular R²
R2m <- r.squaredGLMM(model_veg_1)[1]  # marginal R²
R2m
n <- nobs(model_veg_1)
p <- length(fixef(model_veg_1)) - 1

R2_adj <- 1 - (1 - R2m) * ((n - 1) / (n - p - 1))
R2_adj
###################

model_veg_2 <- glmer(alpha_geral_ZE ~ garbage_RG + simpson_inside_native + esforco_PT_ZE + (1|ID_RG),
                     family = poisson,
                     data = models_data_2) 

summary(model_veg_2)
####calcular R²
R2m <- r.squaredGLMM(model_veg_2)[1]  # marginal R²
R2m
n <- nobs(model_veg_2)
p <- length(fixef(model_veg_2)) - 1

R2_adj <- 1 - (1 - R2m) * ((n - 1) / (n - p - 1))
R2_adj
###################

### Graphic
# Predictions
pred_veg <- ggpredict(model_veg_2, terms = c("garbage_RG[all]"))

#FUNÇÃO PARA AJUSTAR O VALOR DO EIXO Y EM FUNÇÃO DO INTERVALO DE CONFIANÇA
y_min_veg <- min(pred_veg$conf.low, na.rm = TRUE)
y_max_veg <- max(pred_veg$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_veg <- ggplot(pred_veg, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Presence of anthropogenic litter", y = "Insect alpha diversity (effect zone)") +
  my_theme +
  scale_x_continuous(breaks = c(0, 1)) +
  ylim(y_min_veg, y_max_veg)

plot(plot_veg)

ggsave("DIRECTORY/plot_alpha_insect_veg_ZE_v2.png", 
       plot = plot_veg, 
       width = 12, height = 9.5, dpi = 300)


#####################################################################################################################################################################################################################NMDS

# Pacotes necessários
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(vegan)

setwd("DIRECTORY")
# --------------------------------------------------
# 1. Ler os dados
# --------------------------------------------------
dados <- read.csv("dados_NMDS.csv", sep=";", fileEncoding = "latin1")

# --------------------------------------------------
# 2. Criar o identificador do site (ID_JC + zona)
# --------------------------------------------------
dados <- dados %>%
  mutate(site = paste(ID_JC, zona, sep = "_"))

# --------------------------------------------------
# 3. Função para criar matriz site x morfoespécie
#    (presença / ausência)
# --------------------------------------------------
criar_matriz_pa <- function(df) {
  df %>%
    distinct(site, morfoespecie) %>%   # evita duplicações
    mutate(presenca = 1) %>%
    pivot_wider(
      names_from  = morfoespecie,
      values_from = presenca,
      values_fill = 0
    )%>%
    as.data.frame() %>%
    `rownames<-`(.$site) %>%
    select(-site)
}

criar_matriz_abund <- function(df) {
  df %>%
    count(site, morfoespecie) %>%   # conta indivíduos
    pivot_wider(
      names_from  = morfoespecie,
      values_from = n,
      values_fill = 0
    ) %>%
    as.data.frame() %>%
    `rownames<-`(.$site) %>%
    select(-site)
}
# --------------------------------------------------
# 4. Conjunto de dados: TODOS os insetos
# --------------------------------------------------
comm_all <- criar_matriz_pa(dados)
comm_all <- criar_matriz_abund(dados)
# --------------------------------------------------
# 5. Conjunto de dados: APENAS polinizadores
# --------------------------------------------------
dados_pol <- dados %>%
  filter(polinizador == "sim")

comm_poll <- criar_matriz_pa(dados_pol)
comm_poll <- criar_matriz_abund(dados_pol)
# --------------------------------------------------
# 6. Verificação rápida
# --------------------------------------------------
dim(comm_all)    # sites x morfoespécies (todos)
dim(comm_poll)   # sites x morfoespécies (polinizadores)

head(comm_all[, 1:5])
head(comm_poll[, 1:5])


# --------------------------------------------------
# 7. Extrair zona
# --------------------------------------------------
metadata_all <- tibble(site = rownames(comm_all)) %>%
  separate(site, into = c("ID_JC", "zona"), sep = "_", remove = FALSE)

metadata_poll <- tibble(site = rownames(comm_poll)) %>%
  separate(site, into = c("ID_JC", "zona"), sep = "_", remove = FALSE)

# --------------------------------------------------
# 8. Rodar NMDS
# --------------------------------------------------
##TODOS OS INSETOS
set.seed(123)  # reprodutibilidade

nmds_all <- metaMDS(
  comm_all,
  distance = "bray",
  k = 2,
  trymax = 100
)

nmds_all$stress #Interpretação do stress
#< 0.10 → excelente
#0.10–0.20 → aceitável
#0.20 → cuidado na interpretação

#POLINIZADORES

nmds_poll <- metaMDS(
  comm_poll,
  distance = "bray",
  k = 2,
  trymax = 100
)

nmds_poll$stress

#PREPARACAO PARA PLOTAGEM 
scores_all <- scores(nmds_all, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("site") %>%
  left_join(metadata_all, by = "site")

scores_poll <- scores(nmds_poll, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("site") %>%
  left_join(metadata_poll, by = "site")

##PLOTAGEM
##TODOS OS INSETOS
library(ggplot2)

ggplot(scores_all, aes(x = NMDS1, y = NMDS2, color = zona)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95) +
  theme_classic() +
  labs(
    title = "NMDS – All insects",
    subtitle = paste("Stress =", round(nmds_all$stress, 3)),
    color = "Zone"
  )

##POLINIZADORES
ggplot(scores_poll, aes(x = NMDS1, y = NMDS2, color = zona)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95) +
  theme_classic() +
  labs(
    title = "NMDS – Pollinators",
    subtitle = paste("Stress =", round(nmds_poll$stress, 3)),
    color = "Zone"
  )

#######PERMANOVA
# Metadados
# -----------------------------
# Metadados completos
metadata_all <- data.frame(site = rownames(comm_all)) %>%
  separate(site, into = c("ID_JC", "zona"), sep = "_", remove = FALSE)

# Metadados apenas dos sites com polinizadores
metadata_poll <- metadata_all %>%
  filter(site %in% rownames(comm_poll))

# -----------------------------
# Distâncias (Jaccard)
# -----------------------------
dist_all  <- vegdist(comm_all,  method = "jaccard", binary = TRUE)
dist_poll <- vegdist(comm_poll, method = "jaccard", binary = TRUE)

# -----------------------------
# PERMANOVA (estratificada por jardim)
# -----------------------------
perm_all <- adonis2(
  dist_all ~ zona,
  data = metadata_all,
  permutations = 999,
  strata = metadata$ID_JC
)

perm_poll <- adonis2(
  dist_poll ~ zona,
  data = metadata_poll,
  permutations = 999,
  strata = metadata_poll$ID_JC
)

# -----------------------------
# Homogeneidade de dispersão
# -----------------------------
disp_all  <- betadisper(dist_all,  metadata_all$zona)
disp_poll <- betadisper(dist_poll, metadata_poll$zona)

anova(disp_all)
anova(disp_poll)

# -----------------------------
# Resultados
# -----------------------------
perm_all
perm_poll
