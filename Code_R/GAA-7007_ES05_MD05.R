#
# ===========================================================================================
#            GAA-7007 (Analyse et modélisation des agroécosystèmes)
#                             [Evaluation Sommative 5]
#
#                         Analyse des séries temporelles
#
#                             Par : Marc Duchemin
#                               (26 avril 2020)
# ================================================================================
# ================================================================================
#
##Chargement des librairies
#
library(csv)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(caret)
library(forecast)
library(fpp2)
library(ggpubr)
library(e1071)
#
#======================================================================
#
##Chargement des données et création d'une série temporelle
#-----
#a)Établir un lien absolu dans un Working Directory
#setwd("C:\\Users\\Utilisateur\\Documents\\#_ULAVAL\\GAA-7007_ES_5\\Code_R")
# ou bien avec ...
#b)Établir un lien relatif avec : Session\Set Working Directory\To Source File Location
# OK pour l'option b) ...
#-----
#Création d'un dataframe (tableau initial)
HawaiData <- read.csv("Data/hawai_MD.csv")
colnames(HawaiData) <- c("Date", "CO2")
#Transformation en Time-Series
HawaiData_ts <- ts(HawaiData, frequency = 12, start = c(1958,3))
#
#======================================================================
#
##Exploration des données
#
#Statistiques descriptives du CO₂ atmosphérique
summary(HawaiData_ts)
#
##Histogramme du CO₂ atmosphériqe
q1 <- ggplot(HawaiData, aes(CO2)) +  
  geom_histogram(binwidth = 1, color = "brown", fill = "green")
countminC = min(ggplot_build(q1)$data[[1]]$count)
countmaxC = max(ggplot_build(q1)$data[[1]]$count)
ggplot(HawaiData, aes(CO2)) +  
  geom_histogram(binwidth = 1, color = "brown", fill = "green") +
  geom_vline(aes(xintercept = mean(HawaiData$CO2)),
             color ="red", linetype ="dashed", size = 1) +
  geom_text(x = (mean(HawaiData$CO2)),
            y = ((countmaxC-countminC)/2),
            aes(fontface=2),label = "concentration moyenne", 
            color = "black", size = 4,srt = 90) +
  scale_x_continuous(breaks = seq(trunc(min(HawaiData$CO2)), ceiling(max(HawaiData$CO2)), by = 5)) +
  scale_y_continuous(breaks = seq(0, countmaxC, by = 2)) +
  labs(title = "Histogramme du CO₂ atmosphérique au Mauna Loa Observatory (Hawaï) entre 1958 et 2001",
       x = "Concentration en CO₂ atmosphérique (ppm)",
       y = "Fréquence")
ggsave("Images/Figure_HistoCO2.png", width = 9, height = 7, dpi = 300)
#
#Vérification de la normalité du CO₂ atmosphérique
#a)Q-Q Plot
ggqqplot(HawaiData$CO2,
         main = "Graphique Q-Q plot du CO₂ atmosphérique",
         xlab = "Théorique",
         ylab = "Concentration en CO₂ atmosphérique (ppm)")
ggsave("Images/Figure_QQplot.png", width = 9, height = 7, dpi = 300)
#b)test de normalité de Shapiro-Wilk (si p-value > 0.05 : normalité)
shapiro.test(HawaiData$CO2)
#
#Graphique de la variation temporelle du  CO₂ 
HawaiData %>%
  ggplot(aes(x = Date, y = CO2)) +
  geom_line(color = "darkgreen", linetype = 1, size = 1) +
  scale_x_continuous(breaks = seq(trunc(min(HawaiData)), ceiling(max(HawaiData)), by = 4)) +
  scale_y_continuous(breaks = seq(trunc(min(HawaiData)), ceiling(max(HawaiData)), by = 10)) +
  labs(title = "Variation des concentrations de CO₂ au Mauna Loa Observatory (Hawaï) entre 1958 et 2001",
       x = "Temps (année)",
       y = "Concentration en CO₂ atmosphérique (ppm)")
ggsave("Images/Figure_VarCO2a.png", width = 9, height = 7, dpi = 300)
#
#Extraction d'un sous-ensemble de dates (pour faire ressortir l'effet cyclique)
HawaiData_part <- HawaiData[which(HawaiData$Date > 1976 & HawaiData$Date < 1984),]
HawaiData_part %>%
  ggplot(aes(x = Date, y = CO2)) +
  geom_line(color = "darkgreen", linetype = 1, size = 1) +
  scale_y_continuous(breaks = seq(trunc(min(HawaiData[2])), ceiling(max(HawaiData[2])), by = 2)) +
  labs(title = "Variation des concentrations de CO₂ au Mauna Loa Observatory (Hawaï) entre 1976 et 1984",
       x = "Temps (année)",
       y = "Concentration en CO₂ atmosphérique (ppm)")
ggsave("Images/Figure_VarCO2b.png", width = 9, height = 7, dpi = 300)
#
#vérification de l'autocorrélation du CO₂
HD_ggAcf <- ggAcf(HawaiData_ts[,2], ci = 0.95, type = c("correlation"), plot = TRUE) +
  labs(title = "CO₂ : autocorrélation VS retardement",
       x = "Retardement (Lag)",
       y = "Coefficient d'autocorrélation (ACF)")
HD_ggAcf
ggsave("Images/Figure_AutoCorCO2.png", width = 9, height = 7, dpi = 300)
#
#Le test de Ljung-Box permet de vérifier si la série temporelle entière
#peut être différenciée d’un bruit blanc.
#Vérification d'un bruit blanc et d'une structure dans les données temporelles
#si p-value <0,05 : faible possibilité de bruit blanc et présence d'une structure
Box.test(HawaiData_ts[,2], lag = max(HD_ggAcf$data$lag), type = "Ljung-Box")
#note: si p<<<<0,05 : probabilité presque nulle d'un bruit blanc.
#
#======================================================================
#
##Séparation des données en deux jeux : entraînement (70%) et test (30%)
#
HawaiData_tsX <- HawaiData_ts[,2]
#jeu d'entraînement
HawaiTrain <- window(HawaiData_tsX, start = 1958.16667, end = 1988.917)
#jeu de test
HawaiTest <- window(HawaiData_tsX, start = 1989.00, end = 2001.91667)
#
#Statistiques descriptives des données Train et Test
summary(HawaiTrain) #données d'entraînement
summary(HawaiTest) #données de test
#
##BOXPLOT des données Train et Test (avec V-notch)
#A)préparation des données
#copies des jeux 
HawaiTrainX <- HawaiData[which(HawaiData$Date >= min(HawaiData$Date) 
                               & HawaiData$Date < 1989.0),]
HawaiTestX <- HawaiData[which(HawaiData$Date >= 1989.0 
                               & HawaiData$Date <= max(HawaiData$Date)),]
#codification (factor) des données Train et Test
#-ajouter une variable Code-facteur "Entraînement" aux données Train
HawaiTrainX$Code <- (HawaiTrainX / HawaiTrainX)
HawaiTrainX$Code <- as.factor(HawaiTrainX$Code)
HawaiTrainX$Code = "Entraînement"
HawaiTrainX$Code <- as.factor(HawaiTrainX$Code)
#-ajouter une variable Code-facteur "Test" aux données Test
HawaiTestX$Code <- (HawaiTestX[2] / HawaiTestX[2])
HawaiTestX$Code <- as.factor(HawaiTestX$Code)
HawaiTestX$Code = "Test"
HawaiTestX$Code <- as.factor(HawaiTestX$Code)
#-réunion des deux tables Train et Test avec Codification (factor)
HawaiTrainTest <- bind_rows(HawaiTrainX, HawaiTestX)
#B)Création des boxplots
ggplot(HawaiTrainTest, aes(Code, CO2, fill = Code, colour = Code)) +
  geom_boxplot(show.legend = FALSE,
               color = "brown", alpha=0.25,
               notch = TRUE, notchwidth = 0.5,
               coef = 1.5, outlier.colour = "red", 
               outlier.fill = "red", outlier.shape = 20,
               outlier.size = 6, outlier.alpha = 0.5, na.rm = TRUE) +
  geom_jitter(show.legend = FALSE, fill = "black") +
  stat_summary(fun = mean, geom = "point",
               shape = 22, fill = "darkgreen", 
               color = "white", size = 5, 
               alpha = 0.5, na.rm = TRUE) +
  theme(axis.text.x = element_text(size=12)) +
  scale_y_continuous(breaks = seq(trunc(min(HawaiData[2])), ceiling(max(HawaiData[2])), by = 10)) +
  labs(title = "Diagrammes en boîte à moustaches du CO₂ atmosphérique",
       x = "Groupes de données",
       y = "Concentration en CO₂ atmosphérique (ppm)")
ggsave("Images/Figure_BoxPlotJeux.png", width = 9, height = 7, dpi = 300)
#
#======================================================================
#
##Création et projection d'un modèle ETS (Erreur, Tendance, Saison)
#
#a)modélisation ETS
CO2_ets <- HawaiTrain %>% ets(model = "ZZZ", damped = NULL)
CO2_ets %>% autoplot()
ggsave("Images/Figure_ETSComposants.png", width = 9, height = 7, dpi = 300)
#b)projection (forecast)
CO2_forecast <- CO2_ets %>% forecast(h = length(HawaiTest))
CO2_forecast %>% autoplot()
ggsave("Images/Figure_ETSPrédiction.png", width = 9, height = 7, dpi = 300)
#
summary(CO2_ets)
#
#======================================================================
#
##Effectuer une analyse des résidus
#-test de Ljung-Box (si p-value > 0.05 : présence de bruit blanc)
CO2_ets %>% checkresiduals()
ggsave("Images/Figure_Residus.png", width = 9, height = 7, dpi = 300)
#sauvegarder avec [Export\Save as Image] de RStudio:""Images/Figure_Residusx.png" [W:900,H:750]
#
#-test de Shapiro-Wilk (si p-value > 0.05 : normalité)
shapiro.test(residuals(CO2_ets))
#
#======================================================================
#
##Sauvegarder un fichier(.RData) [si désiré]
#save.image("GAA-7007_ES05_MD05.RData")
#=======================================================================
#