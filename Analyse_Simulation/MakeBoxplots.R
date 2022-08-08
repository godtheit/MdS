library(MdS2)
#setwd("H:/Documents/MdS/Analyse_Simulation/")
setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation/")
load("PlottingData.RData")

boxplot(ForPlots_MdS_Oracle$`Hamming-Distance`)

library(ggplot2)


MdS_Oracle <- ForPlots_MdS_Oracle[,1:3]
MdS_Oracle$Method <- "MdS-Oracle"
MdS_Aic <- ForPlots_MdS_Aic[,1:3]
MdS_Aic$Method <- "MdS-Aikaike-Kriterium"
Glasso <- ForPlots_glasso[,1:3]
Glasso$Method <- "Glasso-Oracle"
JGL <- ForPlots_JGL[41:60,1:3]
JGL$Method <- "JGL"

Allthedata <- rbind(MdS_Oracle, MdS_Aic, Glasso, JGL)

Alldata_noJGL <- rbind(MdS_Oracle, MdS_Aic, Glasso)

bp <- ggplot(Allthedata, aes(x = Allthedata$Method, y = Allthedata$`Hamming-Distance`, fill = Method)) +
  geom_boxplot()+
  labs(title="Hamming-Distanzen",x="Methode", y = "Hamming-Distanz")
bp + theme_classic()

bp <- ggplot(Alldata_noJGL, aes(x = Alldata_noJGL$Method, y = Alldata_noJGL$`Hamming-Distance`, fill = Method)) +
  geom_boxplot()+
  labs(title="Hamming-Distanzen, ohne JGL",x="Methode", y = "Hamming-Distanz")+
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#C77CFF"))
bp + theme_classic()



bp <- ggplot(Allthedata, aes(x = Allthedata$Method, y = Allthedata$Precision,fill = Method)) +
  geom_boxplot()+
  labs(title="Precision",x="Methode", y = "Precision = TP/(TP+FP)")
bp + theme_classic()


bp <- ggplot(Allthedata, aes(x = Allthedata$Method, y = Allthedata$Recall, fill = Method)) +
  geom_boxplot()+
  labs(title="Recall",x="Methode", y = "Recall = TP/(TP+FN)")
bp + theme_classic()
