library(Hmisc)
data <- read.csv("Data.csv")
data <- data[,c(1:10)]
table1 <- rcorr(as.matrix(data))
a <- table1$r

library(evir)

gev_fit <- gev(data, type = "GEV", location =~ long + lat + alt+popdensity+Temperature+
                  Precipitaion+Vapour.Pressure)
