#######################################################################
# Codigo esquema paper soja Frontiers in Plants Science               ##
# Datos  Soja estres hidrico en Cámara de crecimiento                 ##
# Datos tomados por Esteban Casaretto                                 ##
# Sebastian Simondi - Gaston Quero                                    ##
# 13-07-2018                                                          ##
########################################################################

getwd ()
setwd ("R:/wtr")

# Paquetes 
library (lme4)
library (emmeans)
library ("car")
library ("nlmrt")
library ("easynls")
library (tidyverse)
library ("ggplot2")       
library ("lattice")
library ("latticeExtra")
library (multcompView)
library(reshape2)
library("dplyr")

########  se carga los datos 
# cargar datos 
SPM_data <- read.table ("./Data/rawdata/SPM_data.txt" ,
                           header = TRUE, sep = "",dec = ".",
                           na.strings = "NA" )
str (SPM_data)

summary (SPM_data)

##############
470 -453
mean.pss      <- mean (c(470,458,453,468))   # Esto es el peso seco del sustrato
sd.pss        <- sd   (c(470,458,453,468))

pm <- 41                           # Esto es el peso de la maceta

S <- mean.pss + pm                    # Este es el peso del sustrato mas 
                                      # el peso de la maceta

# el peso del SPM al tiempo cero
P0  <- SPM_data %>% 
       dplyr::filter (time==0) %>% 
       dplyr::select (w.SPM)

# genero la matriz de el agua del sustrato

SPM_data_1 <- cbind ( SPM_data, 
                      "A" = SPM_data$w.SPM - S,
                      "ET.SPM" = unique(P0$w.SPM) - SPM_data$w.SPM)

# filtro lps datos por maceta 
pot.1 <- SPM_data_1 %>% filter (pot== "pot_1")
pot.2 <- SPM_data_1 %>% filter (pot== "pot_2")
pot.3 <- SPM_data_1 %>% filter (pot== "pot_3")
pot.4 <- SPM_data_1 %>% filter (pot== "pot_4")
pot.5 <- SPM_data_1 %>% filter (pot== "pot_5")


####################  Figura 3A ##############

svg (filename="./Figures/Fig.3/fig_3A.svg", 
     width=7, 
     height=5, 
     pointsize=12)

plot (gw.SPM ~ w.SPM , pch=21, cex=1, 
      ylim=c(0,200),
      #xlim=c(-100,800),
      axes=F,
      las=2,
      xlab="", ylab="", 
      type="n",
      main="Gw vs weight",
      data=SPM_data_1)

points (gw.SPM ~ w.SPM, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=SPM_data_1)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("GW (mmol m-2 s-1 )",side=2,line=2.5)
axis (1,seq(550,750, 50))
mtext("W (g)",side=1,line=2.5)
box()

mod.0 <- lm (gw.SPM ~ w.SPM, data=SPM_data_1)
abline (mod.0 , lwd=1, col="black")

summary(mod.0)
coef.mod.0 <- coef (mod.0)

a <- coef.mod.0 [1]
b <- coef.mod.0 [2]

dev.off()

####################  Figura 3B  ##############

svg (filename="./Figures/Fig.3/fig_3B.svg", 
       width=7, 
       height=5, 
       pointsize=12)

plot (gw.SPM ~ time , pch=21, cex=1, 
        ylim=c(0,200),
        #xlim=c(-100,800),
        axes=F,
        las=2,
        xlab="", ylab="", 
        type="n",
        main="Gw vs time",
        data=SPM_data_1)
  
points (gw.SPM ~ time, pch=16, type="p", 
          lty=3, lwd=1.5, cex=1,
          col="black", bg="black",
          data=SPM_data_1)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext ("GW (mmol m-2 s-1 )",side=2,line=2.5)
axis (1,0:10)
mtext("time (d)",side=1,line=2.5)
box() 

AR.S <- as.numeric (unique (A0.eq.07))

curve  (b * B.eq.07 * exp (- k.eq.07  * x) + b * AR.S + a , 
        from=0, to=10, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="blue",lty= 1)

Gw.tm <- b * B.eq.07 * exp (- k.eq.07  * tm.1) + b * AR.S + a 

abline (h= Gw.tm  , col="red", lty=2)

tm.1 <- 1/k.eq.07*log(2)
abline (v=tm.1, col="red", lty=2)

GW.prim <-  2 * k.eq.07 * b * B.eq.07

abline (h= GW.prim , col="red", lty=2)

m.Best <- mean (Best [,1])
m.kest <- mean (kest [,1])
tm.1a <- 1/m.kest*log(2)

abline (v=tm.1a, col="red", lty=2) 
GW.prim.1 <-  2*m.kest *b*m.Best
abline (h= GW.prim.1 , col="blue", lty=2)

dev.off()

################### ACA termina el codigo para la figura 3A y 3B ###############################