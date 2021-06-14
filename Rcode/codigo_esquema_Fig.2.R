#######################################################################
# Codigo esquema paper modelo de consumo de agua en soja              ##
# Datos  Soja estres hidrico en Cámara de crecimiento                 ##
# Datos tomados por Esteban Casaretto                                 ##
# Sebastian Simondi - Gaston Quero                                    ##
# 14-06-2020                                                          ##
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

pm <- 41                                     # Esto es el peso de la maceta

S <- mean.pss + pm                    # Este es el peso del sustrato mas 
                                      # el peso de la maceta


#SPM_data$rem <- SPM_data$w.SPM - S
#


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

##### Figura 2A #### 

svg (filename="./Figures/Fig.2/fig_2A.svg", 
     width=7, 
     height=5, 
     pointsize=12)

plot (w.SPM ~ time, pch=21, cex=1, 
      ylim=c(500,760),
      xlim = c(0, 10),
      axes=F,
      las=2,
      xlab="", ylab="", 
      type="n",
      main="peso.SPM vs tiempo",
      data=SPM_data_1)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.1)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.2)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.3)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.4)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.5)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("peso (g)",side=2,line=2.5)
axis (1,0:50)
mtext("time (d)",side=1,line=2.5)

box()

##### Esta ecuacion la necesito para calcular AR

eq.05 <-  nls (SPM_data$w.SPM  ~ B * exp (-k*time) + A0 ,
                          start = list (k = 0.2, B = 235, A0=515),
                          trace = TRUE , 
                          data= SPM_data_1)

coef.eq.05 <- rbind (coef (eq.05))

k.eq.05  <- coef.eq.05  [1,1]
B.eq.05  <- coef.eq.05  [1,2]
A0.eq.05 <- coef.eq.05  [1,3]

# Aca calculo en AR 
AR <- A0.eq.05 - S

curve  (B.eq.05 *(exp (-k.eq.05 * x)) + AR + S, 
        from=0, to=50, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="red",lty= 1)

# el peso del SPM al tiempo cero
abline (h=unique(P0), col="red", lty=2) # esto es P(0)

# el peso del sustrato mas el peso de la maceta 
abline (h=S, col="blue", lty=2)               # esto es S
abline (h= A0.eq.05, col="darkgreen", lty=3)  # Esto estima el delta 

dev.off()



######### Figura 2B #### 

svg (filename="./Figures/Fig.2/fig_2B.svg", 
     width=7, 
     height=5, 
     pointsize=12)

plot (ET.SPM ~ time, pch=21, cex=1, 
      ylim=c(0,240),
      axes=F,
      las=2,
      xlab="", ylab="", 
      type="n",
      main="ET vs tiempo",
      data=SPM_data_1)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.1)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.2)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.3)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.4)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=0.7,
        col="black", bg="black",
        data=pot.5)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("ET (g)",side=2,line=2.5)
axis (1,0:10)
mtext("time (d)",side=1,line=2.5)
box()

#### 
x  <- SPM_data_1
j1 <- 4
W <- "w.SPM"

#ET.model <- function(x, j1, W) {
  W0  <- x %>% 
         dplyr::filter (time==0)  %>% 
         dplyr::select(all_of(W))
  
  Wn  <- x %>% 
         dplyr::filter (time==j1) %>% 
         dplyr::select(all_of(W))
  
  W2n <- x %>% 
        dplyr::filter (time==2*j1) %>% 
        dplyr::select(all_of(W))
  
  n <- 2*j1 - j1
  Best <-  (Wn - W0)^2 /( W2n - 2*Wn + W0 )
  kest <-  -(1/n) * log ( (W2n-Wn)/ (Wn-W0))
  par.est <- cbind (k.est=kest[1,1],B.est=Best[1,1] )
  
##### esta es la ecuacion 0.7 

eq.07 <- nls (ET.SPM ~ B * (1 - exp (-k*time)),
                 start = list (k = kest[1,1], B = Best[1,1]),
                 trace = TRUE , 
                 data= x)
  
coef.eq.07 <- rbind (coef (eq.07))
 # par.mod <- data.frame (pot=x[1,1],zz, par.est)
  #print(par.mod)

k.eq.07 <- coef.eq.07[1,1] 
  
B.eq.07 <- coef.eq.07[1,2]
  
########## usando los B y los k de la ecuacion [0.7]
curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
          from=0, to=10, xlab="x", ylab="y",lwd=3,
          add = TRUE, col="blue",lty= 1)
  
### AO.eq.07
#AR.eq.07 <- P0 - B.eq.07 - S
A0.eq.07 <- P0 - B.eq.07
#abline (h= unique(P0)- S, col="red", lty=2) 
  
tm.1 <- 1/k.eq.07*log(2)
abline (v=tm.1, col="red", lty=2)
  
  
  
ET.tm.1 <- B.eq.07 *(1- exp (- k.eq.07 * tm.1))
abline (h=ET.tm.1, col="red", lty=2) 
  
abline (h= B.eq.07, col="red", lty=2)
 

dev.off()

################ ACA termina el codigo para la figura 2 ##################