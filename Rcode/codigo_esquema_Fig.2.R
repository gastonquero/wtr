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

############## datos de los pesos de las macetas y el sustrato ########
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

#data= SPM_data_1 

##### ecuacion 05 ##########
run_curve.eq05 <- function ( data = NULL){
   
  list.pot <- unique (data$pot)
  
  coefic <- bind_rows (lapply (list.pot, function (filt.pot){
    
    #filt.pot ="pot_1"
    
    pot.x <- data %>%
             dplyr::filter (pot == filt.pot )
    
    eq.05 <-  nls (w.SPM  ~ B * exp (-k*time) + A0 ,
                   start = list (k = 0.2, B = 235, A0=515),
                   trace = TRUE , 
                   data= pot.x)
    
    coefi <- rbind (coef (eq.05))
    
    A0.eq.05 <- coefi  [1,3]
    
    AR <- A0.eq.05 - S
    
    coefic <- data.frame (pot=filt.pot,  coefi, AR = AR [[1]])
  
  }))
  
  eq.05x <-  nls (w.SPM  ~ B * exp (-k*time) + A0 ,
                  start = list (k = 0.2, B = 235, A0=515),
                  trace = TRUE , 
                  data= data)
  
  coefix <- rbind (coef (eq.05x))
  
  A0.eq.05x <- coefix  [1,3]
  
  ARx <- A0.eq.05x - S
  
  coefix <- data.frame (pot="gral",  coefix, AR = ARx [[1]])
  
  par.model <- bind_rows (coefic, coefix)
  
  par.model <- par.model %>%
               dplyr::mutate ( k  = round (k, 3)) %>%
               dplyr::mutate ( B  = round (B, 1)) %>%
               dplyr::mutate ( A0 = round (A0, 1)) %>%
               dplyr::mutate ( AR = round (AR, 1)) 
  
  return( par.model)

}

parSPM.eq05 <- run_curve.eq05 (data= SPM_data_1 )
parSPM.eq05

################ figura 2A ############
svg (filename="./Figures/Fig.2/fig.2A.svg", 
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
      data=SPM_data_1 )

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.1)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.2)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.4)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.5)

points (w.SPM ~ time, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.3)

curve  (224.6  *(exp (-0.263 * x)) + 29.6 + 503.25, 
        from=0, to=50, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="black",lty= 1)

curve  (234.7 *(exp (-0.225 * x)) +19.0 + 503.25, 
        from=0, to=50, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="red",lty= 1)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("peso (g)",side=2,line=2.5)
axis (1,0:50)
mtext("time (d)",side=1,line=2.5)

box()

# el peso del SPM al tiempo cero
abline (h=unique(P0), col="red", lty=2) # esto es P(0)

# el peso del sustrato mas el peso de la maceta 
abline (h=S, col="blue", lty=2)               # esto es S
abline (h= 532.8, col="black", lty=2)  # Esto estima el delta 
abline (h= 522.3, col="red", lty=2)

dev.off()

########## ecuacion 07 ##########
run_curve.eq07 <- function ( data = NULL){
  
  list.pot <- unique (data$pot)
  
  coefic.07 <- bind_rows (lapply (list.pot, function (filt.pot){
    
    #filt.pot ="pot_3"
    x <- data %>%
             dplyr::filter (pot == filt.pot )
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
    
    #### esta es la ecuacion 0.7 
    
    eq.07 <- nls (ET.SPM ~ B * (1 - exp (-k*time)),
                  start = list (k = kest[1,1], B = Best[1,1]),
                  trace = TRUE , 
                  data= x)
    
    coef.eq.07 <- rbind (coef (eq.07)) 
    
    k.eq.07 <- coef.eq.07[1,1] 
    
    B.eq.07 <- coef.eq.07[1,2]
    
    tm.1 <- 1/k.eq.07*log(2)

    ET.tm.1 <- B.eq.07 *(1- exp (- k.eq.07 * tm.1))
    
    coefic.07 <- data.frame (pot=filt.pot,  coef.eq.07, tm =tm.1[[1]], Btm = ET.tm.1 [[1]])

    
  }))
  
  eq.07x <- nls (ET.SPM ~ B * (1 - exp (-k*time)),
                start = list (k = kest[1,1], B = Best[1,1]),
                trace = TRUE , 
                data= data)
  
  coef.eq.07x <- rbind (coef (eq.07x)) 
  
  k.eq.07x <- coef.eq.07x[1,1] 
  
  B.eq.07x <- coef.eq.07x[1,2]
  
  tm.1x <- 1/k.eq.07x*log(2)
  
  ET.tm.1x <- B.eq.07x *(1- exp (- k.eq.07x * tm.1x))
  
  coefix.07 <- data.frame (pot="gral",  coef.eq.07x, tm =tm.1x[[1]], Btm = ET.tm.1x [[1]])
  
  par.model.07 <- bind_rows (coefic.07, coefix.07)
  
  
  par.model.07 <- par.model.07 %>%
                  dplyr::mutate ( k  = round (k, 3)) %>%
                  dplyr::mutate ( B  = round (B, 1)) %>%
                  dplyr::mutate ( tm = round (tm, 2)) %>%
                  dplyr::mutate ( Btm = round (Btm, 1)) 
  
  return (par.model.07)
  
}

parSPM.eq07 <- run_curve.eq07 (data= SPM_data_1 )
parSPM.eq07




#### ##### Figura 2B #### 
svg (filename="./Figures/Fig.2/fig.2B.svg", 
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
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.1)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.2)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.4)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="gray38", bg="gray38",
        data=pot.5)

points (ET.SPM ~ time , pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.3)


########## usando los B y los k de la ecuacion [0.7]

curve  ( 220.8 *(1- exp (- 0.2457  * x)), 
        from=0, to=10, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="black",lty= 1)

curve  ( 232.9 *(1- exp (- 0.208  * x)), 
         from=0, to=10, xlab="x", ylab="y",lwd=3,
         add = TRUE, col="red",lty= 1)

abline (v=2.83, col="black", lty=2)
abline (v=3.34, col="red", lty=2)

abline (h= 110.4, col="black", lty=2) 
abline (h= 116.4, col="red", lty=2)


abline (h= 220.8, col="black", lty=2) 
abline (h= 232.9, col="red", lty=2)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("ET (g)",side=2,line=2.5)
axis (1,0:10)
mtext("time (d)",side=1,line=2.5)
box()

dev.off()



############ ACA termina el codigo de la figura 2B ############
