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
library(ggpubr)

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

###########################################################
########### codigo para estimaar con menos datos ###########
data = SPM_data_1
time = c(2,3,4,5)

run_curve.eq07.time <- function ( data = NULL, time = c(2,3,4,5)){
  
    
  list.time <- unique (time)
   list.pot <- unique (data$pot)

     coefic.07.t <- bind_rows (lapply (list.pot, function (filt.pot){
       
 # filt.pot ="pot_3"
       
       pot.x <- data %>%
                dplyr::filter (pot == filt.pot )
       
   #list.time <- unique (time)
   
 dt.est <- bind_rows (lapply (list.time , function (filt.time) {
   
 # filt.time =  2
    
    y1 <- pot.x %>%
          dplyr::filter (time ==filt.time )%>%
          dplyr::select (ET.SPM)
    
    y2 <- pot.x %>%
          dplyr::filter  (time == 2* filt.time )%>%
          dplyr::select (ET.SPM)
    
    n <- 2* filt.time - filt.time
    
    
    Best <-  ( y1[[1]] ^ 2) / (2 * y1[[1]] - y2 [[1]])
    
    kest <-  (-1/n) * log (( y2 [[1]] / y1 [[1]]) - 1)
    
    
    tm.est <- 1/kest*log(2)
    
    ET.tm.est <- Best *(1- exp (-kest * tm.est))
    
    di <- str_c(filt.time,":", 2*filt.time)
    
    par.est <- data.frame (pot = filt.pot, dias =di, 
                           k.e=kest, B.e= Best, 
                           tm.e=tm.est, Btm.e=ET.tm.est)
    
    eq.07 <- nls (ET.SPM ~ B * (1 - exp (-k*time)),
                  start = list (k = kest [[1]], B = Best [[1]]),
                  trace = TRUE , 
                  data=  pot.x)
    
    coef.eq.07 <- rbind (coef (eq.07)) 
    
    k.eq.07 <- coef.eq.07[1,1] 
    
    B.eq.07 <- coef.eq.07[1,2]
    
    tm.1 <- 1/k.eq.07*log(2)
    
    ET.tm.1 <- B.eq.07 *(1- exp (- k.eq.07 * tm.1))
    
    coefic.07 <- data.frame ( par.est,  coef.eq.07, tm =tm.1[[1]], Btm = ET.tm.1 [[1]])
  
 }))
       
     }))    
       
     
     coefic.07.tx <- bind_rows (lapply (list.time , function (filt.time) {
       
       #filt.time =  2
       y1x <- data %>%
         dplyr::filter (time ==filt.time )%>%
         dplyr::select (ET.SPM)
       
       y1x <- mean (y1x$ET.SPM)
       
       y2x <- data %>%
         dplyr::filter (time ==2* filt.time )%>%
         dplyr::select(ET.SPM)
       
       y2x <- mean (y2x$ET.SPM)
       
       nx <- 2* filt.time - filt.time
       Bestx <-  (y1x [[1]] ^2) / (2* y1x [[1]] - y2x [[1]])
       kestx <-  (-1/nx) * log((y2x [[1 ]]/ y1x [[1]]) - 1)
       
       
       tm.estx <- 1/kestx*log(2)
       
       ET.tm.estx <- Bestx *(1- exp (-kestx * tm.estx))
       
       dix <- str_c(filt.time,":", 2*filt.time)
       
       par.estx <- data.frame (pot = "gral", dias =dix, 
                               k.e=kestx, B.e= Bestx, 
                               tm.e=tm.estx, Btm.e=ET.tm.estx)
       
       eq.07x <- nls (ET.SPM ~ B * (1 - exp (-k*time)),
                      start = list (k = kestx [[1]], B = Bestx [[1]]),
                      trace = TRUE , 
                      data=  data)
       
       coef.eq.07x <- rbind (coef (eq.07x)) 
       
       k.eq.07x <- coef.eq.07x[1,1] 
       
       B.eq.07x <- coef.eq.07x[1,2]
       
       tm.1x <- 1/k.eq.07x*log(2)
       
       ET.tm.1x <- B.eq.07x *(1- exp (- k.eq.07x * tm.1x))
       
       coefic.07x <- data.frame ( par.estx,  coef.eq.07x, tm =tm.1x[[1]], Btm = ET.tm.1x [[1]])
       
     })) 
     
param.time <- bind_rows( coefic.07.t,coefic.07.tx  )

param.time <- param.time %>%
              dplyr::mutate (k.e = round (k.e, 2)) %>%
              dplyr::mutate (B.e = round (B.e, 1)) %>%
              dplyr::mutate (tm.e = round (tm.e, 1)) %>%
              dplyr::mutate (Btm.e = round (Btm.e, 1)) %>%
              dplyr::mutate (k = round (k, 3)) %>%
              dplyr::mutate (B = round (B, 1)) %>%
              dplyr::mutate (tm = round (tm, 1)) %>%
              dplyr::mutate (Btm = round (Btm, 1))

return (param.time)    
       
       
     }


gauss.simondi <- run_curve.eq07.time (data =SPM_data_1, time = c(2,3,4,5))
#gral <- gauss.simondi  %>%
 #       dplyr::filter (pot == "gral")

my.reg <- function(x) 0 + 1*x 

ggscatter (gauss.simondi, x = "k", y = "k.e", color = "dias",
           palette = c("darkorange" , "navyblue", "red", "gray48"),
           size=3,
           #add="reg.line"
           ) + 
           stat_function( fun = my.reg, lwd=1) 

ggscatter (gauss.simondi,  x = "B", y = "B.e", color = "dias",
           palette = c("darkorange" , "navyblue", "red", "gray48"),
           size=3) + 
           stat_function( fun = my.reg, lwd=1) 






ggscatter (gral, x = "B", y = "B.e", color = "dias",
           palette = c("darkorange" , "navyblue","red" , "gray48"),
           size=3,
           #add="reg.line"
) + 
  stat_function( fun = my.reg, lwd=1)





#### funcionnes

myfun <- function(x) 232.8534 *(1- exp ( - 0.2077901  * x))
myfun.24 <- function(x) 547.9348  * (1- exp ( -0.06938826  * x))
myfun.36 <- function(x) 368.9362  * (1- exp ( -0.10787136  * x))
myfun.48 <- function(x) 233.2783  * (1- exp ( -0.21057258  * x))
myfun.510 <- function(x) 206.1552 * (1- exp ( -0.29167244  * x))


ggscatter (SPM_data_1, x = "time", y = "ET.SPM", 
           ylim=c(0,240),
           color = "gray38") + 
  stat_function(fun = myfun, lwd=1) +
  stat_function(fun = myfun.24, color = "darkorange" , lwd=1, lty=2) +
  stat_function(fun = myfun.36, color = "navyblue", lwd=1, lty=2) +
  stat_function(fun = myfun.48, color = "red", lwd=1.5, lty=2) +
  stat_function(fun = myfun.510, color = "gray48", lwd=1, lty=2)




# define function
gral <- gauss.simondi %>%
        dplyr::filter (pot== "gral")






curve  (232.8534 *(1- exp ( - 0.2077901  * x)), 
        from=0, to=10, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="gray38",lty= 1)













svg (filename="./Figures/Fig.2/fig.2B.svg", 
     width=7, 
     height=5, 
     pointsize=12)

plot (ET.SPM ~ time, pch=21, cex=1, 
      ylim=c(0,240),
      axes=T,
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
        col="gray38", bg="gray38",
        data=pot.3)

curve  ( 232.9 *(1- exp (- 0.208  * x)), 
         from=0, to=10, xlab="x", ylab="y",lwd=3,
         add = TRUE, col="red",lty= 1)

gauss.simondi 

curve  (gauss.simondi [1,3] *(1- exp ( - gauss.simondi [1,3]  * x)), 
         from=0, to=10, xlab="x", ylab="y",lwd=3,
         add = TRUE, col="gray38",lty= 1)



mean (gauss.simondi$ks_24)
