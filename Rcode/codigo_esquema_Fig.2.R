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

################ figura 1B 

### ver los datos 

##### Figura 2A #### 
svg (filename="./Figures/Fig.1/Pt.svg", 
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



#### ##### Figura 1C #### 

svg (filename="./Figures/Fig.1/Et.svg", 
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
  W0  <- x %>% filter (time==0)  %>% select(W)
  Wn  <- x %>% filter (time==j1) %>% select(W)
  W2n <- x %>% filter (time==2*j1) %>% select(W)
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


####################  Figura 1D ##############

svg (filename="./Figures/Fig.1/Gw_vs_Wg.svg", 
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

####################  Figura 1E  ##############

svg (filename="./Figures/Fig.1/Gw_vs_time.svg", 
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


#curve  (b * B.eq.07 * exp (- k.eq.07  * x) + b * AR.S + a , 
        #from=0, to=10, xlab="x", ylab="y",lwd=3,
        #add = TRUE, col="blue",lty= 1)

m.Best <- mean (Best [,1])
m.kest <- mean (kest [,1])
tm.1a <- 1/m.kest*log(2)

abline (v=tm.1a, col="red", lty=2) 
GW.prim.1 <-  2*m.kest *b*m.Best
abline (h= GW.prim.1 , col="blue", lty=2)



curve  (b.2 * m.Best * exp (- m.kest  * x) + b.2 * AR.S + a.2 , 
        from=0, to=10, xlab="x", ylab="y",lwd=3,
        add = TRUE, col="red",lty= 1)



#dev.off()
###########################################################
########### codigo para estimaar con menos datos ###########

ee.1 <- SPM_data_1 %>% filter (pot== "pot_1")
ee.2 <- SPM_data_1 %>% filter (pot== "pot_2")
ee.3 <- SPM_data_1 %>% filter (pot== "pot_3")
ee.4 <- SPM_data_1 %>% filter (pot== "pot_4")
ee.5 <- SPM_data_1 %>% filter (pot== "pot_5")
head(ee.5)

et.model.2 <- function(x, y) {
  y1 <- filter (x, time==y)    %>% select(ET.SPM)
  y2 <- filter (x, time==2*y)  %>% select(ET.SPM)
  n <- 2*y - y
  Aest <-  (y1^2) / (2*y1 - y2)
  kest <-  (-1/n) * log((y2/y1) - 1)
  par.est <- cbind (k.est=kest[1,1],A.est=Aest[1,1] )
  ET.par <- nls (ET.SPM ~ A * (1 - exp (-k*time)),
                 start = list (k = kest[1,1], A = Aest[1,1]),
                 trace = TRUE , 
                 data= x)
  zz <- rbind (coef (ET.par))
  par.mod <- data.frame (pot=x[1,1],zz, par.est)
  print(par.mod)
}




rp1.2 <- et.model.2 (ee.1, y=2)
rp2.2 <- et.model.2 (ee.2, y=2)
rp3.2 <- et.model.2 (ee.3, y=2)
rp4.2 <- et.model.2 (ee.4, y=2)
rp5.2 <- et.model.2 (ee.5, y=2)

rp1.3 <- et.model.2 (ee.1, y=3)
rp2.3 <- et.model.2 (ee.2, y=3)
rp3.3 <- et.model.2 (ee.3, y=3)
rp4.3 <- et.model.2 (ee.4, y=3)
rp5.3 <- et.model.2 (ee.5, y=3)

rp1.4 <- et.model.2 (ee.1, y=4)
rp2.4 <- et.model.2 (ee.2, y=4)
rp3.4 <- et.model.2 (ee.3, y=4)
rp4.4 <- et.model.2 (ee.4, y=4)
rp5.4 <- et.model.2 (ee.5, y=4)

rp1.5 <- et.model.2 (ee.1, y=5)
rp2.5 <- et.model.2 (ee.2, y=5)
rp3.5 <- et.model.2 (ee.3, y=5)
rp4.5 <- et.model.2 (ee.4, y=5)
rp5.5 <- et.model.2 (ee.5, y=5)

ee.nls.par.24 <- rbind (rp1.2 ,rp2.2,rp3.2,rp4.2,rp5.2)
ee.nls.par.36 <- rbind (rp1.3 ,rp2.3,rp3.3,rp4.3,rp5.3)
ee.nls.par.48 <- rbind (rp1.4 ,rp2.4,rp3.4,rp4.4,rp5.4)
ee.nls.par.510 <- rbind (rp1.5 ,rp2.5,rp3.5,rp4.5,rp5.5)

ee.nls.par.510 <- rbind (rp1.5 ,rp2.5,rp3.5,rp4.5,rp5.5)

gauss.simondi <- ee.nls.par.24 %>% 
  rename (B =A) %>%
  rename (ks_24 = k.est) %>%
  rename (Bs_24 = A.est) %>%
  mutate (ks_36 = ee.nls.par.36$k.est)%>%
  mutate (Bs_36 = ee.nls.par.36$A.est)%>%
  mutate (ks_48 = ee.nls.par.48$k.est)%>%
  mutate (Bs_48 = ee.nls.par.48$A.est)%>%
  mutate (ks_510 = ee.nls.par.510$k.est)%>%
  mutate (Bs_510 = ee.nls.par.510$A.est)

# Hago la matriz de correlaciones 
cormat <- round(cor(gauss.simondi[,-1]),2)
melted_cormat <- melt(cormat)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

#Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

###################### minimos para Gw ############################
View (SPM_data_1)
gw0 <- SPM_data_1 %>% filter  (time == 0)
gw2 <- SPM_data_1 %>% filter  (time == 2)
gw3 <- SPM_data_1 %>% filter  (time == 3)
gw4 <- SPM_data_1 %>% filter  (time == 4)
gw5 <- SPM_data_1 %>% filter  (time == 5)
gw6 <- SPM_data_1 %>% filter  (time == 6)
gw8 <- SPM_data_1 %>% filter  (time == 8)
gw10 <- SPM_data_1 %>% filter (time == 10)

##### selecciono tiempo y pot 
ee.1    <- SPM_data_1 %>% filter (pot== "pot_1")
ee.1.0  <- gw0  %>% filter (pot== "pot_1")
ee.1.2  <- gw2  %>% filter (pot== "pot_1")
ee.1.4  <- gw4  %>% filter (pot== "pot_1")
ee.1.8  <- gw8  %>% filter (pot== "pot_1")

#### modelo completo #####
mod.ee.1 <- lm (gw.SPM ~ w.SPM, data=ee.1)
summary(mod.ee.1)
coef.mod.ee.1 <- coef (mod.ee.1)

a.ee.1 <- coef.mod.ee.1 [1]
b.ee.1  <- coef.mod.ee.1 [2]

######## modelo con 0, 2, 4, 8 ########
ee.1.248     <- rbind (ee.1.0, ee.1.2,ee.1.4, ee.1.8 )
mod.ee.1.248 <- lm (gw.SPM ~ w.SPM, data=ee.1.248)
summary(mod.ee.1.248)
coef.mod.ee.1.248 <- coef (mod.ee.1.248)

a.ee.1.248 <- coef.mod.ee.1.248 [1]
b.ee.1.248  <- coef.mod.ee.1.248 [2]

coef.ee.1 <- cbind (pot=1, 
                    a.lm = a.ee.1,
                    b.lm = b.ee.1,
                    a.248 = a.ee.1.248,
                    b.248 = b.ee.1.248)
#### 
##### selecciono tiempo y pot 
ee.2    <- SPM_data_1 %>% filter (pot== "pot_2")
ee.2.0  <- gw0  %>% filter (pot== "pot_2")
ee.2.2  <- gw2  %>% filter (pot== "pot_2")
ee.2.4  <- gw4  %>% filter (pot== "pot_2")
ee.2.8  <- gw8  %>% filter (pot== "pot_2")

#### modelo completo #####
mod.ee.2 <- lm (gw.SPM ~ w.SPM, data=ee.2)
summary(mod.ee.2)
coef.mod.ee.2 <- coef (mod.ee.2)

a.ee.2 <- coef.mod.ee.2 [1]
b.ee.2  <- coef.mod.ee.2 [2]

######## modelo con 0, 2, 4, 8 ########
ee.2.248     <- rbind (ee.2.0, ee.2.2,ee.2.4, ee.2.8 )
mod.ee.2.248 <- lm (gw.SPM ~ w.SPM, data=ee.2.248)
summary(mod.ee.2.248)
coef.mod.ee.2.248 <- coef (mod.ee.2.248)

a.ee.2.248 <- coef.mod.ee.2.248 [1]
b.ee.2.248  <- coef.mod.ee.2.248 [2]

coef.ee.2 <- cbind (pot=2, 
                    a.lm = a.ee.2,
                    b.lm = b.ee.2,
                    a.248 = a.ee.2.248,
                    b.248 = b.ee.2.248)

######## ee.3 
ee.3    <- SPM_data_1 %>% filter (pot== "pot_3")
ee.3.0  <- gw0  %>% filter (pot== "pot_3")
ee.3.2  <- gw2  %>% filter (pot== "pot_3")
ee.3.4  <- gw4  %>% filter (pot== "pot_3")
ee.3.8  <- gw8  %>% filter (pot== "pot_3")

#### modelo completo #####
mod.ee.3 <- lm (gw.SPM ~ w.SPM, data=ee.3)
summary(mod.ee.3)
coef.mod.ee.3 <- coef (mod.ee.3)

a.ee.3 <- coef.mod.ee.3 [1]
b.ee.3  <- coef.mod.ee.3 [2]

######## modelo con 0, 2, 4, 8 ########
ee.3.248     <- rbind (ee.3.0, ee.3.2,ee.3.4, ee.3.8 )
mod.ee.3.248 <- lm (gw.SPM ~ w.SPM, data=ee.3.248)
summary(mod.ee.3.248)
coef.mod.ee.3.248 <- coef (mod.ee.3.248)

a.ee.3.248 <- coef.mod.ee.3.248 [1]
b.ee.3.248  <- coef.mod.ee.3.248 [2]

coef.ee.3 <- cbind (pot=3, 
                    a.lm = a.ee.3,
                    b.lm = b.ee.3,
                    a.248 = a.ee.3.248,
                    b.248 = b.ee.3.248)

######## ee.4 
ee.4    <- SPM_data_1 %>% filter (pot== "pot_4")
ee.4.0  <- gw0  %>% filter (pot== "pot_4")
ee.4.2  <- gw2  %>% filter (pot== "pot_4")
ee.4.4  <- gw4  %>% filter (pot== "pot_4")
ee.4.8  <- gw8  %>% filter (pot== "pot_4")

#### modelo completo #####
mod.ee.4 <- lm (gw.SPM ~ w.SPM, data=ee.4)
summary(mod.ee.4)
coef.mod.ee.4 <- coef (mod.ee.4)

a.ee.4 <- coef.mod.ee.4 [1]
b.ee.4  <- coef.mod.ee.4 [2]

######## modelo con 0, 2, 4, 8 ########
ee.4.248     <- rbind (ee.4.0, ee.4.2,ee.4.4, ee.4.8 )
mod.ee.4.248 <- lm (gw.SPM ~ w.SPM, data=ee.4.248)
summary(mod.ee.4.248)
coef.mod.ee.4.248 <- coef (mod.ee.4.248)

a.ee.4.248 <- coef.mod.ee.4.248 [1]
b.ee.4.248  <- coef.mod.ee.4.248 [2]

coef.ee.4 <- cbind (pot=4, 
                    a.lm = a.ee.4,
                    b.lm = b.ee.4,
                    a.248 = a.ee.4.248,
                    b.248 = b.ee.4.248)
###########
######## ee.5 
ee.5    <- SPM_data_1 %>% filter (pot== "pot_5")
ee.5.0  <- gw0  %>% filter (pot== "pot_5")
ee.5.2  <- gw2  %>% filter (pot== "pot_5")
ee.5.4  <- gw4  %>% filter (pot== "pot_5")
ee.5.8  <- gw8  %>% filter (pot== "pot_5")

#### modelo completo #####
mod.ee.5 <- lm (gw.SPM ~ w.SPM, data=ee.5)
summary(mod.ee.5)
coef.mod.ee.5 <- coef (mod.ee.5)

a.ee.5 <- coef.mod.ee.5 [1]
b.ee.5  <- coef.mod.ee.5 [2]

######## modelo con 0, 2, 4, 8 ########
ee.5.248     <- rbind (ee.5.0, ee.5.2,ee.5.4, ee.5.8 )
mod.ee.5.248 <- lm (gw.SPM ~ w.SPM, data=ee.5.248)
summary(mod.ee.5.248)
coef.mod.ee.5.248 <- coef (mod.ee.5.248)

a.ee.5.248 <- coef.mod.ee.5.248 [1]
b.ee.5.248  <- coef.mod.ee.5.248 [2]

coef.ee.5 <- cbind (pot=5, 
                    a.lm = a.ee.5,
                    b.lm = b.ee.5,
                    a.248 = a.ee.5.248,
                    b.248 = b.ee.5.248)

coef.GW <- rbind (coef.ee.1,coef.ee.2,coef.ee.3,coef.ee.4,coef.ee.5)

dim (coef.GW)

coef.GW [,-1]
          
cormat.Gw <- round(cor(coef.GW [,-1]),2)
melted_cormat.Gw <- melt(cormat.Gw)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat.Gw){
  cormat.Gw[upper.tri(cormat.Gw)] <- NA
  return(cormat.Gw)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat.Gw){
  cormat.Gw[lower.tri(cormat.Gw)]<- NA
  return(cormat.Gw)
}

upper_tri <- get_upper_tri(cormat.Gw)
upper_tri

#Melt the correlation matrix
melted_cormat.Gw <- melt (upper_tri, na.rm = TRUE)
# Heatmap

reorder_cormat.Gw <- function(cormat.Gw){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat.Gw)/2)
  hc <- hclust(dd)
  cormat.Gw <-cormat.Gw[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat.Gw <- reorder_cormat.Gw(cormat.Gw)
upper_tri <- get_upper_tri(cormat.Gw)
# Melt the correlation matrix
melted_cormat.Gw <- melt(upper_tri, na.rm = TRUE)
ggheatmap.Gw <- ggplot(melted_cormat.Gw, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap.Gw)
ggheatmap.Gw + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



ee.1.gw.todos <- SPM_data_1 %>% filter (pot== "pot_1")

mod.0 <- lm (gw.SPM ~ w.SPM, data=SPM_data_1)
abline (mod.0 , lwd=1, col="black")


summary(mod.0)
coef.mod.0 <- coef (mod.0)

a <- coef.mod.0 [1]
b <- coef.mod.0 [2]







gw0 <- SPM_data_1 %>%
       filter (time == 0)

#gw1 <- SPM_data_1 %>%
      #filter (time == 1)

gw2 <- SPM_data_1 %>%
  filter (time == 2)


gw4 <- SPM_data_1 %>%
       filter (time == 4)

gw8 <- SPM_data_1 %>%
       filter (time == 8)


############
gw0248 <- rbind (gw0,gw2, gw4, gw8)

mod.0248 <- lm (gw.SPM ~ w.SPM, data=gw0248)
#abline (mod.0248, lwd=2, col="red", tty =3)


summary (mod.0248)
coef.mod.0248 <- coef (mod.0248)

a.2 <- coef.mod.0248 [1]
b.2 <- coef.mod.0248 [2]




gw048 <- rbind (gw0, gw4, gw8)

mod.048 <- lm (gw.SPM ~ w.SPM, data=gw048)
abline (mod.048, lwd=1.3, col="blue")


summary (mod.048)
coef.mod.048 <- coef (mod.048)

a.1 <- coef.mod.048 [1]
b.1 <- coef.mod.048 [2]
a - a.1 

b - b.1



gw0148 <- rbind (gw0,gw1, gw4, gw8)

mod.0148 <- lm (gw.SPM ~ w.SPM, data=gw0148)
abline (mod.0148, lwd=2, col="red")


summary (mod.0148)
coef.mod.0148 <- coef (mod.0148)

a.2 <- coef.mod.0148 [1]
b.2 <- coef.mod.0148 [2]
a - a.2

b - b.2

##############

gw148 <- rbind (gw1, gw4, gw8)

mod.148 <- lm (gw.SPM ~ w.SPM, data=gw148)
abline (mod.148, lwd=2, col="green")


summary (mod.048)
coef.mod.048 <- coef (mod.048)

a.1 <- coef.mod.048 [1]
b.1 <- coef.mod.048 [2]
a - a.1 

b - b.1


############
gw0248 <- rbind (gw0,gw2, gw4, gw8)

mod.0248 <- lm (gw.SPM ~ w.SPM, data=gw0248)
abline (mod.0248, lwd=2, col="red", tty =3)


summary (mod.0248)
coef.mod.0248 <- coef (mod.0248)

a.2 <- coef.mod.0248 [1]
b.2 <- coef.mod.0248 [2]
a - a.2

b - b.2

gw0148 <- rbind (gw0,gw1, gw4, gw8)

mod.0148 <- lm (gw.SPM ~ w.SPM, data=gw0148)
abline (mod.0148, lwd=2, col="red")


summary (mod.0148)
coef.mod.0148 <- coef (mod.0148)

a.2 <- coef.mod.0148 [1]
b.2 <- coef.mod.0148 [2]
a - a.2

b - b.2






#### figura Gw vs ET.diaria

plot (ET.diaria ~ gw.SPM , pch=21, cex=1, 
      #ylim=c(0,200),
      #xlim=c(-100,800),
      axes=T,
      las=2,
      #xlab="", 
      #ylab="", 
      type="n",
      main="ET vs Gw",
      data=SPM_data_1)

points (ET.diaria ~ gw.SPM, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.1)



points (ET.diaria ~ gw.SPM, pch=17, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.2)

points (ET.diaria ~ gw.SPM, pch=19, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.3)


points (ET.diaria ~ gw.SPM, pch=15, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.4)

points (ET.diaria ~ gw.SPM, pch=25, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="black", bg="black",
        data=pot.5)


abline (h=0 , lwd=1, col="black")
abline (v=0 , lwd=1, col="black")

3^0



mod.1 <- lm (gw.SPM ~ w.SPM, data=pot.1 )
summary(mod.1)
abline (mod.1 , lwd=1, col="black")



mod.2 <- lm (gw.SPM ~ w.SPM, data=pot.2 )
summary(mod.2)
abline (mod.2 , lwd=1, col="blue")

points (gw.SPM ~ w.SPM, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="red", bg="red",
        data=pot.3)

mod.3 <- lm (gw.SPM ~ w.SPM, data=pot.3 )
summary(mod.3)
abline (mod.3 , lwd=1, col="red")







mod.4 <- lm (gw.SPM ~ w.SPM, data=pot.4 )
summary(mod.4)
abline (mod.4 , lwd=1, col="darkgreen")

points (gw.SPM ~ w.SPM, pch=16, type="p", 
        lty=3, lwd=1.5, cex=1,
        col="orange", bg="orange",
        data=pot.5)

mod.5 <- lm (gw.SPM ~ w.SPM, data=pot.5 )
summary(mod.5)
abline (mod.5 , lwd=1, col="orange")
        
axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("GW (mmol m-2 s-1 )",side=2,line=2.5)
#axis (1,500:800)
mtext("water (g)",side=1,line=2.5)
box()



mod.1 <- lm (gw.SPM ~ w.SPM, data=pot.1 )
summary(mod.1)
abline (mod.1 , lwd=2, col="black")

et.model.2 <- function(x, y) {
  y1 <- filter (x, time==y)    %>% select(ET)
  y2 <- filter (x, time==2*y)  %>% select(ET)
  n <- 2*y - y
  Aest <-  (y1^2) / (2*y1 - y2)
  kest <-  (-1/n) * log((y2/y1) - 1)
  par.est <- cbind (k.est=kest[1,1],A.est=Aest[1,1] )
  ET.par <- nls (ET ~ A * (1 - exp (-k*time)),
                 start = list (k = kest[1,1], A = Aest[1,1]),
                 trace = TRUE , 
                 data= x)
  zz <- rbind (coef (ET.par))
  par.mod <- data.frame (pot=x[1,1],zz, par.est)
  print(par.mod)
}
