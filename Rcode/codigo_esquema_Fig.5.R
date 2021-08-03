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
library("dplyr")
library("ggridges")
library("ggjoy")
library("viridis")
library("lmerTest")
library("lubridate")
library(ggcorrplot)

########  se carga los datos 
## ensayo 1
pred.model.geno.DMSO.1 <- read.table ("./Data/rawdata/model_geno_DMSO_1.txt", 
                                    header = TRUE, sep = "\t",dec = ".",
                                    na.strings = "NA")
str(pred.model.geno.DMSO.1)

pred.model.geno.DMSO.1 <- pred.model.geno.DMSO.1%>%
                          rename (genotype = lineas)


## plot las correlaciones entre las variables estimadas 
ggplot (data = pred.model.geno.DMSO.1 )+
  geom_point(mapping = aes (x= k.eq.07,y=ks))

cor (pred.model.geno.DMSO.1$k.eq.07, pred.model.geno.DMSO.1$ks)

ggplot (data = pred.model.geno.DMSO.1 )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.DMSO.1$B.eq.07, pred.model.geno.DMSO.1$Bs)

px1 <- pred.model.geno.DMSO.1 %>%
       filter (B.eq.07 > 1000) 

pred.model.geno.DMSO.1a <- pred.model.geno.DMSO.1 %>%
                          filter (pot!=99)

ggplot (data = pred.model.geno.DMSO.1a )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.DMSO.1a$B.eq.07, pred.model.geno.DMSO.1a$Bs)

px2 <- pred.model.geno.DMSO.1a %>%
      filter (B.eq.07 > 400)

pred.model.geno.DMSO.1b <- pred.model.geno.DMSO.1a %>%
                           filter (pot!=174)

ggplot (data = pred.model.geno.DMSO.1b )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.DMSO.1b$B.eq.07, pred.model.geno.DMSO.1b$Bs)


## ensayo 2
pred.model.geno.DMSO.2 <- read.table ("./Data/rawdata/model_geno_DMSO_2.txt", 
                                      header = TRUE, sep = "\t",dec = ".",
                                      na.strings = "NA")
summary (pred.model.geno.DMSO.2)

pred.model.geno.DMSO.1b <- pred.model.geno.DMSO.1b %>%
                           arrange (genotype)

pred.model.geno.DMSO.2 <- pred.model.geno.DMSO.2 %>%
                          arrange (genotype)

pred.model.geno.DMSO.1b$genotype <-  as.character (pred.model.geno.DMSO.1b$genotype)
pred.model.geno.DMSO.2$genotype <-  as.character(pred.model.geno.DMSO.2$genotype)

setdiff (pred.model.geno.DMSO.1b$genotype, pred.model.geno.DMSO.2$genotype)  
setdiff (pred.model.geno.DMSO.2$genotype, pred.model.geno.DMSO.1b$genotype)

pred.model.geno.DMSO.2a <- pred.model.geno.DMSO.2 %>%
                          filter ( genotype !="pl.107" &
                                   genotype !="pl.117" &
                                   genotype !="pl.121" &
                                   genotype !="pl.124" &
                                   genotype !="pl.140" &
                                   genotype !="pl.145" &
                                   genotype !="pl.150" &
                                   genotype !="pl.23"  &
                                   genotype !="pl.74"  &
                                   genotype !="pl.80"  &
                                   genotype !="pl.95" )


pred.model.geno.DMSO.1b$genotype == pred.model.geno.DMSO.2a$genotype

head(pred.model.geno.DMSO.1b)
head(pred.model.geno.DMSO.2a)

pred.model.geno.DMSO.12 <- pred.model.geno.DMSO.1b %>%
                           mutate ( pred.a = pred.model.geno.DMSO.2a$a)%>%
                           mutate ( pred.b = pred.model.geno.DMSO.2a$b)%>%
                           mutate ( pred.rsq = pred.model.geno.DMSO.2a$rsq)

pred.model.geno.DMSO.12.a <- pred.model.geno.DMSO.12 %>%
                             filter (pred.rsq > 0.51)

####### plot de loas datos crudos ###################
ggplot (data = pred.model.geno.DMSO.12.a )+
        geom_point(mapping = aes (x= k.eq.07,y=ks))

cor (pred.model.geno.DMSO.12.a$k.eq.07, pred.model.geno.DMSO.12.a$ks)

ggplot (data =pred.model.geno.DMSO.12.a )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.DMSO.12.a$B.eq.07, pred.model.geno.DMSO.12.a$Bs)

cor (pred.model.geno.DMSO.12.a$ks, pred.model.geno.DMSO.12.a$Bs)

#################  xy plot par K k.est ################
head (pred.model.geno.DMSO.12.a)

pred.model.geno.DMSO.12.a <- pred.model.geno.DMSO.12.a %>%
                             mutate (tm = (1/pred.model.geno.DMSO.12.a$ks) * log(2)) %>%
                             mutate (ab= -pred.model.geno.DMSO.12.a$pred.a/pred.model.geno.DMSO.12.a$pred.b)


pred.model.geno.DMSO.12.b <- pred.model.geno.DMSO.12.a %>%
                             filter (tm < 11)
summary (pred.model.geno.DMSO.12.b)

pred.model.geno.DMSO.12.c <- pred.model.geno.DMSO.12.b %>%
                             filter (ab > 400 )

pred.model.geno.DMSO.12.d <-  pred.model.geno.DMSO.12.c %>%
                              mutate (dGwtm = (2*pred.model.geno.DMSO.12.c$ks*
                                               pred.model.geno.DMSO.12.c$pred.b* 
                                               pred.model.geno.DMSO.12.c$Bs))

head(pred.model.geno.DMSO.12.c)


summary (pred.model.geno.DMSO.12.c)

## ensayo 3
pred.model.geno.DMSO.3 <- read.table ("./Data/rawdata/model_geno_DMSO_3.txt", 
                                      header = TRUE, sep = "\t",dec = ".",
                                      na.strings = "NA")

pred.model.geno.DMSO.3 <- pred.model.geno.DMSO.3 %>%
                          arrange (genotype)



pred.model.geno.DMSO.12.d <- pred.model.geno.DMSO.12.d %>%
                             arrange (genotype)


dim(pred.model.geno.DMSO.3)
dim(pred.model.geno.DMSO.12.d)

pred.model.geno.DMSO.3$genotype <-  as.character(pred.model.geno.DMSO.3$genotype)


setdiff (pred.model.geno.DMSO.12.d$genotype, pred.model.geno.DMSO.3$genotype)  
setdiff (pred.model.geno.DMSO.3$genotype, pred.model.geno.DMSO.12.d$genotype)

pred.model.geno.DMSO.12.e <- pred.model.geno.DMSO.12.d %>%
                             filter ( genotype !="pl.105" & 
                                    genotype !="pl.108" & 
                                    genotype !="pl.112" & 
                                    genotype !="pl.115" & 
                                    genotype !="pl.126" & 
                                    genotype !="pl.130" & 
                                    genotype !="pl.132" & 
                                    genotype !="pl.133" & 
                                    genotype !="pl.138" & 
                                    genotype !="pl.142" & 
                                    genotype !="pl.153" & 
                                    genotype !="pl.154" & 
                                    genotype !="pl.155" &
                                    genotype !="pl.158" & 
                                    genotype !="pl.163" & 
                                    genotype !="pl.166" & 
                                    genotype != "pl.168" & 
                                    genotype !="pl.170"  & 
                                    genotype !="pl.172"  &
                                    genotype !="pl.174"  &
                                    genotype !="pl.176"  & 
                                    genotype !="pl.177"  & 
                                    genotype !="pl.179"  & 
                                    genotype !="pl.183"  &
                                    genotype !="pl.184"  &
                                    genotype !="pl.19"   & 
                                    genotype !="pl.21"    & 
                                    genotype !="pl.24"   & 
                                    genotype !="pl.31"   & 
                                    genotype !="pl.32"   & 
                                    genotype !="pl.33"  & 
                                    genotype !="pl.36"  & 
                                    genotype !="pl.37" & 
                                    genotype !="pl.46"  & 
                                    genotype !="pl.51" & 
                                    genotype !="pl.52" & 
                                    genotype !="pl.54"  &
                                    genotype !="pl.56" & 
                                    genotype !="pl.61"  & 
                                    genotype !="pl.62"  & 
                                    genotype != "pl.64" & 
                                    genotype !="pl.67" & 
                                    genotype !="pl.68"  & 
                                    genotype !="pl.7"   & 
                                    genotype !="pl.73" & 
                                    genotype !="pl.78" & 
                                    genotype !="pl.83" & 
                                    genotype !="pl.85" & 
                                    genotype != "pl.86" & 
                                    genotype !="pl.87" & 
                                    genotype !="pl.88" & 
                                    genotype !="pl.90" & 
                                    genotype != "pl.93" & 
                                    genotype !="pl.94" & 
                                    genotype !="pl.96" & 
                                    genotype !="pl.97")


pred.model.geno.DMSO.3a <-  pred.model.geno.DMSO.3 %>%
                              filter ( genotype !="pl.1" &
                                       genotype !="pl.101"  &
                                       genotype !="pl.114"  &
                                       genotype !="pl.139"  &
                                       genotype !="pl.144" &
                                       genotype !="pl.152" &
                                       genotype !="pl.159" &
                                       genotype !="pl.27" & 
                                       genotype != "pl.3" &
                                       genotype !="pl.39"  &
                                       genotype !="pl.63" &
                                       genotype !="pl.65"  &
                                       genotype != "pl.75"  & 
                                       genotype !="SO76557")


pred.model.geno.DMSO.3a <- pred.model.geno.DMSO.3a %>% 
                           arrange (genotype)  
str(pred.model.geno.DMSO.3a)
                                         
pred.model.geno.DMSO.12.e <-pred.model.geno.DMSO.12.e %>%
                            arrange (genotype)
                                         
pred.model.geno.DMSO.12.e$genotype == pred.model.geno.DMSO.3a$genotype

str(pred.model.geno.DMSO.12.e)


########### armo la matriz para graficar 
b <- pred.model.geno.DMSO.12.e$pred.b
a <- pred.model.geno.DMSO.12.e$pred.a
Bs <-  pred.model.geno.DMSO.12.e$Bs
ks <-  pred.model.geno.DMSO.12.e$ks
tm <-  pred.model.geno.DMSO.12.e$tm
ARS <- pred.model.geno.DMSO.3a$ARS

#### aca se calcula la conductancia en t0.5

pred.model.geno.DMSO.123 <- pred.model.geno.DMSO.12.e %>%
                            mutate (k.eq.xy = pred.model.geno.DMSO.3a$k.eq.xy) %>%
                            mutate (B.eq.xy = pred.model.geno.DMSO.3a$B.eq.xy) %>%
                            mutate (AR =pred.model.geno.DMSO.3a$AR) %>%
                            mutate (tm.AR = (1/pred.model.geno.DMSO.3a$k.eq.xy) * log(2)) %>%
                            mutate (Gwtm = b*Bs*(exp (-ks *tm)) + b * ARS + a )

#View (pred.model.geno.DMSO.123)

write.table (pred.model.geno.DMSO.123, file = "./Data/procdata/pred.model.geno.DMSO.123", 
             append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

summary (pred.model.geno.DMSO.123$genotype)

head(pred.model.geno.DMSO.123)



ggplot (data =pred.model.geno.DMSO.123)+
  geom_point(mapping = aes (x=tm.AR, y=tm))

cor (pred.model.geno.DMSO.123$tm.AR, pred.model.geno.DMSO.123$tm)

##################
# PRUEBA DE NORMALIDAD
svg (filename="./Figures/Fig.5/fig.5a.hist.tm.svg", 
     width=7, 
     height=5, 
     pointsize=12)

tm <- pred.model.geno.DMSO.123$tm

q25.tm <-  quantile (tm, .25, na.rm = TRUE)
q75.tm <-  quantile (tm, .75, na.rm = TRUE)
q50.tm <-  quantile (tm, .50, na.rm = TRUE)

#### aca se agrupa por el tm ####
tm.nq.q1 <- pred.model.geno.DMSO.123 %>%
            dplyr::filter (tm < q25.tm) %>%
            dplyr::arrange (tm) %>%
            dplyr::mutate (group.tm = "tm.q1" )


tm.nq.q2 <- pred.model.geno.DMSO.123 %>%
            filter (tm >= q25.tm) %>%
            filter (tm <= q75.tm)%>%
            arrange (tm) %>%
            dplyr::mutate (group.tm = "tm.q2" )


tm.nq.q3 <- pred.model.geno.DMSO.123 %>%
            filter (tm > q75.tm)%>%
            arrange (tm) %>%
            dplyr::mutate (group.tm = "tm.q3" )


pred.model.geno.DMSO.123.qtm <- bind_rows(tm.nq.q1, tm.nq.q2, tm.nq.q3)


h.tm  <- hist (tm , breaks ="Sturges")
h.tm$density <- h.tm$counts/sum (h.tm$counts)*100


plot (h.tm, freq=FALSE, xlab="em.tm",main = "",
      ylab="Percent of total", col=NULL
      #,xlim = c(400,1000)
      ,xaxt='l')

par(new=TRUE)

hist (tm, main="tm.geno", 
      freq=F, xlab="", ylab = "",
      axes=FALSE,col="gray85")

lines (density (tm, na.rm = T), lwd=2,
       col="navyblue", lty=1)

dens.tm <- density (tm, na.rm = T)
x1.tm <- min (which (dens.tm$x >= 0))
x2.tm <- max (which (dens.tm$x <= q25.tm))

with (dens.tm, polygon(x=c(x[c(x1.tm,x1.tm:x2.tm,x2.tm)]),
                         y= c(0, y[x1.tm:x2.tm], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue', 0.5)))

x3.tm  <- min (which (dens.tm$x >= q75.tm))
x4.tm  <- max (which (dens.tm$x > 0))
with (dens.tm, polygon(x=c(x[c(x3.tm,x3.tm:x4.tm,x4.tm)]),
                         y= c(0, y[x3.tm:x4.tm], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',1)))

abline (v=q50.tm, col="navyblue", lwd=1.5, lty =2)


box()
dev.off()


########### 
svg (filename="./Figures/Fig.5/fig.5b.hist.Gwtm.svg", 
     width=7, 
     height=5, 
     pointsize=12)

Gwtm <- pred.model.geno.DMSO.123$Gwtm
shapiro.test (Gwtm)

summary (Gwtm)
q25.Gwtm <-  quantile (Gwtm, .25, na.rm = TRUE)
q75.Gwtm <-  quantile (Gwtm, .75, na.rm = TRUE)
q50.Gwtm <-  quantile (Gwtm, .50, na.rm = TRUE)

Gwtm.nq.q1 <- pred.model.geno.DMSO.123 %>%
              dplyr::filter (Gwtm < q25.Gwtm)%>%
              dplyr::arrange (Gwtm) %>%
              dplyr::mutate (group.Gwtm = "Gwtm.q1" )


Gwtm.nq.q2 <- pred.model.geno.DMSO.123 %>%
              dplyr::filter (Gwtm  >= q25.Gwtm) %>%
              dplyr::filter (Gwtm  <=  q75.Gwtm)%>%
              dplyr::arrange (Gwtm) %>%
              dplyr::mutate (group.Gwtm = "Gwtm.q2" )

Gwtm.nq.q3 <- pred.model.geno.DMSO.123 %>%
              dplyr::filter (Gwtm > q75.Gwtm)%>%
              dplyr::arrange (Gwtm)%>%
              dplyr::mutate (group.Gwtm = "Gwtm.q3" )


pred.model.geno.DMSO.123.qGw <- bind_rows(Gwtm.nq.q1, Gwtm.nq.q2, Gwtm.nq.q3)

h.Gwtm  <- hist (Gwtm , breaks="Sturges")
h.Gwtm$density <- h.Gwtm$counts/sum (h.Gwtm$counts)*100


plot (h.Gwtm, freq=FALSE, xlab="Gwtm",main = "",
      ylab="Percent of total",
      col=NULL
      #,xlim = c(400,1000)
      ,xaxt='l')

par(new=TRUE)

hist (Gwtm, main="Gwtm.geno", 
      freq=F, xlab="", ylab = "",
      axes=FALSE,col="gray85")

lines (density (Gwtm, na.rm = T), lwd=2,
       col="navyblue", lty=1)

dens.Gwtm <- density (Gwtm, na.rm = T)
x1.Gwtm <- min (which (dens.Gwtm$x >= 0))
x2.Gwtm <- max (which (dens.Gwtm$x <= q25.Gwtm))

with (dens.Gwtm, polygon (x=c(x[c(x1.Gwtm,x1.Gwtm:x2.Gwtm,x2.Gwtm)]),
                          y= c(0, y[x1.Gwtm:x2.Gwtm], 0),border ='navyblue',
                          col=scales::alpha( 'navyblue',.5)))

x3.Gwtm  <- min (which (dens.Gwtm$x >= q75.Gwtm))
x4.Gwtm  <- max (which (dens.Gwtm$x > 0))
with (dens.Gwtm, polygon(x=c(x[c(x3.Gwtm,x3.Gwtm:x4.Gwtm,x4.Gwtm)]),
                          y= c(0, y[x3.Gwtm:x4.Gwtm], 0),border ='navyblue',
                          col=scales::alpha( 'navyblue',1)))

abline (v=q50.Gwtm, col="navyblue", lwd=1.5, lty =2)

box()
dev.off()

pred.model.geno.DMSO.123.qGw.1 <- pred.model.geno.DMSO.123.qGw %>%
                                dplyr::select (c(pot, group.Gwtm))

group.feno.DMSO <- pred.model.geno.DMSO.123.qtm %>%
                   dplyr::inner_join(pred.model.geno.DMSO.123.qGw.1, by="pot" )

write.table (group.feno.DMSO, file = "./Data/procdata/pred.model.geno.DMSO.123", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

write_delim (group.feno.DMSO, file = "./Data/procdata/group.feno.DMSO.txt",
             delim = ";",na = "NA")

###### ACA TERMINA EL CODIGO PARA LA FIGURA 5A Y 5B ##############################

########  se carga los datos #########
# cargar datos crudos
ec2brep.d.model.4 <- read.table ("./Data/rawdata/ec2brep.d.model.4.txt" ,
                                 header = TRUE, sep = ",",dec = ".",
                                 na.strings = "NA" )
head(ec2brep.d.model.4)

xx <- ec2brep.d.model.4 %>%
      dplyr::filter (str_detect(linea.disenio , "testigo")) %>%
      dplyr::select (genotype,linea.disenio)


distinct(xx, genotype, .keep_all = TRUE)



## esta son las predichos de los parametros de los predichos de las macetas 

pred.model.geno.EC.1 <- read.table ("./Data/rawdata/model_geno_EC_1.txt", 
                                    header = TRUE, sep = "\t",dec = ".",
                                    na.strings = "NA")
head(pred.model.geno.EC.1)

pred.model.geno.EC.1a <- pred.model.geno.EC.1 %>%
                         dplyr::filter ( genotype !="111_4"  &
                                         genotype != "133_TARD_38" &
                                         genotype !="133_TARD_8")

pred.model.geno.EC.2 <- read.table ("./Data/rawdata/model_geno_EC_2.txt", 
                                    header = TRUE, sep = "\t",dec = ".",
                                    na.strings = "NA")

pred.model.geno.EC.2a <- pred.model.geno.EC.2 %>%
                         dplyr::filter (genotype != "133_TARD_38")




setdiff ( pred.model.geno.EC.1a$genotype, pred.model.geno.EC.2a$genotype)  
setdiff ( pred.model.geno.EC.2a$genotype, pred.model.geno.EC.1a$genotype)  


pred.model.geno.EC.1a$genotype <-  as.character(pred.model.geno.EC.1a$genotype)
pred.model.geno.EC.2a$genotype  <-  as.character(pred.model.geno.EC.2a$genotype)


pred.model.geno.EC.1a$genotype ==  pred.model.geno.EC.2a$genotype


pred.model.geno.EC.12 <- pred.model.geno.EC.1a %>%
                         dplyr::mutate ( pred.a = pred.model.geno.EC.2a$a)%>%
                         dplyr::mutate ( pred.b = pred.model.geno.EC.2a$b)%>%
                         dplyr::mutate ( pred.rsq = pred.model.geno.EC.2a$rsq)

summary (pred.model.geno.EC.12)

pred.model.geno.EC.12a <- pred.model.geno.EC.12 %>%
                          dplyr::filter ( pred.rsq > 0.51) 

summary (pred.model.geno.EC.12a)


######## plot de loas datos crudos ###################
ggplot (data = pred.model.geno.EC.12a )+
  geom_point(mapping = aes (x= k.eq.07,y=ks))

cor (pred.model.geno.EC.12a$k.eq.07, pred.model.geno.EC.12a$ks)

ggplot (data =pred.model.geno.EC.12a )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.EC.12a$B.eq.07, pred.model.geno.EC.12a$Bs)


pred.model.geno.EC.12b <- pred.model.geno.EC.12a %>%
  filter (B.eq.07 < 1100)

ggplot (data =pred.model.geno.EC.12b )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.EC.12b$B.eq.07, pred.model.geno.EC.12b$Bs)

px1 <- pred.model.geno.EC.12b %>%
  filter (B.eq.07 > 1000) %>%
  filter (Bs < 600) 

pred.model.geno.EC.12c <- pred.model.geno.EC.12b %>%
  filter (genotype !="101_39")


ggplot (data =pred.model.geno.EC.12c )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.EC.12c$B.eq.07, pred.model.geno.EC.12c$Bs)


px2 <- pred.model.geno.EC.12c %>%
  filter (B.eq.07 > 750) %>%
  filter (Bs < 550) 

pred.model.geno.EC.12d <- pred.model.geno.EC.12c %>%
  filter (genotype !="85_2_8_20")


ggplot (data =pred.model.geno.EC.12d )+
  geom_point(mapping = aes (x=B.eq.07, y=Bs))

cor (pred.model.geno.EC.12d$B.eq.07, pred.model.geno.EC.12d$Bs)

cor(pred.model.geno.EC.12d$ks, pred.model.geno.EC.12d$Bs)
head(pred.model.geno.EC.12d)


pred.model.geno.EC.12d  <- pred.model.geno.EC.12d %>%
                           dplyr::mutate (ab= -pred.model.geno.EC.12d$pred.a/pred.model.geno.EC.12d$pred.b)%>%
                           dplyr:: mutate (dGwtm = (2*pred.model.geno.EC.12d$ks*
                                                     pred.model.geno.EC.12d$pred.b* 
                                                     pred.model.geno.EC.12d$Bs))


################### Aca voy a cargar los datos de agua remanente #
pred.model.geno.EC.3 <- read.table ("./Data/rawdata/model_geno_EC_3.txt", 
                                    header = TRUE, sep = "\t",dec = ".",
                                    na.strings = "NA")

pred.model.geno.EC.3 <- pred.model.geno.EC.3 %>%
                        dplyr::arrange (genotype)



pred.model.geno.EC.12d <- pred.model.geno.EC.12d %>% 
                          dplyr::arrange (genotype)


dim(pred.model.geno.EC.3)
dim(pred.model.geno.EC.12d)

pred.model.geno.EC.3$genotype <-  as.character(pred.model.geno.EC.3$genotype)


setdiff (pred.model.geno.EC.12d$genotype, pred.model.geno.EC.3$genotype)  
setdiff (pred.model.geno.EC.3$genotype, pred.model.geno.EC.12d$genotype)

pred.model.geno.EC.12e <- pred.model.geno.EC.12d %>%
                          dplyr::filter (genotype !="118_26" &
                                         genotype !="133_TARD_21" &
                                         genotype !="YM10_36005")

pred.model.geno.EC.3a <- pred.model.geno.EC.3 %>%
                         dplyr::filter (genotype !="100_30" &
                                        genotype !="101_3"  &
                                        genotype !="101_31" &
                                        genotype !="111_13" &
                                        genotype !="117_23" &
                                        genotype != "117_9" &
                                        genotype != "118_23" &
                                        genotype != "124_19" &
                                        genotype != "127_9" &
                                        genotype != "133_TARD_32" &
                                        genotype != "133_TARD_38" &
                                        genotype != "139_30" &
                                        genotype != "85_2_8_20")                               


pred.model.geno.EC.3a$genotype == pred.model.geno.EC.12e$genotype                             

########### armo la matriz paragraficar 
head(pred.model.geno.EC.12e)
b <- pred.model.geno.EC.12e$pred.b
a <- pred.model.geno.EC.12e$pred.a
Bs <-  pred.model.geno.EC.12e$Bs
ks <-  pred.model.geno.EC.12e$ks
tm <-  pred.model.geno.EC.12e$tm
ARS <- pred.model.geno.EC.3a$ARS


pred.model.geno.ES.123 <- pred.model.geno.EC.12e %>%
  mutate ( k.eq.xy = pred.model.geno.EC.3a$k.eq.xy) %>%
  mutate ( B.eq.xy = pred.model.geno.EC.3a$B.eq.xy) %>%
  mutate ( AR =pred.model.geno.EC.3a$AR) %>%
  mutate (tm.AR = (1/pred.model.geno.EC.3a$k.eq.xy) * log(2)) %>%
  mutate (Gwtm = b*Bs*(exp (-ks *tm)) + b * ARS + a )


write.table (pred.model.geno.ES.123, file = "./Data/procdata/pred.model.geno.ES.123.txt",
              append = FALSE, quote = TRUE, sep = ",",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE)




# PRUEBA DE NORMAILDAD
svg (filename="./Figures/Fig.5/fig.5c.hist.tm.svg", 
     width=7, 
     height=5, 
     pointsize=12)

tm <- pred.model.geno.ES.123$tm
shapiro.test (tm)

#q25.kest <-  quantile (kest, .25, na.rm = TRUE)
#q75.kest <-  quantile (kest, .75, na.rm = TRUE)

q25.tm <-  quantile (tm, .25, na.rm = TRUE)
q75.tm <-  quantile (tm, .75, na.rm = TRUE)
q50.tm <-  quantile (tm, .50, na.rm = TRUE)
###############################################33
tm.nq.q1 <- pred.model.geno.ES.123  %>%
  dplyr::filter (tm < q25.tm) %>%
  dplyr::arrange (tm) %>%
  dplyr::mutate (group.tm = "tm.q1" )


tm.nq.q2 <- pred.model.geno.ES.123 %>%
  filter (tm >= q25.tm) %>%
  filter (tm <= q75.tm)%>%
  arrange (tm) %>%
  dplyr::mutate (group.tm = "tm.q2" )


tm.nq.q3 <- pred.model.geno.ES.123 %>%
  filter (tm > q75.tm)%>%
  arrange (tm) %>%
  dplyr::mutate (group.tm = "tm.q3" )


pred.model.geno.Elite.123.qtm <- bind_rows(tm.nq.q1, tm.nq.q2, tm.nq.q3)

##############################3
h.tm  <- hist (tm , breaks ="Sturges")
h.tm$density <- h.tm$counts/sum (h.tm$counts)*100


plot (h.tm, freq=FALSE, xlab="em.tm",main = "",
      ylab="Percent of total", col=NULL
      #,xlim = c(400,1000)
      ,xaxt='l')

par(new=TRUE)

hist (tm, main="tm.geno", 
      freq=F, xlab="", ylab = "",
      axes=FALSE,col="gray85")
#dev.off()
lines (density (tm, na.rm = T), lwd=2,
       col="navyblue", lty=1)

dens.tm <- density (tm, na.rm = T)
x1.tm <- min (which (dens.tm$x >= 0))
x2.tm <- max (which (dens.tm$x <= q25.tm))

with (dens.tm, polygon(x=c(x[c(x1.tm,x1.tm:x2.tm,x2.tm)]),
                       y= c(0, y[x1.tm:x2.tm], 0),border ='navyblue',
                       col=scales::alpha( 'navyblue')))

x3.tm  <- min (which (dens.tm$x >= q75.tm))
x4.tm  <- max (which (dens.tm$x > 0))
with (dens.tm, polygon(x=c(x[c(x3.tm,x3.tm:x4.tm,x4.tm)]),
                       y= c(0, y[x3.tm:x4.tm], 0),border ='navyblue',
                       col=scales::alpha( 'navyblue')))

abline (v=q50.tm, col="navyblue", lwd=1.5, lty =2)
box()

dev.off()


####################
svg (filename="./Figures/Fig.5/fig.5d.hist.Gwtm.svg", 
     width=7, 
     height=5, 
     pointsize=12)

Gwtm <- pred.model.geno.ES.123$Gwtm
shapiro.test (Gwtm)

q25.Gwtm <-  quantile (Gwtm, .25, na.rm = TRUE)
q75.Gwtm <-  quantile (Gwtm, .75, na.rm = TRUE)
q50.Gwtm <-  quantile (Gwtm, .50, na.rm = TRUE)

Gwtm.nq.q1 <- pred.model.geno.ES.123  %>%
               dplyr::filter (Gwtm < q25.Gwtm ) %>%
               dplyr::arrange (Gwtm) %>%
               dplyr::mutate (group.Gwtm = "Gwtm.q1" )


Gwtm.nq.q2 <- pred.model.geno.ES.123 %>%
              dplyr::filter (Gwtm >= q25.Gwtm) %>%
              dplyr::filter (Gwtm <= q75.Gwtm) %>%
              dplyr::arrange (Gwtm) %>%
              dplyr::mutate (group.Gwtm = "Gwtm.q2")


Gwtm.nq.q3 <- pred.model.geno.ES.123 %>%
            dplyr::filter (Gwtm > q75.Gwtm)%>%
            dplyr::arrange (Gwtm) %>%
            dplyr::mutate (group.Gwtm = "Gwtm.q3" )


pred.model.geno.Elite.123.qGwtm <- bind_rows(Gwtm.nq.q1, Gwtm.nq.q2, Gwtm.nq.q3)


h.Gwtm  <- hist (Gwtm , breaks="Sturges")
h.Gwtm$density <- h.Gwtm$counts/sum (h.Gwtm$counts)*100


plot (h.Gwtm, freq=FALSE, xlab="Gwtm",main = "",
      ylab="Percent of total",
      col=NULL
      #,xlim = c(400,1000)
      ,xaxt='l')

par(new=TRUE)

hist (Gwtm, main="Gwtm.geno", 
      freq=F, xlab="", ylab = "",
      axes=FALSE,col="gray85")

lines (density (Gwtm, na.rm = T), lwd=2,
       col="navyblue", lty=1)

dens.Gwtm <- density (Gwtm, na.rm = T)
x1.Gwtm <- min (which (dens.Gwtm$x >= 0))
x2.Gwtm <- max (which (dens.Gwtm$x <= q25.Gwtm))

with (dens.Gwtm, polygon(x=c(x[c(x1.Gwtm,x1.Gwtm:x2.Gwtm,x2.Gwtm)]),
                         y= c(0, y[x1.Gwtm:x2.Gwtm], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))

x3.Gwtm  <- min (which (dens.Gwtm$x >= q75.Gwtm))
x4.Gwtm  <- max (which (dens.Gwtm$x > 0))
with (dens.Gwtm, polygon(x=c(x[c(x3.Gwtm,x3.Gwtm:x4.Gwtm,x4.Gwtm)]),
                         y= c(0, y[x3.Gwtm:x4.Gwtm], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))

abline (v=q50.Gwtm, col="navyblue", lwd=1.5, lty =2)

box()

dev.off()

pred.model.geno.Elite.123.qGw.1 <- pred.model.geno.Elite.123.qGwtm %>%
                                   dplyr::select (c(pot, group.Gwtm))

group.feno.Elite <- pred.model.geno.Elite.123.qtm %>%
                    dplyr::inner_join(pred.model.geno.Elite.123.qGw.1, by="pot" )



write_delim (group.feno.Elite, file = "./Data/procdata/group.feno.Elite.txt",
             delim = ";",na = "NA")



############### ACA termina el codigo para la figura 5 ##########################