#######################################################################
# Codigo esquema paper soja Frontiers in Plants Science               #
# Datos  Soja estres hidrico en Cámara de crecimiento                 #
# Datos tomados por Esteban Casaretto                                 #
# Modelizacion de ET y Gw de la poblacion de GWAS                     #
#                                                                     #
# Sebastian Simondi - Gaston Quero                                    #
# 6-09-2018                                                           #
#######################################################################

getwd ()
#setwd ("E:/Paper_Soja_FPS")

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
library ("dplyr")
library (ggjoy)
library ("ggridges")
library (hrbrthemes)
library (forcats)
library ("viridis")
library ("lmerTest")
library (nlme)
library (xtable)
library (stringr)
library (data.table)
library (svMisc)
library (ggpubr)
library ("ggsci")
library ("FactoMineR")
library ("factoextra")
library ("corrplot")
library ("readr")


########  se carga los datos 

SPM_data <- read_delime (file ="./Data/rawdata/SPM_data.txt" ,
                        delim="\t", na="NA",
                        escape_backslash = FALSE,
                        escape_double = TRUE,
                        col_names = TRUE,
                        col_types = NULL)








str (SPM_data)

summary (SPM_data)





# cargar datos crudos
# ec2b.rep
ec2b.rep <- read.table ("./Data/rawdata/input/ec2b.rep.txt" ,
                        header = TRUE, sep = "",dec = ".",
                        na.strings = "NA" )

###### control de la matriz #########3
str (ec2b.rep)
length(unique(ec2b.rep$pot))
length(unique(ec2b.rep$genotype))


####### esto es para setear la salida grafica
##########  argumentos ########
dt = ec2b.rep 
pss = 1970 
pm = 215 
t1 =4 
W="w.SPM"
essay="ec2b.rep1"
ue ="pot"
genotype = "genotype"

run_eq.07 <- function (dt = NULL,  genotype=NULL, t1=NULL, W=NULL, pss=NULL, pm = NULL, essay=NULL, ue=NULL) {
 
   dir.create (file.path ("Data"), showWarnings = FALSE)
   dir.create (file.path ("Data", "outpout.ET"), showWarnings = FALSE)
   dir.create (file.path ("Figures"), showWarnings = FALSE)
   dir.create (file.path ("Figures", "Plots.ET"), showWarnings = FALSE)

   dt1 <- dt %>% 
          dplyr::rename (genotype = all_of (genotype)) %>%
          dplyr::rename (W = all_of (W)) %>%
          dplyr::rename (pot = all_of (ue))
#getwd()
#setwd (file.path ("Data", "outpout"))
#listapot <- unique(ec2b.rep1$pot)  ### aca podria ser ue y a la reemplazo por dt

resultados <- NULL
datos_ET <- NULL

listapot <- unique(dt1$pot)

dt.ET <- lapply (listapot, function(filtro){
         
         datos <- filter (dt1, pot== filtro)
         S <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
         datos$water.SPM.t <-  datos$W -S 
         W0 = datos$W[which(datos$time == 0)]
         datos$ET.t <- W0 - datos$W
  
         write.table (datos, file =paste("./Data/outpout.ET/dt.ET","pot_",filtro,".txt",sep=""),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = FALSE,
               col.names = TRUE)
          datos_ET <- datos 
          print (datos_ET) 
})

datos.ET <- do.call(rbind.data.frame, dt.ET)

write.table (datos.ET, file =paste("./Data/outpout.ET/dt.ET",essay,".csv",sep=""),
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)


et.model  <- sapply (listapot, function (filtro) {        # aca reeplazaria por listaue
                  datos <- filter (dt1, pot == filtro) # aca reemplaze por dt ec2b.rep
                  S     <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
                  datos$water.SPM.t <-  datos$W - S 
                  W0 <- datos$W[which(datos$time == 0)]
                  datos$ET.t <- W0 - datos$W

#W0 = datos$w.SPM[which(datos$time == 0)]
Wn =  datos$W[which(datos$time == t1)]
W2n = datos$W[which(datos$time == 2*t1)]
n = t1
Bs <-  (Wn - W0)^2 /(W2n - 2*Wn + W0)                     # este es el "B"  estimado por Simondi
ks <-  -(1/n) * log ((W2n-Wn)/ (Wn-W0))                   # este es el "k"  estimado por Simondi

if (is.na (ks) == FALSE) {

if(ks > 0.05){ 
  
resultados <- cbind(ks,Bs)

eq.07 <- nls (ET.t ~ B * (1 - exp (-k*time)),               # esta es la estimacion de los parametros por Gauss-Newton
              start = list (k = ks, B = Bs),
              trace = FALSE , 
              data= datos, nls.control(warnOnly=TRUE))

k.eq.07 <- round (coef (eq.07)[1],3)                            # este es el k de Gauss-Newton
B.eq.07 <- round (coef (eq.07)[2],3)    
tm <- (1/ k.eq.07) * log(2)    

resultados <- cbind(B.eq.07,resultados)
resultados <- cbind(k.eq.07,resultados)

es   <-  data.frame (datos$essay[1])
geno <-  data.frame(datos$genotype[1])
pot  <-  filtro

resultados <- cbind(k.eq.07,B.eq.07)

es   <-  data.frame(datos$essay[1])
geno <-  data.frame(datos$genotype[1])
pot  <-  filtro

resultados <- data.frame (es,geno,pot, resultados)
resultados <- resultados %>%
  rename (essay = datos.essay.1.)%>%
  rename (genotype = datos.genotype.1.) %>%
  mutate (tm = tm)



write.table (resultados, file =paste("./Data/outpout.ET/par_ET_model","pot_",filtro,".txt",sep=""),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = FALSE,
               col.names = TRUE)
print(resultados)

  }
 }
})

resultado.global <- do.call (rbind.data.frame, et.model)

write.table (resultado.global, file =paste("./Data/outpout.ET/resultado.global",essay,".csv",sep=""),
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

################3
plot.et.model  <- sapply (listapot, function (filtro) {        # aca reeplazaria por listaue
  datos <- filter (dt, pot == filtro) # aca reemplaze por dt ec2b.rep
  S     <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
  datos$water.SPM.t <-  datos$w.SPM - S 
  W0 <- datos$w.SPM[which(datos$time == 0)]
  datos$ET.t <- W0 - datos$w.SPM
  datos_ET <- datos 
  #W0 = datos$w.SPM[which(datos$time == 0)]
  Wn =  datos$w.SPM[which(datos$time == t1)]
  W2n = datos$w.SPM[which(datos$time == 2*t1)]
  n = t1
  Bs <-  (Wn - W0)^2 /(W2n - 2*Wn + W0)                     # este es el "B"  estimado por Simondi
  ks <-  -(1/n) * log ((W2n-Wn)/ (Wn-W0))                   # este es el "k"  estimado por Simondi
  if (is.na (ks) == FALSE) {
    if(ks > 0.05){ 
      resultados <- cbind(ks,Bs)
      
      eq.07 <- nls (ET.t ~ B * (1 - exp (-k*time)),               # esta es la estimacion de los parametros por Gauss-Newton
                    start = list (k = ks, B = Bs),
                    trace = FALSE , 
                    data= datos, nls.control(warnOnly=TRUE))
      
      k.eq.07 <- coef (eq.07)[1]                             # este es el k de Gauss-Newton
      B.eq.07 <- coef (eq.07)[2]                             # este es el B de Gauss-Newton
      tm <- (1/ k.eq.07) * log(2)
    
      resultados <- cbind(k.eq.07,B.eq.07)
      
      es   <-  data.frame(datos$essay[1])
      geno <-  data.frame(datos$genotype[1])
      pot  <-  filtro
      
      resultados <- data.frame (es,geno,pot, resultados)
      resultados <- resultados %>%
        rename (essay = datos.essay.1.)%>%
        rename (genotype = datos.genotype.1.) %>%
        mutate (tm = tm)
    

svg (filename= paste ("./Figures/Plots.ET/pot_",filtro,".svg", sep=""), 
           width=7, 
           height=5, 
          pointsize=12)   
    
       main <- paste ("G_",datos$genotype [1])
     
       plot (ET.t ~ time, pch=16,
            col="Black", 
            cex=1.5, 
            ylim=c(0,max(datos$ET.t)+10),
            axes=F,
            las=2,
            xlab="", ylab="", 
            type="p",
            main=paste (main,"_","pot_",filtro),
            data=datos_ET)
      
      points (ET.t ~ time , pch=16, type="p", 
             lty=3, lwd=1.5, cex=0.7,
             col="black", bg="black",
             data=datos)
      
      curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
              from=0, to = max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
              add = TRUE, col="red",lty= 2)
      
      curve  (Bs *(1- exp (- ks  * x)), 
              from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
              add = TRUE, col="black",lty= 1)
      
      abline (h=B.eq.07, lwd=1, lty =2, col="red")
      segments(x0=-5, y0=B.eq.07/2, x1 = tm, y1 =B.eq.07/2,col="red",lty= 2)
      segments(x0= tm, y0=-10, x1 = tm, y1 =  B.eq.07/2,col="red",lty= 2)
      
      axis(2,col="black",las=1)  ## las=1 makes horizontal labels
      mtext("ETt (g)",side=2,line=2.5)
      axis (1,0:max(dt$time))
      mtext("time (d)",side=1,line=2.5)
      legend ("topleft", legend =c("gauss.newton","simondi" )
              , col =c("red", "black" ), ncol = 1,lwd=2, pch=15, bty = "n")
      box()
      
      dev.off()
      
plot (ET.t ~ time, pch=16,
            col="Black", 
            cex=1.5, 
            ylim=c(0,max(datos_ET$ET.t)+10),
            axes=F,
            las=2,
            xlab="", ylab="", 
            type="p",
            main= paste ("pot_",filtro),
            data=datos_ET)
      
      points (ET.t ~ time , pch=16, type="p", 
              lty=3, lwd=1.5, cex=0.7,
              col="black", bg="black",
              data=datos)
      
      curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
              from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
              add = TRUE, col="red",lty= 2)
      
      curve  (Bs *(1- exp (- ks  * x)), 
              from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
              add = TRUE, col="black",lty= 1)
      
      abline (h=B.eq.07, lwd=1, lty =2, col="red")
      segments(x0=-5, y0=B.eq.07/2, x1 = tm, y1 =B.eq.07/2,col="red",lty= 2)
      segments(x0= tm, y0=-10, x1 = tm, y1 =  B.eq.07/2,col="red",lty= 2)
      
      axis(2,col="black",las=1)  ## las=1 makes horizontal labels
      mtext("ETt (g)",side=2,line=2.5)
      axis (1,0:max(dt$time))
      mtext("time (d)",side=1,line=2.5)
      legend ("topleft", legend =c("gauss.newton","simondi" )
              , col =c("red", "black" ), ncol = 1,lwd=2, pch=15, bty = "n")
      
      box()
     }
   }
 })
return (resultado.global)
} 




ET.ec2rep <- ET ( dt = ec2b.rep,
                  genotype = "genotype",t1 =4,
                  W="w.SPM",pss = 1970 , pm = 215, 
                  essay="ec2b.rep1",
                  ue ="pot")


ET.1 <- function (dt = NULL, t1=NULL, W=NULL, pss=NULL, pm = NULL, essay=NULL, ue=NULL) {
  dir.create (file.path ("Data"))
  dir.create (file.path ("Data", "outpout.ET"))
  
  dir.create (file.path ("Figures"))

  dir.create (file.path ("Figures", "Plots.ET"), showWarnings = FALSE)
  
  resultados <- NULL
  datos_ET <- NULL
  
########### si haces por pot #############################  
  if (ue == "pot") {
    listapot <- unique(dt$pot)
    
    dt.ET <- lapply (listapot, function(filtro){
      datos <- filter (dt, pot== filtro)
      S <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
      datos$water.SPM.t <-  datos$w.SPM -S 
      W0 = datos$w.SPM[which(datos$time == 0)]
      datos$ET.t <- W0 - datos$w.SPM
      
      write.table (datos, file =paste("./Data/outpout.ET/dt.ET","pot_",filtro,".txt",sep=""),
                   append = FALSE, quote = TRUE, sep = ",",
                   eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                   col.names = TRUE)
      datos_ET <- datos 
      print (datos_ET) 
    })
    
    datos.ET <- do.call(rbind.data.frame, dt.ET)
    
    write.table (datos.ET, file =paste("./Data/outpout.ET/dt.ET",essay,".csv",sep=""),
                 append = FALSE, quote = TRUE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE)
    
    
    et.model  <- sapply (listapot, function (filtro) {        # aca reeplazaria por listaue
      datos <- filter (dt, pot == filtro) # aca reemplaze por dt ec2b.rep
      S     <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
      datos$water.SPM.t <-  datos$w.SPM - S 
      W0 <- datos$w.SPM[which(datos$time == 0)]
      datos$ET.t <- W0 - datos$w.SPM
      
      #W0 = datos$w.SPM[which(datos$time == 0)]
      Wn =  datos$w.SPM[which(datos$time == t1)]
      W2n = datos$w.SPM[which(datos$time == 2*t1)]
      n = t1
      Bs <-  (Wn - W0)^2 /(W2n - 2*Wn + W0)                     # este es el "B"  estimado por Simondi
      ks <-  -(1/n) * log ((W2n-Wn)/ (Wn-W0))                   # este es el "k"  estimado por Simondi
      
      if (is.na (ks) == FALSE) {
        
        if(ks > 0.05){ 
          
          resultados <- cbind(ks,Bs)
          
          eq.07 <- nls (ET.t ~ B * (1 - exp (-k*time)),               # esta es la estimacion de los parametros por Gauss-Newton
                        start = list (k = ks, B = Bs),
                        trace = FALSE , 
                        data= datos, nls.control(warnOnly=TRUE))
          
          k.eq.07 <- coef (eq.07)[1]                             # este es el k de Gauss-Newton
          B.eq.07 <- coef (eq.07)[2]                             # este es el B de Gauss-Newton
          
          resultados <- cbind(B.eq.07,resultados)
          resultados <- cbind(k.eq.07,resultados)
          pot <-  filtro
          es  <-  data.frame(datos$essay[1])
          resultados <- cbind(pot,resultados)
          resultados <- cbind(es,resultados)
          
          write.table (resultados, file =paste("./par_ET_model","pot_",filtro,".txt",sep=""),
                       append = FALSE, quote = TRUE, sep = ",",
                       eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                       col.names = TRUE)
          print(resultados)
          
        }
      }
    })
    
    resultado.global <- do.call(rbind.data.frame, et.model)
    
    write.table (resultado.global, file =paste("./Data/outpout.ET/resultado.global.pot",essay,".csv",sep=""),
                 append = FALSE, quote = TRUE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE)
    
    ################3
    plot.et.model  <- sapply (listapot, function (filtro) {        # aca reeplazaria por listaue
      datos <- filter (dt, pot == filtro) # aca reemplaze por dt ec2b.rep
      S     <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
      datos$water.SPM.t <-  datos$w.SPM - S 
      W0 <- datos$w.SPM[which(datos$time == 0)]
      datos$ET.t <- W0 - datos$w.SPM
      datos_ET <- datos 
      #W0 = datos$w.SPM[which(datos$time == 0)]
      Wn =  datos$w.SPM[which(datos$time == t1)]
      W2n = datos$w.SPM[which(datos$time == 2*t1)]
      n = t1
      Bs <-  (Wn - W0)^2 /(W2n - 2*Wn + W0)                     # este es el "B"  estimado por Simondi
      ks <-  -(1/n) * log ((W2n-Wn)/ (Wn-W0))                   # este es el "k"  estimado por Simondi
      if (is.na (ks) == FALSE) {
        if(ks > 0.05){ 
          resultados <- cbind(ks,Bs)
          
          eq.07 <- nls (ET.t ~ B * (1 - exp (-k*time)),               # esta es la estimacion de los parametros por Gauss-Newton
                        start = list (k = ks, B = Bs),
                        trace = FALSE , 
                        data= datos, nls.control(warnOnly=TRUE))
          
          k.eq.07 <- coef (eq.07)[1]                             # este es el k de Gauss-Newton
          B.eq.07 <- coef (eq.07)[2]                             # este es el B de Gauss-Newton
          
          resultados <- cbind(B.eq.07,resultados)
          resultados <- cbind(k.eq.07,resultados)
          pot <-  filtro
          es  <-  data.frame(datos$essay[1])
          resultados <- cbind(pot,resultados)
          resultados <- cbind(es,resultados)
          
          svg (filename= paste ("./Figures/Plots.ET/pot_",filtro,".svg", sep=""), 
               width=7, 
               height=5, 
               pointsize=12)   
          
          plot (ET.t ~ time, pch=16,
                col="Black", 
                cex=1.5, 
                ylim=c(0,max(datos_ET$ET.t)+10),
                axes=F,
                las=2,
                xlab="", ylab="", 
                type="p",
                main= paste ("pot_",filtro),
                data=datos_ET)
          
          points (ET.t ~ time , pch=16, type="p", 
                  lty=3, lwd=1.5, cex=0.7,
                  col="black", bg="black",
                  data=datos)
          
          curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="red",lty= 2)
          
          curve  (Bs *(1- exp (- ks  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="black",lty= 1)
          
          axis(2,col="black",las=1)  ## las=1 makes horizontal labels
          mtext("ETt (g)",side=2,line=2.5)
          axis (1,0:max(dt$time))
          mtext("time (d)",side=1,line=2.5)
          legend ("topleft", legend =c("em.linea","pred.linea" )
                  , col =c("black", "blue" ), ncol = 1,lwd=2, pch=15, bty = "n")
          
          box()
          dev.off()
          
          plot (ET.t ~ time, pch=16,
                col="Black", 
                cex=1.5, 
                ylim=c(0,max(datos_ET$ET.t)+10),
                axes=F,
                las=2,
                xlab="", ylab="", 
                type="p",
                main= paste ("pot_",filtro),
                data=datos_ET)
          
          points (ET.t ~ time , pch=16, type="p", 
                  lty=3, lwd=1.5, cex=0.7,
                  col="black", bg="black",
                  data=datos)
          
          curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="red",lty= 2)
          
          curve  (Bs *(1- exp (- ks  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="black",lty= 1)
          
          axis(2,col="black",las=1)  ## las=1 makes horizontal labels
          mtext("ETt (g)",side=2,line=2.5)
          axis (1,0:max(dt$time))
          mtext("time (d)",side=1,line=2.5)
          legend ("topleft", legend =c("em.linea","pred.linea" )
                  , col =c("black", "blue" ), ncol = 1,lwd=2, pch=15, bty = "n")
          
          box()
          
        }
      }
    })
  } 
  
  if (ue == "genotype") {
    listageno <- unique(dt$genotype)
    
    dt.ET <- lapply (listageno, function (filtro){
      datos <- filter (dt, genotype == filtro)
      S <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
      datos$water.SPM.t <-  datos$w.SPM -S 
      W0 = datos$w.SPM[which(datos$time == 0)]
      datos$ET.t <- W0 - datos$w.SPM
      
      write.table (datos, file =paste("./Data/outpout.ET/dt.ET","geno_",filtro,".txt",sep=""),
                   append = FALSE, quote = TRUE, sep = ",",
                   eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                   col.names = TRUE)
      datos_ET <- datos 
      print (datos_ET) 
    })
    
    datos.ET <- do.call(rbind.data.frame, dt.ET)
    
    write.table (datos.ET, file =paste("./Data/outpout.ET/dt.ET.geno.",essay,".csv",sep=""),
                 append = FALSE, quote = TRUE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE)
    
    
    et.model  <- sapply (listageno, function (filtro) {        # aca reeplazaria por listaue
      datos <- filter (dt, genotype == filtro) # aca reemplaze por dt ec2b.rep
      S     <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
      datos$water.SPM.t <-  datos$w.SPM - S 
      W0 <- datos$w.SPM[which(datos$time == 0)]
      datos$ET.t <- W0 - datos$w.SPM
      
      #W0 = datos$w.SPM[which(datos$time == 0)]
      Wn =  datos$w.SPM[which(datos$time == t1)]
      W2n = datos$w.SPM[which(datos$time == 2*t1)]
      n = t1
      Bs <-  (Wn - W0)^2 /(W2n - 2*Wn + W0)                     # este es el "B"  estimado por Simondi
      ks <-  -(1/n) * log ((W2n-Wn)/ (Wn-W0))                   # este es el "k"  estimado por Simondi
      
      if (is.na (ks) == FALSE) {
        
        if(ks > 0.05){ 
          
          resultados <- cbind(ks,Bs)
          
          eq.07 <- nls (ET.t ~ B * (1 - exp (-k*time)),               # esta es la estimacion de los parametros por Gauss-Newton
                        start = list (k = ks, B = Bs),
                        trace = FALSE , 
                        data= datos, nls.control(warnOnly=TRUE))
          
          k.eq.07 <- coef (eq.07)[1]                             # este es el k de Gauss-Newton
          B.eq.07 <- coef (eq.07)[2]                             # este es el B de Gauss-Newton
          
          resultados <- cbind(B.eq.07,resultados)
          resultados <- cbind(k.eq.07,resultados)
          genotype <-  filtro
          es  <-  data.frame(datos$essay[1])
          resultados <- cbind(genotype,resultados)
          resultados <- cbind(es,resultados)
          
          write.table (resultados, file =paste("./par_ET_model","geno_",filtro,".txt",sep=""),
                       append = FALSE, quote = TRUE, sep = ",",
                       eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                       col.names = TRUE)
          print(resultados)
          
        }
      }
    })
    
    resultado.global <- do.call(rbind.data.frame, et.model)
    
    write.table (resultado.global, file =paste("./Data/outpout.ET/resultado.global.geno",essay,".csv",sep=""),
                 append = FALSE, quote = TRUE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE)
    
    ################3
    plot.et.model  <- sapply (listageno, function (filtro) {        # aca reeplazaria por listaue
      datos <- filter (dt, genotype == filtro) # aca reemplaze por dt ec2b.rep
      S     <- pss + pm # Este es el peso del sustrato mas el peso de la maceta
      datos$water.SPM.t <-  datos$w.SPM - S 
      W0 <- datos$w.SPM[which(datos$time == 0)]
      datos$ET.t <- W0 - datos$w.SPM
      datos_ET <- datos 
      #W0 = datos$w.SPM[which(datos$time == 0)]
      Wn =  datos$w.SPM[which(datos$time == t1)]
      W2n = datos$w.SPM[which(datos$time == 2*t1)]
      n = t1
      Bs <-  (Wn - W0)^2 /(W2n - 2*Wn + W0)                     # este es el "B"  estimado por Simondi
      ks <-  -(1/n) * log ((W2n-Wn)/ (Wn-W0))                   # este es el "k"  estimado por Simondi
      if (is.na (ks) == FALSE) {
        if(ks > 0.05){ 
          resultados <- cbind(ks,Bs)
          
          eq.07 <- nls (ET.t ~ B * (1 - exp (-k*time)),               # esta es la estimacion de los parametros por Gauss-Newton
                        start = list (k = ks, B = Bs),
                        trace = FALSE , 
                        data= datos, nls.control(warnOnly=TRUE))
          
          k.eq.07 <- coef (eq.07)[1]                             # este es el k de Gauss-Newton
          B.eq.07 <- coef (eq.07)[2]                             # este es el B de Gauss-Newton
          
          resultados <- cbind(B.eq.07,resultados)
          resultados <- cbind(k.eq.07,resultados)
          geno <-  filtro
          es  <-  data.frame(datos$essay[1])
          resultados <- cbind(geno,resultados)
          resultados <- cbind(es,resultados)
          
          svg (filename= paste ("./Figures/Plots.ET/geno_",filtro,".svg", sep=""), 
               width=7, 
               height=5, 
               pointsize=12)   
          
          plot (ET.t ~ time, pch=16,
                col="Black", 
                cex=1.5, 
                ylim=c(0,max(datos_ET$ET.t)+10),
                axes=F,
                las=2,
                xlab="", ylab="", 
                type="p",
                main= paste ("geno_",filtro),
                data=datos_ET)
          
          points (ET.t ~ time , pch=16, type="p", 
                  lty=3, lwd=1.5, cex=0.7,
                  col="black", bg="black",
                  data=datos)
          
          curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="red",lty= 2)
          
          curve  (Bs *(1- exp (- ks  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="black",lty= 1)
          
          axis(2,col="black",las=1)  ## las=1 makes horizontal labels
          mtext("ETt (g)",side=2,line=2.5)
          axis (1,0:max(dt$time))
          mtext("time (d)",side=1,line=2.5)
          legend ("topleft", legend =c("em.linea","pred.linea" )
                  , col =c("black", "blue" ), ncol = 1,lwd=2, pch=15, bty = "n")
          
          box()
          dev.off()
          
          plot (ET.t ~ time, pch=16,
                col="Black", 
                cex=1.5, 
                ylim=c(0,max(datos_ET$ET.t)+10),
                axes=F,
                las=2,
                xlab="", ylab="", 
                type="p",
                main= paste ("geno_",filtro),
                data=datos_ET)
          
          points (ET.t ~ time , pch=16, type="p", 
                  lty=3, lwd=1.5, cex=0.7,
                  col="black", bg="black",
                  data=datos)
          
          curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="red",lty= 2)
          
          curve  (Bs *(1- exp (- ks  * x)), 
                  from=0, to=max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
                  add = TRUE, col="black",lty= 1)
          
          axis(2,col="black",las=1)  ## las=1 makes horizontal labels
          mtext("ETt (g)",side=2,line=2.5)
          axis (1,0:max(dt$time))
          mtext("time (d)",side=1,line=2.5)
          legend ("topleft", legend =c("em.linea","pred.linea" )
                  , col =c("black", "blue" ), ncol = 1,lwd=2, pch=15, bty = "n")
          
          box()
          
        }
      }
    })
  } 
  
  
  
  
  
  
  
}

ET (dt = ec2b.rep, pss = 1970, pm = 215, t1 =4 , W="w.SPM",essay="ec2b.rep1", ue ="pot")







#if (ue =="geno"){

#}


  #ET.model(dt = ec2b.rep1.1, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2b.rep1", ue="pot") )

ec2b.rep1.2 <- filter (ec2b.rep1, pot=="2")
ec2b.rep1.3 <- filter (ec2b.rep1, pot=="3")
ec2b.rep1.4 <- filter (ec2b.rep1, pot=="4")
ec2b.rep1.5 <- filter (ec2b.rep1, pot=="5")
ec2b.rep1.6 <- filter (ec2b.rep1, pot=="6")
#ec2b.rep1.7 <- filter (ec2b.rep1, pot=="7")
ec2b.rep1.8 <- filter (ec2b.rep1, pot=="8")
ec2b.rep1.9 <- filter (ec2b.rep1, pot=="9")
ec2b.rep1.10 <- filter (ec2b.rep1, pot=="10")
ec2b.rep1.11 <- filter (ec2b.rep1, pot=="11")
ec2b.rep1.12 <- filter (ec2b.rep1, pot=="12")
ec2b.rep1.13 <- filter (ec2b.rep1, pot=="13")
ec2b.rep1.14 <- filter (ec2b.rep1, pot=="14")
ec2b.rep1.15 <- filter (ec2b.rep1, pot=="15")
ec2b.rep1.16 <- filter (ec2b.rep1, pot=="16")
ec2b.rep1.17 <- filter (ec2b.rep1, pot=="17")
ec2b.rep1.18 <- filter (ec2b.rep1, pot=="18")
ec2b.rep1.19 <- filter (ec2b.rep1, pot=="19")
ec2b.rep1.20 <- filter (ec2b.rep1, pot=="20")
ec2b.rep1.21 <- filter (ec2b.rep1, pot=="21")
ec2b.rep1.22 <- filter (ec2b.rep1, pot=="22")
ec2b.rep1.23 <- filter (ec2b.rep1, pot=="23")
ec2b.rep1.24 <- filter (ec2b.rep1, pot=="24")
ec2b.rep1.25 <- filter (ec2b.rep1, pot=="25")
ec2b.rep1.26 <- filter (ec2b.rep1, pot=="26")
ec2b.rep1.27 <- filter (ec2b.rep1, pot=="27")
ec2b.rep1.28 <- filter (ec2b.rep1, pot=="28")
ec2b.rep1.29 <- filter (ec2b.rep1, pot=="29")
ec2b.rep1.30 <- filter (ec2b.rep1, pot=="30")
ec2b.rep1.31 <- filter (ec2b.rep1, pot=="31")
ec2b.rep1.32 <- filter (ec2b.rep1, pot=="32")
ec2b.rep1.33 <- filter (ec2b.rep1, pot=="33")
ec2b.rep1.34 <- filter (ec2b.rep1, pot=="34")
ec2b.rep1.35 <- filter (ec2b.rep1, pot=="35")
ec2b.rep1.36 <- filter (ec2b.rep1, pot=="36")
ec2b.rep1.37 <- filter (ec2b.rep1, pot=="37")
ec2b.rep1.38 <- filter (ec2b.rep1, pot=="38")
ec2b.rep1.39 <- filter (ec2b.rep1, pot=="39")
ec2b.rep1.40 <- filter (ec2b.rep1, pot=="40")
ec2b.rep1.41 <- filter (ec2b.rep1, pot=="41")
ec2b.rep1.42 <- filter (ec2b.rep1, pot=="42")
ec2b.rep1.43 <- filter (ec2b.rep1, pot=="43")
ec2b.rep1.44 <- filter (ec2b.rep1, pot=="44")
ec2b.rep1.45 <- filter (ec2b.rep1, pot=="45")
ec2b.rep1.46 <- filter (ec2b.rep1, pot=="46")

#ec2b.rep1.47 <- filter (ec2b.rep1, pot=="47")
#ec2b.rep1.48 <- filter (ec2b.rep1, pot=="48")
ec2b.rep1.49 <- filter (ec2b.rep1, pot=="49")
ec2b.rep1.50 <- filter (ec2b.rep1, pot=="50")
ec2b.rep1.51 <- filter (ec2b.rep1, pot=="51")
ec2b.rep1.52 <- filter (ec2b.rep1, pot=="52")
ec2b.rep1.53 <- filter (ec2b.rep1, pot=="53")
ec2b.rep1.54 <- filter (ec2b.rep1, pot=="54")
ec2b.rep1.55 <- filter (ec2b.rep1, pot=="55")
ec2b.rep1.56 <- filter (ec2b.rep1, pot=="56")
ec2b.rep1.57 <- filter (ec2b.rep1, pot=="57")
ec2b.rep1.58 <- filter (ec2b.rep1, pot=="58")
ec2b.rep1.59 <- filter (ec2b.rep1, pot=="59")
ec2b.rep1.60 <- filter (ec2b.rep1, pot=="60")
ec2b.rep1.61 <- filter (ec2b.rep1, pot=="61")
ec2b.rep1.62 <- filter (ec2b.rep1, pot=="62")
ec2b.rep1.63 <- filter (ec2b.rep1, pot=="63")
ec2b.rep1.64 <- filter (ec2b.rep1, pot=="64")
ec2b.rep1.65 <- filter (ec2b.rep1, pot=="65")
ec2b.rep1.66 <- filter (ec2b.rep1, pot=="66")
ec2b.rep1.67 <- filter (ec2b.rep1, pot=="67")
ec2b.rep1.68 <- filter (ec2b.rep1, pot=="68")
ec2b.rep1.69 <- filter (ec2b.rep1, pot=="69")
ec2b.rep1.70 <- filter (ec2b.rep1, pot=="70")
ec2b.rep1.71 <- filter (ec2b.rep1, pot=="71")
#ec2b.rep1.72 <- filter (ec2b.rep1, pot=="72")
ec2b.rep1.73 <- filter (ec2b.rep1, pot=="73")
ec2b.rep1.74 <- filter (ec2b.rep1, pot=="74")
ec2b.rep1.75 <- filter (ec2b.rep1, pot=="75")
ec2b.rep1.76 <- filter (ec2b.rep1, pot=="76")
ec2b.rep1.77 <- filter (ec2b.rep1, pot=="77")
ec2b.rep1.78 <- filter (ec2b.rep1, pot=="78")
ec2b.rep1.79 <- filter (ec2b.rep1, pot=="79")
ec2b.rep1.80 <- filter (ec2b.rep1, pot=="80")
ec2b.rep1.81 <- filter (ec2b.rep1, pot=="81")
ec2b.rep1.82 <- filter (ec2b.rep1, pot=="82")
ec2b.rep1.83 <- filter (ec2b.rep1, pot=="83")
ec2b.rep1.84 <- filter (ec2b.rep1, pot=="84")
ec2b.rep1.85 <- filter (ec2b.rep1, pot=="85")
ec2b.rep1.86 <- filter (ec2b.rep1, pot=="86")
ec2b.rep1.87 <- filter (ec2b.rep1, pot=="87")
ec2b.rep1.88 <- filter (ec2b.rep1, pot=="88")
ec2b.rep1.89 <- filter (ec2b.rep1, pot=="89")
ec2b.rep1.90 <- filter (ec2b.rep1, pot=="90")
ec2b.rep1.91 <- filter (ec2b.rep1, pot=="91")
ec2b.rep1.92 <- filter (ec2b.rep1, pot=="92")
ec2b.rep1.93 <- filter (ec2b.rep1, pot=="93")
#ec2b.rep1.94 <- filter (ec2b.rep1, pot=="94")
ec2b.rep1.95 <- filter (ec2b.rep1, pot=="95")
ec2b.rep1.96 <- filter (ec2b.rep1, pot=="96")
ec2b.rep1.97 <- filter (ec2b.rep1, pot=="97")
ec2b.rep1.98 <- filter (ec2b.rep1, pot=="98")
ec2b.rep1.99 <- filter (ec2b.rep1, pot=="99")
ec2b.rep1.100 <- filter (ec2b.rep1, pot=="100")
ec2b.rep1.101 <- filter (ec2b.rep1, pot=="101")
ec2b.rep1.102 <- filter (ec2b.rep1, pot=="102")
ec2b.rep1.103 <- filter (ec2b.rep1, pot=="103")
ec2b.rep1.104 <- filter (ec2b.rep1, pot=="104")

pb5.1 <- ET.model (dt = ec2b.rep1.1, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2b.rep1", ue="pot")
pb5.2 <- ET.model (dt = ec2b.rep1.2, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.3 <- ET.model (dt = ec2b.rep1.3, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.4 <- ET.model (dt = ec2b.rep1.4, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.5 <- ET.model (dt = ec2b.rep1.5, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.6 <- ET.model (dt = ec2b.rep1.6, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.7 <- ET.model (dt = ec2b.rep1.7, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.8 <- ET.model (dt = ec2b.rep1.8, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.9 <- ET.model (dt = ec2b.rep1.9, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.10 <- ET.model (dt = ec2b.rep1.10, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.11 <- ET.model (dt = ec2b.rep1.11, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.12 <- ET.model (dt = ec2b.rep1.12, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.13 <- ET.model (dt = ec2b.rep1.13, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.14 <- ET.model (dt = ec2b.rep1.14, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.15 <- ET.model (dt = ec2b.rep1.15, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.16 <- ET.model (dt = ec2b.rep1.16, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.17 <- ET.model (dt = ec2b.rep1.17, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.18 <- ET.model (dt = ec2b.rep1.18, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.19 <- ET.model (dt = ec2b.rep1.19, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.20 <- ET.model (dt = ec2b.rep1.20, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.21 <- ET.model (dt = ec2b.rep1.21, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.22 <- ET.model (dt = ec2b.rep1.22, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.23 <- ET.model (dt = ec2b.rep1.23, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.24 <- ET.model (dt = ec2b.rep1.24, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.25 <- ET.model (dt = ec2b.rep1.25, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.26 <- ET.model (dt = ec2b.rep1.26, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.27 <- ET.model (dt = ec2b.rep1.27, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.28 <- ET.model (dt = ec2b.rep1.28, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.29 <- ET.model (dt = ec2b.rep1.29, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.30 <- ET.model (dt = ec2b.rep1.30, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.31 <- ET.model (dt = ec2b.rep1.31, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.32 <- ET.model (dt = ec2b.rep1.32, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.33 <- ET.model (dt = ec2b.rep1.33, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.34 <- ET.model (dt = ec2b.rep1.34, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.35 <- ET.model (dt = ec2b.rep1.35, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.36 <- ET.model (dt = ec2b.rep1.36, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.37 <- ET.model (dt = ec2b.rep1.37, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.38 <- ET.model (dt = ec2b.rep1.38, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.39 <- ET.model (dt = ec2b.rep1.39, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.40 <- ET.model (dt = ec2b.rep1.40, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.41 <- ET.model (dt = ec2b.rep1.41, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.42 <- ET.model (dt = ec2b.rep1.42, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.43 <- ET.model (dt = ec2b.rep1.43, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.44 <- ET.model (dt = ec2b.rep1.44, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.45 <- ET.model (dt = ec2b.rep1.45, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.46 <- ET.model (dt = ec2b.rep1.46, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.47 <- ET.model (dt = ec2b.rep1.47, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.48 <- ET.model (dt = ec2b.rep1.48, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.49 <- ET.model (dt = ec2b.rep1.49, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.50 <- ET.model (dt = ec2b.rep1.50, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.51 <- ET.model (dt = ec2b.rep1.51, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.52 <- ET.model (dt = ec2b.rep1.52, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.53 <- ET.model (dt = ec2b.rep1.53, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.54 <- ET.model (dt = ec2b.rep1.54, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.55 <- ET.model (dt = ec2b.rep1.55, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.56 <- ET.model (dt = ec2b.rep1.56, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.57 <- ET.model (dt = ec2b.rep1.57, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.58 <- ET.model (dt = ec2b.rep1.58, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.59 <- ET.model (dt = ec2b.rep1.59, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.60 <- ET.model (dt = ec2b.rep1.60, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.61 <- ET.model (dt = ec2b.rep1.61, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.62 <- ET.model (dt = ec2b.rep1.62, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.63 <- ET.model (dt = ec2b.rep1.63, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.64 <- ET.model (dt = ec2b.rep1.64, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.65 <- ET.model (dt = ec2b.rep1.65, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.66 <- ET.model (dt = ec2b.rep1.66, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.67 <- ET.model (dt = ec2b.rep1.67, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.68 <- ET.model (dt = ec2b.rep1.68, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.69 <- ET.model (dt = ec2b.rep1.69, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.70 <- ET.model (dt = ec2b.rep1.70, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.71 <- ET.model (dt = ec2b.rep1.71, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.72 <- ET.model (dt = ec2b.rep1.72, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.73 <- ET.model (dt = ec2b.rep1.73, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.74 <- ET.model (dt = ec2b.rep1.74, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.75 <- ET.model (dt = ec2b.rep1.75, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.76 <- ET.model (dt = ec2b.rep1.76, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.77 <- ET.model (dt = ec2b.rep1.77, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.78 <- ET.model (dt = ec2b.rep1.78, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.79 <- ET.model (dt = ec2b.rep1.79, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.80 <- ET.model (dt = ec2b.rep1.80, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.81 <- ET.model (dt = ec2b.rep1.81, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.82 <- ET.model (dt = ec2b.rep1.82, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.83 <- ET.model (dt = ec2b.rep1.83, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.84 <- ET.model (dt = ec2b.rep1.84, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.85 <- ET.model (dt = ec2b.rep1.85, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.86 <- ET.model (dt = ec2b.rep1.86, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.87 <- ET.model (dt = ec2b.rep1.87, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.88 <- ET.model (dt = ec2b.rep1.88, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.89 <- ET.model (dt = ec2b.rep1.89, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.90 <- ET.model (dt = ec2b.rep1.90, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.91 <- ET.model (dt = ec2b.rep1.91, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.92 <- ET.model (dt = ec2b.rep1.92, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.93 <- ET.model (dt = ec2b.rep1.93, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
#pb5.94 <- ET.model (dt = ec2b.rep1.94, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.95 <- ET.model (dt = ec2b.rep1.95, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.96 <- ET.model (dt = ec2b.rep1.96, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.97 <- ET.model (dt = ec2b.rep1.97, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.98 <- ET.model (dt = ec2b.rep1.98, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.99 <- ET.model (dt = ec2b.rep1.99, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.100 <- ET.model (dt = ec2b.rep1.100, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.101 <- ET.model (dt = ec2b.rep1.101, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.102 <- ET.model (dt = ec2b.rep1.102, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.103 <- ET.model (dt = ec2b.rep1.103, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")
pb5.104 <- ET.model (dt = ec2b.rep1.104, t1=4, W="w.SPM", pss=1970, pm = 215,  essay="ec2b.rep1", ue="pot")

ec2b.rep1.ETmodel.par <- rbind (pb5.1,pb5.2
                           #,pb5.3,pb5.4
                           ,pb5.5,pb5.6
                           #,pb5.7
                           ,pb5.8,pb5.9,pb5.10
                           ,pb5.11,pb5.12,pb5.13,pb5.14,pb5.15,pb5.16,pb5.17,pb5.18,pb5.19,pb5.20
                           ,pb5.21,pb5.22,pb5.23,pb5.24,pb5.25,pb5.26,pb5.27,pb5.28,pb5.29,pb5.30
                           ,pb5.31,pb5.32,pb5.33,pb5.34,pb5.35,pb5.36,pb5.37,pb5.38,pb5.39,pb5.40
                           ,pb5.41,pb5.42,pb5.43,pb5.44,pb5.45,pb5.46
                           #,pb5.47,pb5.48
                           ,pb5.49
                           ,pb5.50,pb5.51,pb5.52,pb5.53,pb5.54,pb5.55,pb5.56,pb5.57,pb5.58,pb5.59
                           ,pb5.60,pb5.61,pb5.62,pb5.63,pb5.64,pb5.65,pb5.66,pb5.67,pb5.68,pb5.69
                           ,pb5.70,pb5.71
                           #,pb5.72
                           ,pb5.73,pb5.74,pb5.75,pb5.76,pb5.77,pb5.78,pb5.79
                           ,pb5.80,pb5.81,pb5.82,pb5.83,pb5.84,pb5.85,pb5.86,pb5.87,pb5.88,pb5.89
                           ,pb5.90,pb5.91,pb5.92,pb5.93
                           #,pb5.94
                           ,pb5.95,pb5.96,pb5.97,pb5.98,pb5.99
                           ,pb5.100,pb5.101,pb5.102,pb5.103,pb5.104)


length(unique(ec2b.rep1.ETmodel.par$pot))

dim (ec2b.rep1.ETmodel.par)

ec2b.rep2 <- filter (ec2b.rep1, pot!="3",pot!="4",time==0) 

ec2b.rep2a <- ec2b.rep2 %>%
              select(essay:trat )

unique(ec2b.rep2a$pot) == unique(ec2b.rep1.ETmodel.par$pot)
ec2b.rep2a$pot == ec2b.rep1.ETmodel.par$pot


ec2b.rep1a.ETmodel.par <- ec2b.rep1.ETmodel.par %>%
                         select (-c(essay:pot))
################3 esta es la metris de Ec2b.rep para juntar con la otra matriz 

ec2b.rep3 <- data.frame (ec2b.rep2a, ec2b.rep1a.ETmodel.par)


#### ec2d
ec2d <- read.table ("./Data/rawdata/input/ec2d.txt" ,
                    header = TRUE, sep = "",dec = ".",
                    na.strings = "NA" )

###### control de la matriz #########3
#ec2d$pot <- as.factor(ec2d$pot)

length(unique(ec2d$pot))
length(unique(ec2d$genotype))

(chks.ec2d <-  table (ec2d$checks))

#### filtro los NA ##########
ec2d1 <- ec2d %>%
  filter (pot!="30" & pot!="31")

length(unique(ec2d1$pot))
length(unique(ec2d1$genotype))

(chks.ec2d1 <-  table (ec2d1$checks))

######## filtro para la ec2d
ec2d1.1 <- filter (ec2d1, pot=="1")
ec2d1.2 <- filter (ec2d1, pot=="2")
ec2d1.3 <- filter (ec2d1, pot=="3")
ec2d1.4 <- filter (ec2d1, pot=="4")
ec2d1.5 <- filter (ec2d1, pot=="5")
ec2d1.6 <- filter (ec2d1, pot=="6")
ec2d1.7 <- filter (ec2d1, pot=="7")
ec2d1.8 <- filter (ec2d1, pot=="8")
ec2d1.9 <- filter (ec2d1, pot=="9")
ec2d1.10 <- filter (ec2d1, pot=="10")
ec2d1.11 <- filter (ec2d1, pot=="11")
#ec2d1.12 <- filter (ec2d1, pot=="12")
ec2d1.13 <- filter (ec2d1, pot=="13")
ec2d1.14 <- filter (ec2d1, pot=="14")
ec2d1.15 <- filter (ec2d1, pot=="15")
ec2d1.16 <- filter (ec2d1, pot=="16")
ec2d1.17 <- filter (ec2d1, pot=="17")
ec2d1.18 <- filter (ec2d1, pot=="18")
ec2d1.19 <- filter (ec2d1, pot=="19")
ec2d1.20 <- filter (ec2d1, pot=="20")
ec2d1.21 <- filter (ec2d1, pot=="21")
ec2d1.22 <- filter (ec2d1, pot=="22")
ec2d1.23 <- filter (ec2d1, pot=="23")
ec2d1.24 <- filter (ec2d1, pot=="24")
ec2d1.25 <- filter (ec2d1, pot=="25")
ec2d1.26 <- filter (ec2d1, pot=="26")
ec2d1.27 <- filter (ec2d1, pot=="27")
ec2d1.28 <- filter (ec2d1, pot=="28")
ec2d1.29 <- filter (ec2d1, pot=="29")
#ec2d1.30 <- filter (ec2d1, pot=="30")
#ec2d1.31 <- filter (ec2d1, pot=="31")
ec2d1.32 <- filter (ec2d1, pot=="32")
ec2d1.33 <- filter (ec2d1, pot=="33")
ec2d1.34 <- filter (ec2d1, pot=="34")
ec2d1.35 <- filter (ec2d1, pot=="35")
ec2d1.36 <- filter (ec2d1, pot=="36")
ec2d1.37 <- filter (ec2d1, pot=="37")
ec2d1.38 <- filter (ec2d1, pot=="38")
ec2d1.39 <- filter (ec2d1, pot=="39")
ec2d1.40 <- filter (ec2d1, pot=="40")

pd1 <- ET.model (dt = ec2d1.1, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd2 <- ET.model (dt = ec2d1.2, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd3 <- ET.model (dt = ec2d1.3, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd4 <- ET.model (dt = ec2d1.4, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd5 <- ET.model (dt = ec2d1.5, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd6 <- ET.model (dt = ec2d1.6, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd7 <- ET.model (dt = ec2d1.7, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd8 <- ET.model (dt = ec2d1.8, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd9 <- ET.model (dt = ec2d1.9, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd10 <- ET.model (dt = ec2d1.10, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd11 <- ET.model (dt = ec2d1.11, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
#pd12 <- ET.model (dt = ec2d1.12, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd13 <- ET.model (dt = ec2d1.13, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd14 <- ET.model (dt = ec2d1.14, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd15 <- ET.model (dt = ec2d1.15, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd16 <- ET.model (dt = ec2d1.16, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd17 <- ET.model (dt = ec2d1.17, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd18 <- ET.model (dt = ec2d1.18, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd19 <- ET.model (dt = ec2d1.19, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd20 <- ET.model (dt = ec2d1.20, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd21 <- ET.model (dt = ec2d1.21, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd22 <- ET.model (dt = ec2d1.22, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd23 <- ET.model (dt = ec2d1.23, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd24 <- ET.model (dt = ec2d1.24, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd25 <- ET.model (dt = ec2d1.25, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd26 <- ET.model (dt = ec2d1.26, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd27 <- ET.model (dt = ec2d1.27, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd28 <- ET.model (dt = ec2d1.28, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd29 <- ET.model (dt = ec2d1.29, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
#pd30 <- ET.model (dt = ec2d1.30, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
#pd31 <- ET.model (dt = ec2d1.31, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd32 <- ET.model (dt = ec2d1.32, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd33 <- ET.model (dt = ec2d1.33, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd34 <- ET.model (dt = ec2d1.34, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd35 <- ET.model (dt = ec2d1.35, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd36 <- ET.model (dt = ec2d1.36, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd37 <- ET.model (dt = ec2d1.37, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd38 <- ET.model (dt = ec2d1.38, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd39 <- ET.model (dt = ec2d1.39, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")
pd40 <- ET.model (dt = ec2d1.40, t1=4, W="w.SPM", pss=1970, pm = 215, essay="ec2d1", ue="pot")

ec2d1.ETmodel.par <- rbind (pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10
                           ,pd11
                           #,pd12
                           ,pd13,pd14,pd15,pd16,pd17,pd18,pd19,pd20
                           ,pd21,pd22,pd23,pd24,pd25,pd26,pd27,pd28,pd29
                           #,pd30,pd31
                           ,pd32,pd33,pd34,pd35,pd36,pd37,pd38,pd39,pd40)

length(unique(ec2d1.ETmodel.par$pot))

dim (ec2d1.ETmodel.par)

ec2d.2 <- filter (ec2d1,time==0) 

ec2d.2a <- ec2d.2 %>%
  select(essay:trat )

unique(ec2d.2a$pot) == unique(ec2d1.ETmodel.par$pot)
        ec2d.2a$pot == ec2d1.ETmodel.par$pot


ec2d1a.ETmodel.par <- ec2d1.ETmodel.par %>%
                      select (-c(essay:pot))

################ esta es la matriz del ensayo Ec2b.rep para juntar con la otra matriz 

ec2d3 <- data.frame (ec2d.2a, ec2d1a.ETmodel.par)


############## 
##### armo la matriz para modelar ###########
dim (ec2b.rep3)
dim (ec2d3)

# controlo la cantidad de repeticiones 
head(ec2b.rep3)
(ec2b.rep3_teslines.1 <-  table (ec2b.rep3$testlines.1))
(ec2b.rep3_checks <-  table (ec2b.rep3$checks))


### ## esta es la matriz para modelar 
model_GWAS_EC_1 <- rbind (ec2b.rep3, ec2d3)

dim (model_GWAS_EC)

(model_GWAS_EC_check <-  table (model_GWAS_EC$checks))

write.table (model_GWAS_EC_1, file = "./Data/procdata/model_GWAS_EC_1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

######### Cargo la matrix que tiene los gw corregidos ########
ec_pheno.7.1 <- read.table ("./Data/rawdata/input/ec_pheno.7.1.txt" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )
summary (ec_pheno.7.1)


str (ec_pheno.7.1)
#ec_pheno.7.1$pot <- as.factor(ec_pheno.7.1$pot)

################
ec_pheno.7.1_b <- ec_pheno.7.1 %>%
                   filter ( essay == "ec2b.rep")


setdiff (ec2b.rep1$pot, ec_pheno.7.1_b$pot)

ec2b.rep4 <- filter (ec2b.rep1, 
                     pot!=3,
                     pot!=4,
                     pot!=12,
                     pot!=13,
                     pot!=16,
                     pot!=17,
                     pot!=22,
                     pot!=26,
                     pot!=34,
                     pot!=50,
                     pot!=59,
                     pot!=69,
                     pot!=71,
                     pot!=75,
                     pot!=101) 

setdiff (ec2b.rep4$pot, ec_pheno.7.1_b$pot)
setdiff (ec_pheno.7.1_b$pot, ec2b.rep4$pot)

# verifico que los Gw sean iguales 
ec2b.rep4.to <- ec2b.rep4 %>%
                filter (time ==0 )

ec2b.rep4.t4 <- ec2b.rep4 %>%
  filter (time ==4 )

ec2b.rep4.t7 <- ec2b.rep4 %>%
  filter (time ==7 )

ec2b.rep4.t11 <- ec2b.rep4 %>%
  filter (time ==11 )

plot (ec2b.rep4.to$Gw, ec_pheno.7.1_b$gw0)
plot (ec2b.rep4.t4$Gw, ec_pheno.7.1_b$gw4)
plot (ec2b.rep4.t7$Gw, ec_pheno.7.1_b$gw7)
plot (ec2b.rep4.t11$Gw, ec_pheno.7.1_b$gw11)

### acomodo los gw.cal de ec.pheno 
head( ec_pheno.7.1_b)

gw0cal.ec_pheno.7.1_b <- ec_pheno.7.1_b%>%
                        select (essay : linea.disenio.1,gwcal.0 )%>%
                        rename ( Gwcal = gwcal.0 ) %>%
                        mutate (time =0)


gwc4cal.ec_pheno.7.1_b <- ec_pheno.7.1_b%>%
                        select (essay : linea.disenio.1,gwcal.4 )%>%
                        rename (Gwcal = gwcal.4)%>%
                        mutate (time =4)

gwc7cal.ec_pheno.7.1_b <- ec_pheno.7.1_b%>%
                        select (essay : linea.disenio.1,gwcal.7 )%>%
                        rename (Gwcal = gwcal.7 )%>%
                        mutate (time =7)

gwc8cal.ec_pheno.7.1_b <- ec_pheno.7.1_b%>%
                       select (essay : linea.disenio.1)%>%
                       mutate (Gwcal = NA )%>%
                       mutate (time =8)


gwc11cal.ec_pheno.7.1_b <- ec_pheno.7.1_b%>%
                           select (essay : linea.disenio.1,gwcal.11 )%>%
                           rename (Gwcal = gwcal.11)%>%
                           mutate (time =11)


gwc12cal.ec_pheno.7.1_b <- ec_pheno.7.1_b%>%
                           select (essay : linea.disenio.1 )%>%
                           mutate (Gwcal = NA)%>%
                           mutate (time =12)



gwcal.b <- rbind  (gw0cal.ec_pheno.7.1_b,
                   gwc4cal.ec_pheno.7.1_b,
                   gwc7cal.ec_pheno.7.1_b,
                   gwc8cal.ec_pheno.7.1_b,
                   gwc11cal.ec_pheno.7.1_b,
                   gwc12cal.ec_pheno.7.1_b)

gwcal.b <- gwcal.b %>%
           arrange ( pot, time) 

summary(gwcal.b )
str(gwcal.b)

ec2b.rep4$pot <- as.numeric(ec2b.rep4$pot)

ec2b.rep4 <-  ec2b.rep4 %>%
              arrange ( pot, time)


(ec2b.rep4$pot) == (gwcal.b$pot)

ec2b.rep5 <- ec2b.rep4 %>%
             mutate (GWcal =  gwcal.b$Gwcal )



# vamos con ecd

ec_pheno.7.1_d <- ec_pheno.7.1 %>%
  filter ( essay == "ec2d")

setdiff (ec2d1$pot, ec_pheno.7.1_d$pot)

ec2d4 <- filter (ec2d1, 
                     pot!=3,
                     pot!=6,
                     pot!=10,
                     pot!=24,
                     pot!=28,
                     pot!=35,
                     pot!=38)

setdiff (ec2d4$pot, ec_pheno.7.1_d$pot)

 

# verifico que los Gw sean iguales 
ec2d4.to <- ec2d4 %>%
  filter (time ==0 )

ec2d4.t4 <- ec2d4 %>%
  filter (time ==4 )

ec2d4.t7 <- ec2d4 %>%
  filter (time ==7 )

ec2d4.t11 <- ec2d4 %>%
  filter (time ==11 )

plot (ec2d4.to$Gw, ec_pheno.7.1_d$gw0)
plot (ec2d4.t4$Gw, ec_pheno.7.1_d$gw4)
plot (ec2d4.t7$Gw, ec_pheno.7.1_d$gw7)
plot (ec2d4.t11$Gw, ec_pheno.7.1_d$gw11)


gw0cal.ec_pheno.7.1_d <- ec_pheno.7.1_d%>%
  select (essay : linea.disenio.1,gwcal.0 )%>%
  rename ( Gwcal = gwcal.0 ) %>%
  mutate (time =0)


gwc4cal.ec_pheno.7.1_d <- ec_pheno.7.1_d%>%
  select (essay : linea.disenio.1,gwcal.4 )%>%
  rename (Gwcal = gwcal.4)%>%
  mutate (time =4)

gwc7cal.ec_pheno.7.1_d <- ec_pheno.7.1_d%>%
  select (essay : linea.disenio.1,gwcal.7 )%>%
  rename (Gwcal = gwcal.7 )%>%
  mutate (time =7)

gwc8cal.ec_pheno.7.1_d <- ec_pheno.7.1_d%>%
  select (essay : linea.disenio.1)%>%
  mutate (Gwcal = NA )%>%
  mutate (time =8)


gwc11cal.ec_pheno.7.1_d <- ec_pheno.7.1_d%>%
  select (essay : linea.disenio.1,gwcal.11 )%>%
  rename (Gwcal = gwcal.11)%>%
  mutate (time =11)


gwc12cal.ec_pheno.7.1_d <- ec_pheno.7.1_d%>%
  select (essay : linea.disenio.1 )%>%
  mutate (Gwcal = NA)%>%
  mutate (time =12)



gwcal.d <- rbind  (gw0cal.ec_pheno.7.1_d,
                   gwc4cal.ec_pheno.7.1_d,
                   gwc7cal.ec_pheno.7.1_d,
                   gwc8cal.ec_pheno.7.1_d,
                   gwc11cal.ec_pheno.7.1_d,
                   gwc12cal.ec_pheno.7.1_d)

gwcal.d <- gwcal.d %>%
  arrange ( pot, time) 

#gwcal.d$pot <- as.character(gwcal.d$pot)

summary(gwcal.d )
str(gwcal.d)

#ec2d4$pot <- as.character(ec2d4$pot)

ec2d4 <-  ec2d4 %>%
          arrange (pot, time)

dim(ec2d4)
dim(gwcal.d)

(ec2d4$pot) == (gwcal.d$pot)

ec2d5 <- ec2d4 %>%
         mutate (GWcal =  gwcal.d$Gwcal )


## empiezo a modelar GW vs WSPM
############### ESTA ES LA FUNCION ##############

etiq <-function (X, Y, labs, m = c(0, 0), cex = 0.5, offset = 0.8, ...) 
{
  posXposY <- ((X >= m[1]) & ((Y >= m[2])))
  posXnegY <- ((X >= m[1]) & ((Y < m[2])))
  negXposY <- ((X < m[1]) & ((Y >= m[2])))
  negXnegY <- ((X < m[1]) & ((Y < m[2])))
  if (sum(posXposY) > 0) 
    text(X[posXposY], Y[posXposY], labs[posXposY],
         adj = c(0.5 - offset, 0.5 - offset), cex = cex, ...)
  if (sum(posXnegY) > 0) 
    text(X[posXnegY], Y[posXnegY], labs[posXnegY], 
         adj = c(0.5 -  offset, 0.5 + offset), cex = cex, ...)
  if (sum(negXposY) > 0) 
    text(X[negXposY], Y[negXposY], labs[negXposY], 
         adj = c(0.5 +  offset, 0.5 - offset), cex = cex, ...)
  if (sum(negXnegY) > 0) 
    text(X[negXnegY], Y[negXnegY], labs[negXnegY], 
         adj = c(0.5 + offset, 0.5 + offset), cex = cex, ...)
} 


GW.model <- function (dt = NULL, W=NULL, GX=NULL, essay=NULL, ue=NULL) {
  Gw   <-  dt %>% select(GX) # esta es la conductancia estomatica en cada dia
  ue   <- dt %>% select(ue)       # esta es la unidad experimental
  var1 <- as.vector (dt %>% select(W))      # esta es la variable independiente para modelar la conductancia
  essay <- dt %>% select(essay)
  id.ue <- paste (colnames(ue), ue [1,1],sep="_")
  id.out <-  paste ("ET.t",essay [1,1],id.ue,sep="_")

  mod1 <- lm ( Gw [,1] ~ var1 [,1])
  plot (var1 [,1], Gw [,1], 
        xlim=c(min(var1), max(var1)), 
        #ylim = c(0, na.omit(max(Gw))),
        pch=16, col="black",ylab="Gw", xlab="W(g)", 
        main= paste (id.out))
  
  abline(mod1, lwd=3 , col="red")
  # calculate residuals and predicted values
  res <- signif(residuals(mod1), 5)
  pre <- predict(mod1) # plot distances between points and the regression line
  #segments (var1, Gw, var1, pre, col="red")
  #etiq (var1, Gw, res, cex=0.7)
  
  a <- coef (mod1) [1]                             # este es el intercepto
  b <- coef (mod1) [2]                             # esta es la pemdiente
  
  par_GW_model <- data.frame (unique(essay),unique(ue), 
                              a = round (coef (mod1) [1] , 3), 
                              b = round (coef (mod1) [2], 3),
                              rsq = round (summary(mod1)$r.squared,3),
                              rsq.adj = round (summary(mod1)$adj.r.square,3))
  
  print (par_GW_model)
}

################ modelo  GW_w.SPM ##########333
head(ec2b.rep5.1)
ec2b.rep5.1 <- filter (ec2b.rep5, pot=="1")
ec2b.rep5.2 <- filter (ec2b.rep5, pot=="2")
ec2b.rep5.3 <- filter (ec2b.rep5, pot=="3")
ec2b.rep5.4 <- filter (ec2b.rep5, pot=="4")
ec2b.rep5.5 <- filter (ec2b.rep5, pot=="5")
ec2b.rep5.6 <- filter (ec2b.rep5, pot=="6")
#ec2b.rep5.7 <- filter (ec2b.rep5, pot=="7")
ec2b.rep5.8 <- filter (ec2b.rep5, pot=="8")
ec2b.rep5.9 <- filter (ec2b.rep5, pot=="9")
ec2b.rep5.10 <- filter (ec2b.rep5, pot=="10")
ec2b.rep5.11 <- filter (ec2b.rep5, pot=="11")
#ec2b.rep5.12 <- filter (ec2b.rep5, pot=="12")
#ec2b.rep5.13 <- filter (ec2b.rep5, pot=="13")
ec2b.rep5.14 <- filter (ec2b.rep5, pot=="14")
ec2b.rep5.15 <- filter (ec2b.rep5, pot=="15")
#ec2b.rep5.16 <- filter (ec2b.rep5, pot=="16")
#ec2b.rep5.17 <- filter (ec2b.rep5, pot=="17")
ec2b.rep5.18 <- filter (ec2b.rep5, pot=="18")
ec2b.rep5.19 <- filter (ec2b.rep5, pot=="19")
ec2b.rep5.20 <- filter (ec2b.rep5, pot=="20")
ec2b.rep5.21 <- filter (ec2b.rep5, pot=="21")
#ec2b.rep5.22 <- filter (ec2b.rep5, pot=="22")
ec2b.rep5.23 <- filter (ec2b.rep5, pot=="23")
ec2b.rep5.24 <- filter (ec2b.rep5, pot=="24")
ec2b.rep5.25 <- filter (ec2b.rep5, pot=="25")
#ec2b.rep5.26 <- filter (ec2b.rep5, pot=="26")
ec2b.rep5.27 <- filter (ec2b.rep5, pot=="27")
ec2b.rep5.28 <- filter (ec2b.rep5, pot=="28")
ec2b.rep5.29 <- filter (ec2b.rep5, pot=="29")
ec2b.rep5.30 <- filter (ec2b.rep5, pot=="30")
ec2b.rep5.31 <- filter (ec2b.rep5, pot=="31")
ec2b.rep5.32 <- filter (ec2b.rep5, pot=="32")
ec2b.rep5.33 <- filter (ec2b.rep5, pot=="33")
#ec2b.rep5.34 <- filter (ec2b.rep5, pot=="34")
ec2b.rep5.35 <- filter (ec2b.rep5, pot=="35")
ec2b.rep5.36 <- filter (ec2b.rep5, pot=="36")
ec2b.rep5.37 <- filter (ec2b.rep5, pot=="37")
ec2b.rep5.38 <- filter (ec2b.rep5, pot=="38")
ec2b.rep5.39 <- filter (ec2b.rep5, pot=="39")
ec2b.rep5.40 <- filter (ec2b.rep5, pot=="40")
ec2b.rep5.41 <- filter (ec2b.rep5, pot=="41")
ec2b.rep5.42 <- filter (ec2b.rep5, pot=="42")
ec2b.rep5.43 <- filter (ec2b.rep5, pot=="43")
ec2b.rep5.44 <- filter (ec2b.rep5, pot=="44")
ec2b.rep5.45 <- filter (ec2b.rep5, pot=="45")
ec2b.rep5.46 <- filter (ec2b.rep5, pot=="46")
#ec2b.rep5.47 <- filter (ec2b.rep5, pot=="47")
#ec2b.rep5.48 <- filter (ec2b.rep5, pot=="48")
ec2b.rep5.49 <- filter (ec2b.rep5, pot=="49")
#ec2b.rep5.50 <- filter (ec2b.rep5, pot=="50")
ec2b.rep5.51 <- filter (ec2b.rep5, pot=="51")
ec2b.rep5.52 <- filter (ec2b.rep5, pot=="52")
ec2b.rep5.53 <- filter (ec2b.rep5, pot=="53")
ec2b.rep5.54 <- filter (ec2b.rep5, pot=="54")
ec2b.rep5.55 <- filter (ec2b.rep5, pot=="55")
ec2b.rep5.56 <- filter (ec2b.rep5, pot=="56")
ec2b.rep5.57 <- filter (ec2b.rep5, pot=="57")
ec2b.rep5.58 <- filter (ec2b.rep5, pot=="58")
#ec2b.rep5.59 <- filter (ec2b.rep5, pot=="59")
ec2b.rep5.60 <- filter (ec2b.rep5, pot=="60")
ec2b.rep5.61 <- filter (ec2b.rep5, pot=="61")
ec2b.rep5.62 <- filter (ec2b.rep5, pot=="62")
ec2b.rep5.63 <- filter (ec2b.rep5, pot=="63")
ec2b.rep5.64 <- filter (ec2b.rep5, pot=="64")
ec2b.rep5.65 <- filter (ec2b.rep5, pot=="65")
ec2b.rep5.66 <- filter (ec2b.rep5, pot=="66")
ec2b.rep5.67 <- filter (ec2b.rep5, pot=="67")
ec2b.rep5.68 <- filter (ec2b.rep5, pot=="68")
#ec2b.rep5.69 <- filter (ec2b.rep5, pot=="69")
ec2b.rep5.70 <- filter (ec2b.rep5, pot=="70")
#ec2b.rep5.71 <- filter (ec2b.rep5, pot=="71")
#ec2b.rep5.72 <- filter (ec2b.rep5, pot=="72")
ec2b.rep5.73 <- filter (ec2b.rep5, pot=="73")
ec2b.rep5.74 <- filter (ec2b.rep5, pot=="74")
#ec2b.rep5.75 <- filter (ec2b.rep5, pot=="75")
ec2b.rep5.76 <- filter (ec2b.rep5, pot=="76")
ec2b.rep5.77 <- filter (ec2b.rep5, pot=="77")
ec2b.rep5.78 <- filter (ec2b.rep5, pot=="78")
ec2b.rep5.79 <- filter (ec2b.rep5, pot=="79")
ec2b.rep5.80 <- filter (ec2b.rep5, pot=="80")
ec2b.rep5.81 <- filter (ec2b.rep5, pot=="81")
ec2b.rep5.82 <- filter (ec2b.rep5, pot=="82")
ec2b.rep5.83 <- filter (ec2b.rep5, pot=="83")
ec2b.rep5.84 <- filter (ec2b.rep5, pot=="84")
ec2b.rep5.85 <- filter (ec2b.rep5, pot=="85")
ec2b.rep5.86 <- filter (ec2b.rep5, pot=="86")
ec2b.rep5.87 <- filter (ec2b.rep5, pot=="87")
ec2b.rep5.88 <- filter (ec2b.rep5, pot=="88")
ec2b.rep5.89 <- filter (ec2b.rep5, pot=="89")
ec2b.rep5.90 <- filter (ec2b.rep5, pot=="90")
ec2b.rep5.91 <- filter (ec2b.rep5, pot=="91")
ec2b.rep5.92 <- filter (ec2b.rep5, pot=="92")
ec2b.rep5.93 <- filter (ec2b.rep5, pot=="93")
#ec2b.rep5.94 <- filter (ec2b.rep5, pot=="94")
ec2b.rep5.95 <- filter (ec2b.rep5, pot=="95")
ec2b.rep5.96 <- filter (ec2b.rep5, pot=="96")
ec2b.rep5.97 <- filter (ec2b.rep5, pot=="97")
ec2b.rep5.98 <- filter (ec2b.rep5, pot=="98")
ec2b.rep5.99 <- filter (ec2b.rep5, pot=="99")
ec2b.rep5.100 <- filter (ec2b.rep5, pot=="100")
#ec2b.rep5.101 <- filter (ec2b.rep5, pot=="101")
ec2b.rep5.102 <- filter (ec2b.rep5, pot=="102")
ec2b.rep5.103 <- filter (ec2b.rep5, pot=="103")
ec2b.rep5.104 <- filter (ec2b.rep5, pot=="104")


head(ec2b.rep5.1)
pb5.1 <- GW.model (dt = ec2b.rep5.1,  GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.2 <- GW.model (dt = ec2b.rep5.2, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.3 <- GW.model (dt = ec2b.rep5.3, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.4 <- GW.model (dt = ec2b.rep5.4, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.5 <- GW.model (dt = ec2b.rep5.5, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.6 <- GW.model (dt = ec2b.rep5.6, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.7 <- GW.model (dt = ec2b.rep5.7, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.8 <- GW.model (dt = ec2b.rep5.8, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.9 <- GW.model (dt = ec2b.rep5.9, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.10 <- GW.model (dt = ec2b.rep5.10, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.11 <- GW.model (dt = ec2b.rep5.11, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.12 <- GW.model (dt = ec2b.rep5.12, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.13 <- GW.model (dt = ec2b.rep5.13, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.14 <- GW.model (dt = ec2b.rep5.14, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.15 <- GW.model (dt = ec2b.rep5.15, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.16 <- GW.model (dt = ec2b.rep5.16, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.17 <- GW.model (dt = ec2b.rep5.17, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.18 <- GW.model (dt = ec2b.rep5.18, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.19 <- GW.model (dt = ec2b.rep5.19, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.20 <- GW.model (dt = ec2b.rep5.20, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.21 <- GW.model (dt = ec2b.rep5.21, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.22 <- GW.model (dt = ec2b.rep5.22, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.23 <- GW.model (dt = ec2b.rep5.23, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.24 <- GW.model (dt = ec2b.rep5.24, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.25 <- GW.model (dt = ec2b.rep5.25, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.26 <- GW.model (dt = ec2b.rep5.26, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.27 <- GW.model (dt = ec2b.rep5.27, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.28 <- GW.model (dt = ec2b.rep5.28, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.29 <- GW.model (dt = ec2b.rep5.29, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.30 <- GW.model (dt = ec2b.rep5.30, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.31 <- GW.model (dt = ec2b.rep5.31, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.32 <- GW.model (dt = ec2b.rep5.32, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.33 <- GW.model (dt = ec2b.rep5.33, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.34 <- GW.model (dt = ec2b.rep5.34, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.35 <- GW.model (dt = ec2b.rep5.35, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.36 <- GW.model (dt = ec2b.rep5.36, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.37 <- GW.model (dt = ec2b.rep5.37, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.38 <- GW.model (dt = ec2b.rep5.38, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.39 <- GW.model (dt = ec2b.rep5.39, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.40 <- GW.model (dt = ec2b.rep5.40, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.41 <- GW.model (dt = ec2b.rep5.41, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.42 <- GW.model (dt = ec2b.rep5.42, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.43 <- GW.model (dt = ec2b.rep5.43, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.44 <- GW.model (dt = ec2b.rep5.44, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.45 <- GW.model (dt = ec2b.rep5.45, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.46 <- GW.model (dt = ec2b.rep5.46, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.47 <- GW.model (dt = ec2b.rep5.47, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.48 <- GW.model (dt = ec2b.rep5.48, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.49 <- GW.model (dt = ec2b.rep5.49, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.50 <- GW.model (dt = ec2b.rep5.50, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.51 <- GW.model (dt = ec2b.rep5.51, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.52 <- GW.model (dt = ec2b.rep5.52, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.53 <- GW.model (dt = ec2b.rep5.53, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.54 <- GW.model (dt = ec2b.rep5.54, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.55 <- GW.model (dt = ec2b.rep5.55, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.56 <- GW.model (dt = ec2b.rep5.56, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.57 <- GW.model (dt = ec2b.rep5.57, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.58 <- GW.model (dt = ec2b.rep5.58, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.59 <- GW.model (dt = ec2b.rep5.59, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.60 <- GW.model (dt = ec2b.rep5.60, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.61 <- GW.model (dt = ec2b.rep5.61, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.62 <- GW.model (dt = ec2b.rep5.62, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.63 <- GW.model (dt = ec2b.rep5.63, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.64 <- GW.model (dt = ec2b.rep5.64, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.65 <- GW.model (dt = ec2b.rep5.65, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.66 <- GW.model (dt = ec2b.rep5.66, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.67 <- GW.model (dt = ec2b.rep5.67, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.68 <- GW.model (dt = ec2b.rep5.68, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.69 <- GW.model (dt = ec2b.rep5.69, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.70 <- GW.model (dt = ec2b.rep5.70, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.71 <- GW.model (dt = ec2b.rep5.71, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.72 <- GW.model (dt = ec2b.rep5.72, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.73 <- GW.model (dt = ec2b.rep5.73, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.74 <- GW.model (dt = ec2b.rep5.74, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.75 <- GW.model (dt = ec2b.rep5.75, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.76 <- GW.model (dt = ec2b.rep5.76, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.77 <- GW.model (dt = ec2b.rep5.77, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.78 <- GW.model (dt = ec2b.rep5.78, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.79 <- GW.model (dt = ec2b.rep5.79, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.80 <- GW.model (dt = ec2b.rep5.80, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.81 <- GW.model (dt = ec2b.rep5.81, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.82 <- GW.model (dt = ec2b.rep5.82, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.83 <- GW.model (dt = ec2b.rep5.83, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.84 <- GW.model (dt = ec2b.rep5.84, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.85 <- GW.model (dt = ec2b.rep5.85, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.86 <- GW.model (dt = ec2b.rep5.86, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.87 <- GW.model (dt = ec2b.rep5.87, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.88 <- GW.model (dt = ec2b.rep5.88, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.89 <- GW.model (dt = ec2b.rep5.89, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.90 <- GW.model (dt = ec2b.rep5.90, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.91 <- GW.model (dt = ec2b.rep5.91, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.92 <- GW.model (dt = ec2b.rep5.92, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.93 <- GW.model (dt = ec2b.rep5.93, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.94 <- GW.model (dt = ec2b.rep5.94, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.95 <- GW.model (dt = ec2b.rep5.95, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.96 <- GW.model (dt = ec2b.rep5.96, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.97 <- GW.model (dt = ec2b.rep5.97, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.98 <- GW.model (dt = ec2b.rep5.98, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.99 <- GW.model (dt = ec2b.rep5.99, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.100 <- GW.model (dt = ec2b.rep5.100, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
#pb5.101 <- GW.model (dt = ec2b.rep5.101, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.102 <- GW.model (dt = ec2b.rep5.102, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.103 <- GW.model (dt = ec2b.rep5.103, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")
pb5.104 <- GW.model (dt = ec2b.rep5.104, GX="GWcal", W="w.SPM", essay="ec2b.rep5", ue="pot")

ec2b.GW.model.par <- rbind (pb5.1,pb5.2
                            #,pb5.3,pb5.4
                            ,pb5.5,pb5.6
                            #,pb5.7
                            ,pb5.8,pb5.9,pb5.10
                            ,pb5.11,
                            #pb5.12,
                            #pb5.13,
                            pb5.14,pb5.15,
                            #pb5.16,
                            #pb5.17,
                            pb5.18,pb5.19,pb5.20
                            ,pb5.21
                            #pb5.22
                            ,pb5.23,pb5.24,pb5.25
                            #pb5.26
                            ,pb5.27,pb5.28,pb5.29,pb5.30
                            ,pb5.31,pb5.32,pb5.33
                            #pb5.34
                            ,pb5.35,pb5.36,pb5.37,pb5.38,pb5.39,pb5.40
                            ,pb5.41,pb5.42,pb5.43,pb5.44,pb5.45,pb5.46
                            #,pb5.47,pb5.48
                            ,pb5.49
                            #,pb5.50,
                            ,pb5.51,pb5.52,pb5.53,pb5.54,pb5.55,pb5.56,pb5.57,pb5.58
                            #,pb5.59
                            ,pb5.60,pb5.61,pb5.62,pb5.63,pb5.64,pb5.65,pb5.66,pb5.67
                            ,pb5.68
                            #,pb5.69
                            ,pb5.70
                            #pb5.71
                            #,pb5.72
                            ,pb5.73,pb5.74
                            #pb5.75
                            ,pb5.76,pb5.77,pb5.78,pb5.79
                            ,pb5.80,pb5.81,pb5.82,pb5.83,pb5.84,pb5.85,pb5.86,pb5.87,pb5.88,pb5.89
                            ,pb5.90,pb5.91,pb5.92,pb5.93
                            #,pb5.94
                            ,pb5.95,pb5.96,pb5.97,pb5.98,pb5.99
                            ,pb5.100
                            #pb5.101
                            ,pb5.102,pb5.103,pb5.104)

dim (ec2b.GW.model.par)
dim (ec2b.rep5)

ec2b.rep6 <- ec2b.rep5 %>%
             filter(time == 0) %>%
             select (essay:trat)

ec2b.rep6$pot == ec2b.GW.model.par$pot

head(ec2b.GW.model.par)

ec2b.GW.model.par1 <- ec2b.GW.model.par %>%
                      select (a:rsq.adj)


ec2b.rep.GW.1 <- data.frame (ec2b.rep6, ec2b.GW.model.par1)

########### modelo para ec2d #############
# EC2d
ec2d5.1 <- filter (ec2d5, pot=="1")
ec2d5.2 <- filter (ec2d5, pot=="2")
#ec2d5.3 <- filter (ec2d5, pot=="3")
ec2d5.4 <- filter (ec2d5, pot=="4")
ec2d5.5 <- filter (ec2d5, pot=="5")
#ec2d5.6 <- filter (ec2d5, pot=="6")
ec2d5.7 <- filter (ec2d5, pot=="7")
ec2d5.8 <- filter (ec2d5, pot=="8")
ec2d5.9 <- filter (ec2d5, pot=="9")
#ec2d5.10 <- filter (ec2d5, pot=="10")
ec2d5.11 <- filter (ec2d5, pot=="11")
#ec2d5.12 <- filter (ec2d5, pot=="12")
ec2d5.13 <- filter (ec2d5, pot=="13")
ec2d5.14 <- filter (ec2d5, pot=="14")
ec2d5.15 <- filter (ec2d5, pot=="15")
ec2d5.16 <- filter (ec2d5, pot=="16")
ec2d5.17 <- filter (ec2d5, pot=="17")
ec2d5.18 <- filter (ec2d5, pot=="18")
ec2d5.19 <- filter (ec2d5, pot=="19")
ec2d5.20 <- filter (ec2d5, pot=="20")
ec2d5.21 <- filter (ec2d5, pot=="21")
ec2d5.22 <- filter (ec2d5, pot=="22")
ec2d5.23 <- filter (ec2d5, pot=="23")
#ec2d5.24 <- filter (ec2d5, pot=="24")
ec2d5.25 <- filter (ec2d5, pot=="25")
ec2d5.26 <- filter (ec2d5, pot=="26")
ec2d5.27 <- filter (ec2d5, pot=="27")
#ec2d5.28 <- filter (ec2d5, pot=="28")
ec2d5.29 <- filter (ec2d5, pot=="29")
#ec2d5.30 <- filter (ec2d5, pot=="30")
#ec2d5.31 <- filter (ec2d5, pot=="31")
ec2d5.32 <- filter (ec2d5, pot=="32")
ec2d5.33 <- filter (ec2d5, pot=="33")
ec2d5.34 <- filter (ec2d5, pot=="34")
#ec2d5.35 <- filter (ec2d5, pot=="35")
ec2d5.36 <- filter (ec2d5, pot=="36")
ec2d5.37 <- filter (ec2d5, pot=="37")
#ec2d5.38 <- filter (ec2d5, pot=="38")
ec2d5.39 <- filter (ec2d5, pot=="39")
ec2d5.40 <- filter (ec2d5, pot=="40")

pd5.1 <- GW.model (dt = ec2d5.1, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.2 <- GW.model (dt = ec2d5.2, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.3 <- GW.model (dt = ec2d5.3, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.4 <- GW.model (dt = ec2d5.4, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.5 <- GW.model (dt = ec2d5.5, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.6 <- GW.model (dt = ec2d5.6, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.7 <- GW.model (dt = ec2d5.7, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.8 <- GW.model (dt = ec2d5.8, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.9 <- GW.model (dt = ec2d5.9, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.10 <- GW.model (dt = ec2d5.10, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.11 <- GW.model (dt = ec2d5.11, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.12 <- GW.model (dt = ec2d5.12, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.13 <- GW.model (dt = ec2d5.13, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.14 <- GW.model (dt = ec2d5.14, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.15 <- GW.model (dt = ec2d5.15, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.16 <- GW.model (dt = ec2d5.16, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.17 <- GW.model (dt = ec2d5.17, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.18 <- GW.model (dt = ec2d5.18, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.19 <- GW.model (dt = ec2d5.19, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.20 <- GW.model (dt = ec2d5.20, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.21 <- GW.model (dt = ec2d5.21, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.22 <- GW.model (dt = ec2d5.22, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.23 <- GW.model (dt = ec2d5.23, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.24 <- GW.model (dt = ec2d5.24, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.25 <- GW.model (dt = ec2d5.25, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.26 <- GW.model (dt = ec2d5.26, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.27 <- GW.model (dt = ec2d5.27, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.28 <- GW.model (dt = ec2d5.28, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.29 <- GW.model (dt = ec2d5.29, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.30 <- GW.model (dt = ec2d5.30, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.31 <- GW.model (dt = ec2d5.31, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.32 <- GW.model (dt = ec2d5.32, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.33 <- GW.model (dt = ec2d5.33, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.34 <- GW.model (dt = ec2d5.34, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.35 <- GW.model (dt = ec2d5.35, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.36 <- GW.model (dt = ec2d5.36, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.37 <- GW.model (dt = ec2d5.37, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
#pd5.38 <- GW.model (dt = ec2d5.38, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.39 <- GW.model (dt = ec2d5.39, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")
pd5.40 <- GW.model (dt = ec2d5.40, GX="GWcal", W="w.SPM", essay="ec2d5", ue="pot")

ec2d.GW.model.par <- rbind (pd5.1,pd5.2
                            #pd5.3,
                            ,pd5.4,pd5.5
                            #,pd5.6
                            ,pd5.7,pd5.8,pd5.9
                            #pd5.10
                            ,pd5.11
                            #,pd5.12
                            ,pd5.13,pd5.14,pd5.15,pd5.16,pd5.17,pd5.18,pd5.19,pd5.20
                            ,pd5.21,pd5.22,pd5.23
                            #,pd5.24
                            ,pd5.25,pd5.26,pd5.27
                            #,pd5.28
                            ,pd5.29
                            #,pd5.30,pd5.31
                            ,pd5.32,pd5.33,pd5.34
                            #,pd5.35
                            ,pd5.36,pd5.37
                            #,pd5.38
                            ,pd5.39,pd5.40)

dim (ec2d.GW.model.par)
dim (ec2d5)

ec2d6 <- ec2d5 %>%
  filter(time == 0) %>%
  select (essay:trat)

ec2d6$pot == ec2d.GW.model.par$pot

head(ec2d.GW.model.par)

ec2d.GW.model.par1 <- ec2d.GW.model.par %>%
                      select (a:rsq.adj)


ec2d.GW.1 <- data.frame (ec2d6, ec2d.GW.model.par1)


### ## esta es la matriz para modelar 
model_GWAS_EC_2 <- rbind (ec2b.rep.GW.1, ec2d.GW.1)

dim (model_GWAS_EC_2)

(model_GWAS_EC_2_check <-  table (model_GWAS_EC_2$checks))

(model_GWAS_EC_2_geno <-  table (model_GWAS_EC_2$genotype))

write.table (model_GWAS_EC_2, file = "./Data/procdata/model_GWAS_EC_2.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)
