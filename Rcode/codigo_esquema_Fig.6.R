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
library("factoextra")
library("FactoMineR")
library("ggpubr")
library("corrplot")

########  se carga los datos de matriz que generan los histogramas de tm ####

pred.model.geno.DMSO.123 <- read.table ("./Data/procdata/pred.model.geno.DMSO.123", 
                                    header = TRUE, sep = ",",dec = ".",
                                    na.strings = "NA")

### clusterizo DMSO para tm #########
tm.DMSO <- pred.model.geno.DMSO.123$tm
shapiro.test (tm.DMSO)

q25.tm.DMSO <-  quantile (tm.DMSO, .25, na.rm = TRUE)
q75.tm.DMSO <-  quantile (tm.DMSO, .75, na.rm = TRUE)
q50.tm.DMSO <-  quantile (tm.DMSO, .50, na.rm = TRUE)

tm.DMSO.nq.25 <- pred.model.geno.DMSO.123 %>%
            filter (tm < q25.tm.DMSO) %>%
            mutate (clust.tm = "G1") %>%
            arrange (tm)

tm.DMSO.nq.25.75 <- pred.model.geno.DMSO.123 %>%
               filter (tm > q25.tm.DMSO) %>%
               filter (tm < q75.tm.DMSO) %>%
               mutate (clust.tm = "G2") %>%
               arrange (tm)

tm.DMSO.nq.75 <- pred.model.geno.DMSO.123 %>%
            filter (tm > q75.tm.DMSO)%>%
            mutate (clust.tm = "G3") %>%
            arrange (tm)


PCA.DMSO.tm <- rbind (tm.DMSO.nq.25,
                   tm.DMSO.nq.25.75, 
                   tm.DMSO.nq.75)

### exporto los datos ######


write_delim (PCA.DMSO.tm ,file ="./Data/procdata/PCA.DMSO.tm.txt" ,
             delim = ",", na = "NA")



# Selecciono las variables que voy a usar
PCA.DMSO.tm.a <- PCA.DMSO.tm %>%
                  dplyr::select(-c (datos.essay.1., pot, pred.rsq, k.eq.07, B.eq.07, 
                         k.eq.xy, B.eq.xy, tm.AR, clust.tm, tm, ks))

### Hago el PCA  ############
G88.DMSO.pca <- PCA ( PCA.DMSO.tm.a , scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                      quanti.sup = NULL, quali.sup =1, row.w = NULL, 
                      col.w = NULL, graph = FALSE)

## contribucion de las DIM 
eig.val.G88.DMSO <- get_eigenvalue (G88.DMSO.pca)
eig.val.G88.DMSO

## contribucion de las variables a las DIM
# Contributions to the principal components
(contrib <- G88.DMSO.pca$var$contrib)

### Plot de la contribucion de la dimensiones 
fviz_eig (G88.DMSO.pca, addlabels = TRUE, ylim = c(0, 50))

var.pheno.G88.DMSO <- get_pca_var (G88.DMSO.pca)
var.pheno.G88.DMSO

# Coordenadas de las variables
head (var.pheno.G88.DMSO$coord)

# Cos2: quality on the factore map
head (var.pheno.G88.DMSO$cos2)



# importante para enteder el PCA
# 1. Positively correlated variables are grouped together.
# 2. Negatively correlated variables are positioned on opposite sides 
#    of the plot origin (opposed quadrants).
# 3. The distance between variables and the origin measures 
#    the quality of the variables on the factor map.
#    Variables that are away from the origin are well 
#    represented on the factor map.


# Variables that are correlated with PC1 (i.e., Dim.1) and PC2 (i.e., Dim.2) 
#  are the most important in
#  explaining the variability in the data set.

# Variables that do not correlated with any PC or correlated
# with the last dimensions are variables
# with low contribution and might be removed to 
# simplify the overall analysis.

# The larger the value of the contribution, 
# the more the variable contributes to the component.

# Plot Contributions of variables to PC1
fviz_contrib (G88.DMSO.pca.1, choice = "var", axes = 1, top = 19)

# plot Contributions of variables to PC2
fviz_contrib (G88.DMSO.pca.1, choice = "var", axes = 2, top = 19)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)

res.km.G88.DMSO <- kmeans (var.pheno.G88.DMSO$coord, centers = 3 , nstart = 25)

grp.G88.DMSO <- as.factor(res.km.G88.DMSO$cluster)

# Plot la agrupacion de variables Color variables by groups
fviz_pca_var (G88.DMSO.pca, col.var = grp.G88.DMSO,
              repel =TRUE, 
              arrowsize=0.75,
              palette = c("red", "navyblue", "darkgreen"),
              legend.title = "Cluster")


##### individuos
ind.G88.DMSO <- get_pca_ind (G88.DMSO.pca)
ind.G88.DMSO

#fviz_pca_ind (G88.DMSO.pca.1)


#fviz_pca_ind (G88.DMSO.pca.1,geom = c("point"),
 #             palette = "jco", repel = TRUE)


# Clustering, auto nb of clusters:
hc.pheno.G88.DMSO <- HCPC (G88.DMSO.pca, nb.clust=-1)

clust.pca.pheno.G88.DMSO <- hc.pheno.G88.DMSO$data.clust 


# Plot cluster individuos 
fviz_pca_ind (G88.DMSO.pca,
              geom.ind = "point", # show points only (nbut not "text")
              col.ind = clust.pca.pheno.G88.DMSO$clust, # color by groups
              palette = c("darkorange", "navyblue", "red"),
              addEllipses = TRUE, # Concentration ellipses
              legend.title = "Groups")

# ordeno para  fusionar y  exportar

clust.pca.G88.DMSO <- hc.pheno.G88.DMSO$data.clust 

# ordeno para  fusionar y  exportar
clust.pca.G88.DMSO <- clust.pca.G88.DMSO %>%
                      dplyr::arrange(genotype) %>%
                      dplyr::select (c(genotype, clust))


PCA.DMSO <- PCA.DMSO.tm.a %>%
            dplyr::arrange (genotype)

PCA.DMSO$genotype == clust.pca.1.G88.DMSO$genotype

prev.CA.DMSO <- PCA.DMSO %>%
                dplyr::inner_join(clust.pca.G88.DMSO, by= "genotype") %>%
                dplyr::rename (clust.PCA = clust)

write.table (prev.CA.DMSO, file = "./Data/procdata/prev.CA.DMSO.txt", 
             append = FALSE, quote = TRUE, 
             sep = ",",
             eol = "\n", 
             na = "NA", 
             dec = ".", row.names = FALSE,
             col.names = TRUE)

############## PCA elite breeding #########3
pred.model.geno.ES.123 <- read.table ("./Data/procdata/pred.model.geno.ES.123.txt", 
                                      header = TRUE, sep = ",",dec = ".",
                                      na.strings = "NA")

### clusterizo elite #########
tm.elite <- pred.model.geno.ES.123$tm
shapiro.test (tm.elite)

q25.tm.elite <-  quantile (tm.elite, .25, na.rm = TRUE)
q75.tm.elite <-  quantile (tm.elite, .75, na.rm = TRUE)
q50.tm.elite <-  quantile (tm.elite, .50, na.rm = TRUE)

tm.elite.nq.25 <- pred.model.geno.ES.123 %>%
                  dplyr::filter (tm < q25.tm.elite) %>%
                  dplyr:: mutate (clust.tm = "G1") %>%
                  dplyr::arrange (tm)

tm.elite.nq.25.75 <- pred.model.geno.ES.123 %>%
                     dplyr::filter (tm > q25.tm.elite) %>%
                     dplyr::filter (tm < q75.tm.elite) %>%
                     dplyr::mutate (clust.tm = "G2") %>%
                     dplyr::arrange (tm)

tm.elite.nq.75 <- pred.model.geno.ES.123 %>%
                  dplyr::filter (tm > q75.tm.elite)%>%
                  dplyr::mutate (clust.tm = "G3") %>%
                  dplyr::arrange (tm)


PCA.elite <- rbind (tm.elite.nq.25,
                    tm.elite.nq.25.75, 
                    tm.elite.nq.75)

# PCA con la matrices del 
PCA.elite.a <- PCA.elite %>%
               select(-c (essay, pot, pred.rsq, k.eq.07, B.eq.07, 
                k.eq.xy, B.eq.xy, tm.AR, clust.tm, tm, ks))

G59.elite.pca <- PCA ( PCA.elite.a, scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                        quanti.sup = NULL, quali.sup =1, row.w = NULL, 
                        col.w = NULL, graph = FALSE)

plot (G59.elite.pca, choix = c("ind"))

plot (G59.elite.pca, choix = c("var"), title ="PCA_59_elite")


## contribucion de las DIM 
eig.val.G59.elite <- get_eigenvalue (G59.elite.pca)
eig.val.G59.elite

## contribucion de las variables a las DIM
# Contributions to the principal components
(contrib <- G59.elite.pca$var$contrib)

### Plot de la contribucion de la dimensiones 
fviz_eig (G59.elite.pca, addlabels = TRUE, ylim = c(0, 65))

var.pheno.G59.elite <- get_pca_var (G59.elite.pca)
var.pheno.G59.elite

# Coordenadas de las variables
head (var.pheno.G59.elite$coord)

# Cos2: quality on the factore map
head (var.pheno.G59.elite$cos2)

# Plot Contributions of variables to PC1
fviz_contrib (G59.elite.pca, choice = "var", axes = 1, top = 19)

# plot Contributions of variables to PC2
fviz_contrib (G59.elite.pca, choice = "var", axes = 2, top = 19)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)

res.km.G59.elite <- kmeans (var.pheno.G59.elite$coord, centers = 3 , nstart = 25)

grp.G59.elite <- as.factor(res.km.G59.elite$cluster)

# Plot la agrupacion de variables Color variables by groups
fviz_pca_var (G59.elite.pca, col.var = grp.G59.elite,
              repel =TRUE, 
              arrowsize=0.75,
              palette = c("red", "navyblue", "darkgreen"),
              legend.title = "Cluster")


##### individuos
ind.G59.elite <- get_pca_ind (G59.elite.pca)
ind.G59.elite


# Clustering, auto nb of clusters:
hc.pheno.G59.elite <- HCPC (G59.elite.pca, nb.clust=-1)

clust.pca.pheno.G59.elite <- hc.pheno.G59.elite$data.clust 


# Plot cluster individuos 
fviz_pca_ind (G59.elite.pca,
              geom.ind = "point", # show points only (nbut not "text")
              col.ind = clust.pca.pheno.G59.elite$clust, # color by groups
              palette = c("darkgreen", "navyblue", "red", "black" ),
              addEllipses = TRUE, # Concentration ellipses
              legend.title = "Groups")

# ordeno para  fusionar y  exportar

clust.pca.G59.elite <- hc.pheno.G59.elite$data.clust 

# ordeno para  fusionar y  exportar
clust.pca.G59.elite <- clust.pca.G59.elite %>%
                       dplyr::arrange(genotype) %>%
                       dplyr::select (c(genotype, clust))


PCA.elite <- PCA.elite %>%
            dplyr::arrange (genotype)

PCA.elite$genotype == clust.pca.G59.elite$genotype

prev.CA.elite <- PCA.elite %>%
                 dplyr::inner_join(clust.pca.G59.elite, by= "genotype") %>%
                 dplyr::rename (clust.PCA = clust)

write.table (prev.CA.elite, file = "./Data/procdata/prev.CA.elite.txt", 
             append = FALSE, quote = TRUE, 
             sep = ",",
             eol = "\n", 
             na = "NA", 
             dec = ".", row.names = FALSE,
             col.names = TRUE)

########  se cargan los datos #########
## datos el modelo
CA.DMSO <- read.table ("./Data/rawdata/CA.DMSO.q25.txt", 
                      header = TRUE, sep = "\t",dec = ".",
                      na.strings = "NA")


row.names (CA.DMSO) <- CA.DMSO$var.factor

summary (CA.DMSO)

# cluster tm
CA.DMSO.tm <- CA.DMSO %>%
             select (-var.factor)

res.ca.DMSO.tm <- CA (CA.DMSO.tm, row.sup = NULL, col.sup = NULL)

plot (res.ca.DMSO.tm, title="CA.DMSO")


########## 
CA.elite <- read.table ("./Data/rawdata/CA.elite.q25.txt", 
                       header = TRUE, sep = "\t",dec = ".",
                       na.strings = "NA")
str (CA.elite)
row.names (CA.elite) <- CA.elite$var.factor

summary (CA.elite)

# cluster tm
CA.elite.tm <- CA.elite %>%
  select (-var.factor)

res.ca.elite.tm <- CA (CA.elite.tm, row.sup = NULL, col.sup = NULL)

plot (res.ca.elite.tm, title="CA.elite")







row.names(CA.55G.15.id13c) <- CA.55G.15.id13c$var.factor


CA.55G.15.id13c <- CA.55G.15.id13c %>%
  select (-var.factor)


res.ca.id13c.15 <- CA (CA.55G.15.id13c, row.sup = NULL, col.sup = NULL)

summary(res.ca.id13c.15)

svg (filename="./Figures/Fig.5/CA.id13C.2015.svg", 
     width=7, 
     height=5, 
     pointsize=12)



dev.off ()



# cluster GM 
CA.55G.15.GM <- CA.55G.15 %>%
  filter (var.numeric < 10) %>%
  filter ( 4 < var.numeric) %>%
  select (-c(anio,var.numeric))

row.names(CA.55G.15.GM) <- CA.55G.15.GM$var.factor


CA.55G.15.GM <- CA.55G.15.GM %>%
  select (-var.factor)


res.ca.GM.15 <- CA (CA.55G.15.GM, row.sup = NULL, col.sup = NULL)

summary (res.ca.GM.15)


### cluster IS
CA.55G.15.IS <- CA.55G.15 %>%
  filter ( 9 < var.numeric) %>%
  select (-c(anio,var.numeric))

row.names(CA.55G.15.IS) <- CA.55G.15.IS$var.factor


CA.55G.15.IS <- CA.55G.15.IS %>%
  select (- var.factor)


res.ca.IS.15 <- CA (CA.55G.15.IS, row.sup = NULL, col.sup = NULL)

summary (res.ca.IS.15)

svg (filename="./Figures/Fig.5/CA.IS.2015.svg", 
     width=7, 
     height=5, 
     pointsize=12)

plot (res.ca.IS.15, title="CA.IS.2015")

dev.off ()






write.table (x, file = "", append = FALSE, quote = TRUE, sep = " ",
             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
             col.names = TRUE, qmethod = c("escape", "double"),
             fileEncoding = "")








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
x2.tm <- max (which (dens.tm$x <= q10.tm))

with (dens.tm, polygon(x=c(x[c(x1.tm,x1.tm:x2.tm,x2.tm)]),
                         y= c(0, y[x1.tm:x2.tm], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))

x3.tm  <- min (which (dens.tm$x >= q90.tm))
x4.tm  <- max (which (dens.tm$x > 0))
with (dens.tm, polygon(x=c(x[c(x3.tm,x3.tm:x4.tm,x4.tm)]),
                         y= c(0, y[x3.tm:x4.tm], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))

abline (v=q50.tm, col="navyblue", lwd=1.5, lty =2)


box()
dev.off()

########### 
svg (filename="./Figures/Fig.4/hist.ab.svg", 
     width=7, 
     height=5, 
     pointsize=12)

ab <- pred.model.geno.DMSO.123$ab
shapiro.test (ab)

#q25.ab <-  quantile (ab, .25, na.rm = TRUE)
#q75.ab <-  quantile (ab, .75, na.rm = TRUE)

q10.ab <-  quantile (ab, .10, na.rm = TRUE)
q90.ab <-  quantile (ab, .90, na.rm = TRUE)
q50.ab <-  quantile (ab, .50, na.rm = TRUE)

ab.nq.10 <- pred.model.geno.DMSO.123 %>%
           filter (ab < q10.ab)%>%
           arrange (ab)

ab.nq.90 <- pred.model.geno.DMSO.123 %>%
           filter (ab > q90.ab)%>%
           arrange (ab)

h.ab  <- hist (ab , breaks="Sturges")
h.ab$density <- h.ab$counts/sum (h.ab$counts)*100


plot (h.ab, freq=FALSE, xlab="em.ab",main = "",
      ylab="Percent of total",
      col=NULL
      ,xlim = c(450,600)
      ,xaxt='l')

par(new=TRUE)

hist (ab, main="ab.geno", 
      freq=F, xlab="", ylab = "",
      xlim = c(450,600),
      axes=FALSE,col="gray85")

lines (density (ab, na.rm = T), lwd=2,
       col="navyblue", lty=1)

dens.ab <- density (ab, na.rm = T)
x1.ab <- min (which (dens.ab$x >= 0))
x2.ab <- max (which (dens.ab$x <= q10.ab))

with (dens.ab, polygon(x=c(x[c(x1.ab,x1.ab:x2.ab,x2.ab)]),
                         y= c(0, y[x1.ab:x2.ab], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))

x3.ab  <- min (which (dens.ab$x >= q90.ab))
x4.ab  <- max (which (dens.ab$x > 0))
with (dens.ab, polygon(x=c(x[c(x3.ab,x3.ab:x4.ab,x4.ab)]),
                         y= c(0, y[x3.ab:x4.ab], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))

abline (v=q50.ab, col="navyblue", lwd=1.5, lty =2)

box()
dev.off()

########### 
svg (filename="./Figures/Fig.4/hist.Gwtm.svg", 
     width=7, 
     height=5, 
     pointsize=12)

Gwtm <- pred.model.geno.DMSO.123$Gwtm
shapiro.test (Gwtm)

#q25.dGwtm <-  quantile (dGwtm, .25, na.rm = TRUE)
#q75.dGwtm <-  quantile (dGwtm, .75, na.rm = TRUE)

q10.Gwtm <-  quantile (Gwtm, .10, na.rm = TRUE)
q90.Gwtm <-  quantile (Gwtm, .90, na.rm = TRUE)
q50.Gwtm <-  quantile (Gwtm, .50, na.rm = TRUE)

Gwtm.nq.10 <- pred.model.geno.DMSO.123 %>%
  filter (Gwtm < q10.Gwtm)%>%
  arrange (Gwtm)

Gwtm.nq.90 <-pred.model.geno.DMSO.123 %>%
  filter (Gwtm > q90.Gwtm)%>%
  arrange (Gwtm)

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
x2.Gwtm <- max (which (dens.Gwtm$x <= q10.Gwtm))

with (dens.Gwtm, polygon(x=c(x[c(x1.Gwtm,x1.Gwtm:x2.Gwtm,x2.Gwtm)]),
                          y= c(0, y[x1.Gwtm:x2.Gwtm], 0),border ='navyblue',
                          col=scales::alpha( 'navyblue',.5)))

x3.Gwtm  <- min (which (dens.Gwtm$x >= q90.Gwtm))
x4.Gwtm  <- max (which (dens.Gwtm$x > 0))
with (dens.Gwtm, polygon(x=c(x[c(x3.Gwtm,x3.Gwtm:x4.Gwtm,x4.Gwtm)]),
                          y= c(0, y[x3.Gwtm:x4.Gwtm], 0),border ='navyblue',
                          col=scales::alpha( 'navyblue',.5)))

abline (v=q50.Gwtm, col="navyblue", lwd=1.5, lty =2)

box()
dev.off()


#################### derivada en el tm 
svg (filename="./Figures/Fig.4/hist.dGwtm.svg", 
     width=7, 
     height=5, 
     pointsize=12)

dGwtm <- pred.model.geno.DMSO.123$dGwtm
shapiro.test (dGwtm)

#q25.dGwtm <-  quantile (dGwtm, .25, na.rm = TRUE)
#q75.dGwtm <-  quantile (dGwtm, .75, na.rm = TRUE)

q10.dGwtm <-  quantile (dGwtm, .10, na.rm = TRUE)
q90.dGwtm <-  quantile (dGwtm, .90, na.rm = TRUE)
q50.dGwtm <-  quantile (dGwtm, .50, na.rm = TRUE)

dGwtm.nq.10 <- pred.model.geno.DMSO.123 %>%
  filter (dGwtm < q10.dGwtm)%>%
  arrange (dGwtm)

dGwtm.nq.90 <-pred.model.geno.DMSO.123 %>%
  filter (dGwtm > q90.dGwtm)%>%
  arrange (dGwtm)

h.dGwtm  <- hist (dGwtm , breaks="Sturges")
h.dGwtm$density <- h.dGwtm$counts/sum (h.dGwtm$counts)*100


plot (h.dGwtm, freq=FALSE, xlab="em.dGwtm",main = "",
      ylab="Percent of total",
      col=NULL
      #,xlim = c(400,1000)
      ,xaxt='l')

par(new=TRUE)

hist (dGwtm, main="dGwtm.geno", 
      freq=F, xlab="", ylab = "",
      axes=FALSE,col="gray85")

lines (density (dGwtm, na.rm = T), lwd=2,
       col="navyblue", lty=1)

dens.dGwtm <- density (dGwtm, na.rm = T)
x1.dGwtm <- min (which (dens.dGwtm$x >= 0))
x2.dGwtm <- max (which (dens.dGwtm$x <= q10.dGwtm))

with (dens.dGwtm, polygon(x=c(x[c(x1.dGwtm,x1.dGwtm:x2.dGwtm,x2.dGwtm)]),
                       y= c(0, y[x1.dGwtm:x2.dGwtm], 0),border ='navyblue',
                       col=scales::alpha( 'navyblue',.5)))

x3.dGwtm  <- min (which (dens.dGwtm$x >= q90.dGwtm))
x4.dGwtm  <- max (which (dens.dGwtm$x > 0))
with (dens.dGwtm, polygon(x=c(x[c(x3.dGwtm,x3.dGwtm:x4.dGwtm,x4.dGwtm)]),
                       y= c(0, y[x3.dGwtm:x4.dGwtm], 0),border ='navyblue',
                       col=scales::alpha( 'navyblue',.5)))

abline (v=q50.dGwtm, col="navyblue", lwd=1.5, lty =2)

box()
dev.off()

