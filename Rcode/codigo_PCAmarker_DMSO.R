##############################################################################
# Datos SNPs SOJA                                                            #
# Matriz de datos Soja cruza DM68xSO76557                                    #
# Los datos fueron recibidos el 06/03/2018  ????                             #                                                                        #
#                                                                            #                                            
# Datos originales dados Victoria Bonnecarrere                               #
#                                                                            #
# Gaston Quero                                                               #
##############################################################################

# fijar Directorio

getwd()
setwd ("R:/wtr")

# Cargar Paquetes
library("qtl")
library("ggplot2")
library("FactoMineR")
library("dplyr")
library("car")
library (ggjoy)
library (hrbrthemes)
library(tidyverse)
library(forcats)
library("dplyr")
library("viridis")
library (lme4)
library (emmeans)
library ("car")
library ("nlmrt")
library ("easynls")
library ("lattice")
library ("latticeExtra")
library (multcompView)
library("ggridges")
library("ggjoy")
library("lmerTest")
library("lubridate")
library(ggcorrplot)
library("factoextra")
library("FactoMineR")
library("ggpubr")
library("corrplot")


# Cargar datos como data frame
# Estos son los datos originales 

DM.SO.cross <- read.cross (format="csv",
               dir="./Data/rawdata",
               file="previo_DM.SO_cross.csv", 
               na.strings= "NA",sep=";",dec=",",
               genotypes=c("AA","BB","AB"), alleles=c("A","B"),
               estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)

summary(DM.SO.cross)

nind(DM.SO.cross)
nchr(DM.SO.cross)
totmar(DM.SO.cross)
nmar(DM.SO.cross)
nphe(DM.SO.cross)


###### vamos a sacar los misssing data ##########
plotMissing (DM.SO.cross) 
# Eliminar marcadores con mas del 50 % de datos faltantes
n.missing <- nmissing(DM.SO.cross, what="mar")[(nmissing(DM.SO.cross, what="mar"))/sum(summary(DM.SO.cross)$n.ind)> 0.50]  
names.marker <- c (names(n.missing))
length(names.marker)

DM.SO.cross.1 <- drop.markers (DM.SO.cross, names.marker)

plotMissing(DM.SO.cross.1, main="DM.SO.cross.1")
summary(DM.SO.cross.1)
nmar(DM.SO.cross.1)
totmar(DM.SO.cross.1)

# eliminar los individuos con mas de 50 no genotipado
DM.SO.cross.2 <- subset (DM.SO.cross.1, 
                 ind=(ntyped(DM.SO.cross.1) > ((totmar(DM.SO.cross.1) * 50)/100)))

indiv <- subset(DM.SO.cross.1, 
                ind=(ntyped(DM.SO.cross.1)  < ((totmar(DM.SO.cross.1) * 50)/100)))

indiv$pheno$id
plotMissing(DM.SO.cross.2, main="DM.SO.cross.2")

nmar(DM.SO.cross.2)
totmar(DM.SO.cross.2)

###### voy a comparar los genotipos 
cg <- comparegeno(DM.SO.cross.2)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

########## distorcion de segregacion 
gt <- geno.table(DM.SO.cross.2)
gt[gt$P.value < 0.05/totmar(DM.SO.cross.2),]


### calculo de las frecuencias 
total <- nind (DM.SO.cross.2)
gt <- geno.table(DM.SO.cross.2)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4,5)]
Geno.freq <- as.vector(obs.alelo/num.obs)
head(Geno.freq)

##### No polimorficos para AA
no.polimorficos.AA <- Geno.freq[,"AA"] >= 0.8
no.poli.AA <- rownames(Geno.freq [no.polimorficos.AA,])
length(no.poli.AA)
names.marker.AA <- c (rownames(Geno.freq [no.poli.AA,]))
5405

##### No polimorficos para BB
no.polimorficos.BB <- Geno.freq[,"BB"] >= 0.8
no.poli.BB <- rownames(Geno.freq [no.polimorficos.BB,])
length(no.poli.BB)

##### No polimorficos para AB
no.polimorficos.AB <- Geno.freq[,"AB"] >= 0.8
no.poli.AB <- rownames(Geno.freq [no.polimorficos.AB,])
length(no.poli.AB)
names.marker.AB <- c (rownames(Geno.freq [no.poli.AB,]))


DM.SO.cross.2a <- drop.markers (DM.SO.cross.2, names.marker.AA)
totmar(DM.SO.cross.2a)

DM.SO.cross.2b <- drop.markers (DM.SO.cross.2a, names.marker.AB)
totmar(DM.SO.cross.2b)
summary (DM.SO.cross.2b)



total <- nind (DM.SO.cross.2b)
gt <- geno.table(DM.SO.cross.2b)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

SNPs.MAF <-  Geno.freq[,2] < 0.05 
rownames(Geno.freq [SNPs.MAF,])
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)

DM.SO.cross.2b <- drop.markers (DM.SO.cross.2b, names.marker)

summary (DM.SO.cross )

write.cross (DM.SO.cross.2b, format="csv",
             filestem="./Data/procdata/DM.SO.cross.2b")


write.cross (DM.SO.cross.2b, format="tidy",
             filestem="./Data/procdata/DM.SO.cross.2b")


write.cross (DM.SO.cross, format="tidy",
             filestem="./Data/procdata/DM.SO.cross")

### Ahora cargo los datos para hacer el PCA ###

#DM.SO_dat <-  read_delim ( file ="./Data/rawdata/DM.SO_data.txt"  , delim ="\t",
 #                          quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
  #                         col_names = TRUE, col_types = NULL, na = "-")

# cargo la matriz transpuesta 
DM.SO_dat <-  read_delim ( file ="./Data/procdata/DM.SO.cross.2b.csv", 
                           delim =",",quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                           col_names = TRUE, col_types = NULL, na = "-")

DM.SO_dat.x <- DM.SO_dat [-(1:2), -2 ]


DM.SO_dat.1 <- DM.SO_dat.x %>%
               dplyr::mutate (id.1 = str_replace (id, "[-]", ".")) %>%
               dplyr::select (id.1, everything())
  
  
#ids_DM.SO.cross.2b  <-  read_delim ( file ="./Data/procdata/DM.SO.cross.2b_phe.csv", delim =",",
 #                          quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
  #                         col_names = TRUE, col_types = NULL, na = "NA")

#list.ids_DM.SO.cross.2b <- colnames (ids_DM.SO.cross.2b) [-1]


#head (DM.SO_dat.1)

### filtro  segun la matrix filtrada de los marcadores

#DM.SO_dat.2 <- DM.SO_dat.1 %>%
 #              dplyr::filter (id.1 %in% list.ids_DM.SO.cross.2b )


### cargo la correspondencia de id para los datos fenotipicos 

ids.pheno <-  read_delim ( file ="./Data/rawdata/ids_DMSO.txt"  , delim ="\t",
                           quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                           col_names = TRUE, col_types = NULL, na = "-")

ids.pheno.1 <- ids.pheno %>%
               dplyr::mutate (id.1 = str_replace (id, "[-]", ".")) %>%
               dplyr::select (id.1, everything()) %>%
               dplyr::mutate (id.1 = str_c ("G_",id.1)) %>%
               dplyr::select (c(id.1,genotypes, lineas ))

### cargo la matriz de donde hice los agrupamiento de la figur 5
ids.pheno.1$lineas



group.feno.DMSO <- read_delim ( file = "./Data/procdata/group.feno.DMSO.txt",
                                delim =";", quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                                col_names = TRUE, col_types = NULL, na = "NA")

group.feno.DMSO$genotype

##### uno las dos matrices  ##########

ids.pheno.2 <- ids.pheno.1 %>%
               dplyr::filter (lineas %in% group.feno.DMSO$genotype)
  

DM.SO_dat.3 <-  DM.SO_dat.1 %>%
                dplyr::inner_join (ids.pheno.2, by="id.1" ) %>%
                dplyr::select (id.1,  id , genotypes, lineas, everything())


DM.SO_dat.3  [DM.SO_dat.3  == "AA"] <- "0"
DM.SO_dat.3  [DM.SO_dat.3  == "BB"] <- "1"
DM.SO_dat.3  [DM.SO_dat.3  == "AB"] <- "0.5"


 write_delim (DM.SO_dat.3 ,file ="./Data/procdata/DM.SO_dat.4.txt" ,
                            delim = ",", na = "NA")
 
## reingreso los datos 
 DM.SO_dat.4 <- read_delim ( file ="./Data/procdata/DM.SO_dat.4.txt", delim =",",
                            quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                            col_names = TRUE, col_types = NULL, na = "NA")
 
 DM.SO_dat.4.1 <- DM.SO_dat.4 %>%
                  dplyr::select (-genotypes) %>% 
                  dplyr::rename (genotype = lineas ) %>%
                  dplyr::select (-c(id.1, id))
 
 head (DM.SO_dat.4.1)
 
 ### me quedo con las variables que voya usar para el CA
 prev.CA.DMSO <- group.feno.DMSO %>%
                 dplyr::select (genotype, group.tm, group.Gwtm)
 
 
 DM.SO_dat.4.2 <-  DM.SO_dat.4.1 %>%
                   dplyr::inner_join( prev.CA.DMSO, by="genotype") %>%
                   dplyr::select (genotype, group.tm, group.Gwtm, everything())
 
DM.SO_dat.4.2.pca <- PCA (DM.SO_dat.4.2, scale.unit = TRUE, ncp = 5, 
                            ind.sup = NULL, quanti.sup = NULL, 
                            quali.sup = c(1:3), row.w = NULL, 
                            col.w = NULL, graph = TRUE, axes = c(1,2))


 ## contribucion de las DIM 
 eig.val.G88.DMSO.G <- get_eigenvalue ( DM.SO_dat.4.2.pca)
 eig.val.G88.DMSO.G
 
 ## contribucion de las variables a las DIM
 # Contributions to the principal components
 (contrib <-  DM.SO_dat.4.2.pca$var$contrib)

 
 ### Plot de la contribucion de la dimensiones 
 fviz_eig (DM.SO_dat.4.2.pca, addlabels = TRUE, ylim = c(0, 7.5))
 
 row.names (DM.SO_dat.4.2.pca)
 
 
 var.pheno.G88.DMSO.G <- get_pca_var (DM.SO_dat.4.2.pca)
 
 var.pheno.G88.DMSO.G.1 <- var.pheno.G88.DMSO.G$contrib [,1] > 0.1
 
 class ( var.pheno.G88.DMSO.G)


 # Coordenadas de las variables
 head (var.pheno.G88.DMSO.G$coord)
 # Cos2: quality on the factore map
 head (var.pheno.G88.DMSO.G$cos2)
 
 
 
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
 fviz_contrib (DM.SO_dat.4.2.pca, choice = "var", axes = 1, top = 19)
 
 # plot Contributions of variables to PC2
 fviz_contrib (DM.SO_dat.4.2.pca, choice = "var", axes = 2, top = 19)
 
 # Create a grouping variable using kmeans
 # Create 3 groups of variables (centers = 3)
 set.seed(123)
 
 res.mrk.G88.DMSO <- kmeans (var.pheno.G88.DMSO.G$coord, centers = 3 , nstart = 25)
 
 grp.G88.DMSO.G <- as.factor(res.mrk.G88.DMSO$cluster)
 
 # Plot la agrupacion de variables Color variables by groups
 fviz_pca_var (DM.SO_dat.4.2.pca, col.var = grp.G88.DMSO.G,
               repel =TRUE, 
               arrowsize=0.75,
               palette = c("red", "navyblue", "darkgreen"),
               legend.title = "Cluster")
 
 ##### individuos
 ind.G88.DMSO.G <- get_pca_ind (DM.SO_dat.4.2.pca)
 ind.G88.DMSO.G
 
 #fviz_pca_ind (G88.DMSO.pca.1)
 
 
 #fviz_pca_ind (G88.DMSO.pca.1,geom = c("point"),
 #             palette = "jco", repel = TRUE)
 
 
 # Clustering, auto nb of clusters:
 hc.pheno.G88.DMSO.G <- HCPC (DM.SO_dat.4.2.pca, nb.clust=-1)
 
 clust.pca.pheno.G88.DMSO.G <- hc.pheno.G88.DMSO.G$data.clust 
 
 
 # Plot cluster individuos 
 fviz_pca_ind (DM.SO_dat.4.2.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = clust.pca.pheno.G88.DMSO.G$clust, # color by groups
               #addEllipses = TRUE,
               palette = c("black", "darkorange", "navyblue", "red"),
               ellipse.level=0.95,# Concentration ellipses
               legend.title = "Groups")
 

 fviz_pca_ind (DM.SO_dat.4.2.pca,
               geom.ind = "point", # show points only (nbut not "text")
               pointsize = 2,
               col.ind = clust.pca.pheno.G88.DMSO.G$clust, # color by groups
               #addEllipses = TRUE,
               palette = c("black", "darkorange", "navyblue"),
               ellipse.level=0.95,# Concentration ellipses
               legend.title = "Groups")+ 
               theme_minimal()+
               scale_shape_manual(values=c(19,17,18))
 
 
 # ordeno para  fusionar y  exportar
 
 clust.pca.G88.DMSO.G <-  hc.pheno.G88.DMSO.G$data.clust 
 
 # ordeno para  fusionar y  exportar
 clust.pca.G88.DMSO.G <- clust.pca.G88.DMSO.G %>%
                         dplyr::arrange(genotype) %>%
                         dplyr::select (c(genotype, clust))
 
 prev.CA.DMSO.G <- DM.SO_dat.4.2  %>%
                    dplyr::inner_join( clust.pca.G88.DMSO.G, by= "genotype") %>%
                    dplyr::rename (clust.PCA = clust) %>%
                    dplyr::select (c(genotype, group.tm, group.Gwtm, clust.PCA ,everything() ))
 
 
 write.table (prev.CA.DMSO.G, file = "./Data/procdata/prev.CA.DMSO.G.txt", 
              append = FALSE, quote = TRUE, 
              sep = ",",
              eol = "\n", 
              na = "NA", 
              dec = ".", row.names = FALSE,
              col.names = TRUE)
 
 ######## veo de armar el data frame para el CA ###
 
 prev.CA.DMSO.G.1 <- prev.CA.DMSO.G %>%
                     dplyr::select (genotype, group.tm, group.Gwtm, clust.PCA) 
                      
#### para la agrupacion con tm ########3
 prev.CA.DMSO.G.1$clust.PCA <- as.factor (prev.CA.DMSO.G.1$clust.PCA )
 prev.CA.DMSO.G.1$group.tm <- as.factor  (prev.CA.DMSO.G.1$group.tm )
 
tm.PCApheno <-  prev.CA.DMSO.G.1 %>%
                dplyr::group_by (group.tm) %>%
                         count(clust.PCA) %>%
                dplyr::ungroup()
 
CA.DMSO.tm.PCApheno  <- tm.PCApheno %>%
                        tidyr::pivot_wider(names_from = clust.PCA, values_from = n)

colnames (CA.DMSO.tm.PCApheno) <- c("var.factor", "clust.1", "clust.2", "clust.3")

CA.DMSO.tm.PCApheno <- as.data.frame(CA.DMSO.tm.PCApheno )

row.names (CA.DMSO.tm.PCApheno) <- CA.DMSO.tm.PCApheno$var.factor

summary (CA.DMSO.tm.PCApheno)

# cluster tm
CA.DMSO.tm.PCApheno <- CA.DMSO.tm.PCApheno %>%
                       dplyr::select (-var.factor)

CA.DMSO.tm.PCApheno [3,1] <-0
res.CA.DMSO.tm.PCApheno <- CA (CA.DMSO.tm.PCApheno, row.sup = NULL, col.sup = NULL)

plot (res.CA.DMSO.tm.PCApheno, title="CA.x")

#### para la agrupacion con GWtm ########3
prev.CA.DMSO.G.1$clust.PCA <- as.factor (prev.CA.DMSO.G.1$clust.PCA )
prev.CA.DMSO.G.1$group.Gwtm <- as.factor  (prev.CA.DMSO.G.1$group.Gwtm )

Gwtm.PCApheno <- prev.CA.DMSO.G.1 %>%
                 dplyr::group_by (group.Gwtm) %>%
                 count(clust.PCA) %>%
  dplyr::ungroup()

CA.DMSO.Gwtm.PCApheno  <- Gwtm.PCApheno %>%
                          tidyr::pivot_wider(names_from = clust.PCA, values_from = n)

colnames (CA.DMSO.Gwtm.PCApheno) <- c("var.factor", "clust.1", "clust.2", "clust.3")

CA.DMSO.Gwtm.PCApheno <- as.data.frame(CA.DMSO.Gwtm.PCApheno )

row.names (CA.DMSO.Gwtm.PCApheno) <- CA.DMSO.Gwtm.PCApheno$var.factor

summary (CA.DMSO.Gwtm.PCApheno)

# cluster Gwtm
CA.DMSO.Gwtm.PCApheno <- CA.DMSO.Gwtm.PCApheno %>%
                         dplyr::select (-var.factor)

CA.DMSO.Gwtm.PCApheno [1,3] <-0
res.CA.DMSO.Gwtm.PCApheno <- CA (CA.DMSO.Gwtm.PCApheno, row.sup = NULL, col.sup = NULL)

plot (res.CA.DMSO.Gwtm.PCApheno, title="CA.x")



