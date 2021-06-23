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
DM.SO.cross.2 <- subset(DM.SO.cross.1, 
                        ind=(ntyped(DM.SO.cross.1) > ((totmar(DM.SO.cross.1) * 70)/100)))

indiv <- subset(DM.SO.cross.1, 
                ind=(ntyped(DM.SO.cross.1)  < ((totmar(DM.SO.cross.1) * 70)/100)))

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


write.cross (DM.SO.cross.2b, format="tidy",
             filestem="./Data/procdata/DM.SO.cross.2b")


write.cross (DM.SO.cross, format="tidy",
             filestem="./Data/procdata/DM.SO.cross")

### Ahora cargo los datos para hacer el PCA ###

DM.SO_dat <-  read_delim ( file ="./Data/rawdata/DM.SO_data.txt"  , delim ="\t",
                           quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                           col_names = TRUE, col_types = NULL, na = "-")

DM.SO_dat.1 <- DM.SO_dat %>%
               dplyr::mutate (id.1 = str_replace (id, "[-]", ".")) %>%
               dplyr::select (id.1, everything())
  
  
ids_DM.SO.cross.2b  <-  read_delim ( file ="./Data/procdata/DM.SO.cross.2b_phe.csv", delim =",",
                           quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                           col_names = TRUE, col_types = NULL, na = "NA")

list.ids_DM.SO.cross.2b <- colnames (ids_DM.SO.cross.2b) [-1]


head (DM.SO_dat.1)
### filtro  segun la matrix filtrada de los marcadores

DM.SO_dat.2 <- DM.SO_dat.1 %>%
               dplyr::filter (id.1 %in% list.ids_DM.SO.cross.2b )


### cargo la correspondencia de id para los datos fenotipicos 

ids.pheno <-  read_delim ( file ="./Data/rawdata/ids_DMSO.txt"  , delim ="\t",
                           quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                           col_names = TRUE, col_types = NULL, na = "-")


ids.pheno.1 <- ids.pheno %>%
               dplyr::mutate (id.1 = str_replace (id, "[-]", ".")) %>%
               dplyr::select (id.1, everything()) %>%
               dplyr::mutate (id.1 = str_c ("G_",id.1)) %>%
               dplyr::select (c(id.1,genotypes, lineas ))
 

##### uno las dos matrices  ##########
head (DM.SO_dat.2)

DM.SO_dat.3 <-  DM.SO_dat.2 %>%
                dplyr::inner_join(ids.pheno.1, by="id.1" ) %>%
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
 
 ### ingreso los datos de los datos fenotipicos ###
 
 prev.CA.DMSO.pheno <- read_delim ( file = "./Data/procdata/prev.CA.DMSO.txt", delim =",",
                             quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                             col_names = TRUE, col_types = NULL, na = "NA")
 
 head (prev.CA.DMSO.pheno)
 
 prev.CA.DMSO.pheno.1 <- prev.CA.DMSO.pheno %>%
                         dplyr::select (genotype, clust.tm, clust.PCA)
 
 
 DM.SO_dat.4.2 <-  DM.SO_dat.4.1 %>%
                   dplyr::inner_join( prev.CA.DMSO.pheno.1, by="genotype") %>%
                   dplyr::select (genotype, clust.tm, clust.PCA, everything())
 
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
               addEllipses = TRUE,
               palette = c("black", "darkorange", "navyblue", "red"),
               ellipse.level=0.95,# Concentration ellipses
               legend.title = "Groups")
 
 
 fviz_pca_ind (DM.SO_dat.4.2.pca,
               geom.ind = "point", # show points only (nbut not "text")
               pointsize = 2,
               col.ind = clust.pca.pheno.G88.DMSO.G$clust, # color by groups
               addEllipses = TRUE,
               palette = c("black", "darkorange", "navyblue", "red"),
               ellipse.level=0.95,# Concentration ellipses
               legend.title = "Groups")+ 
               theme_minimal()+
               scale_shape_manual(values=c(19,19,19,19))
 
 
 # ordeno para  fusionar y  exportar
 
 clust.pca.G88.DMSO.G <-  hc.pheno.G88.DMSO.G$data.clust 
 
 # ordeno para  fusionar y  exportar
 clust.pca.G88.DMSO.G <- clust.pca.G88.DMSO.G %>%
                         dplyr::arrange(genotype) %>%
                         dplyr::select (c(genotype, clust))
 
 prev.CA.DMSO.G <- DM.SO_dat.4.2  %>%
                    dplyr::inner_join( clust.pca.G88.DMSO.G, by= "genotype") %>%
                    dplyr::rename (clust.PCA.G = clust) %>%
                    dplyr::select (c(genotype,clust.tm , clust.PCA, clust.PCA.G ,everything() ))
 
 
 write.table (prev.CA.DMSO.G, file = "./Data/procdata/prev.CA.DMSO.G.txt", 
              append = FALSE, quote = TRUE, 
              sep = ",",
              eol = "\n", 
              na = "NA", 
              dec = ".", row.names = FALSE,
              col.names = TRUE)
 
 ######## veo de armar el data frame para el CA ###
 
 prev.CA.DMSO.G.1 <- prev.CA.DMSO.G %>%
                     dplyr::select (genotype, clust.tm, clust.PCA, clust.PCA.G ) 
                      
 
 prev.CA.DMSO.G.1$clust.PCA <- as.factor (prev.CA.DMSO.G.1$clust.PCA )
 prev.CA.DMSO.G.1$clust.tm <- as.factor  (prev.CA.DMSO.G.1$clust.tm )
 
tm.PCApheno <-  prev.CA.DMSO.G.1 %>%
                dplyr::group_by (clust.tm) %>%
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

res.CA.DMSO.tm.PCApheno <- CA (CA.DMSO.tm.PCApheno, row.sup = NULL, col.sup = NULL)

plot (res.CA.DMSO.tm.PCApheno, title="CA.x")


### cluster tm geno #######

tm.PCAgeno <- prev.CA.DMSO.G.1 %>%
              dplyr::group_by (clust.tm) %>%
              count(clust.PCA.G) %>%
              dplyr::ungroup()

CA.DMSO.tm.PCAgeno  <- tm.PCAgeno %>%
                        tidyr::pivot_wider(names_from = clust.PCA.G, values_from = n)

colnames (CA.DMSO.tm.PCAgeno) <- c("var.factor", "cmrk.1", "cmrk.2", "cmrk.3", "cmrk.4")

CA.DMSO.tm.PCAgeno <- as.data.frame(CA.DMSO.tm.PCAgeno )

row.names (CA.DMSO.tm.PCAgeno) <- CA.DMSO.tm.PCAgeno$var.factor

summary (CA.DMSO.tm.PCAgeno)

# cluster tm
CA.DMSO.tm.PCAgeno <- CA.DMSO.tm.PCAgeno %>%
                      dplyr::select (-var.factor)

res.CA.DMSO.tm.PCAgeno <- CA (CA.DMSO.tm.PCAgeno, row.sup = NULL, col.sup = NULL)

plot (res.CA.DMSO.tm.PCAgeno, title="CA.x")

### todas las variables 
### cluster pheno geno #######

PCApheno.PCAgeno <- prev.CA.DMSO.G.1 %>%
              dplyr::group_by (clust.PCA) %>%
              count(clust.PCA.G) %>%
              dplyr::ungroup()

CA.DMSO.PCApheno.PCAgeno <- PCApheno.PCAgeno %>%
                              tidyr::pivot_wider(names_from = clust.PCA.G, values_from = n)

CA.DMSO.PCApheno.PCAgeno <- CA.DMSO.PCApheno.PCAgeno %>%
                            dplyr::mutate (clust.PCA = str_c ("phe", clust.PCA))


colnames (CA.DMSO.PCApheno.PCAgeno) <- c("var.factor", "cmrk.1", "cmrk.2", "cmrk.3", "cmrk.4")

CA.DMSO.PCApheno.PCAgeno <- as.data.frame (CA.DMSO.PCApheno.PCAgeno)

row.names (CA.DMSO.PCApheno.PCAgeno) <- CA.DMSO.PCApheno.PCAgeno$var.factor

summary (CA.DMSO.PCApheno.PCAgeno)


CA.DMSO.PCApheno.PCAgeno <- CA.DMSO.PCApheno.PCAgeno %>%
                            dplyr::select (-var.factor)

CA.DMSO.PCApheno.PCAgeno[is.na(CA.DMSO.PCApheno.PCAgeno)] <- 0


res.CA.DMSO.PCApheno.PCAgeno <- CA (CA.DMSO.PCApheno.PCAgeno, row.sup = NULL, col.sup = NULL)

plot (res.CA.DMSO.PCApheno.PCAgeno, title="CA.x")


