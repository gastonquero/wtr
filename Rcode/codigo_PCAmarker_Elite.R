##############################################################################
# Datos SNPs SOJA                                                           #
# Matriz de datos Soja lineas avanzadas                                      #
# Los datos fueron recibidos el 01/03/2017                                  #
#                                                                            #
#                                                                            #                                            
# Datos originales dados por Mariana Menoni                                  #
#                                                                            #
# Gaston Quero                                                               #
##############################################################################


# fijar Directorio
getwd()
#setwd ("E:/Paper_Soja_G3")
#setwd ("C:/Users/Usuario/OneDrive/Documentos/Paper_Crop_Science")

setwd ("R:/wtr")

# Cargar Paquetes
library("lmem.gwaser", lib.loc="~/R/win-library/3.3")
library("lmem.qtler", lib.loc="~/R/win-library/3.3")
library("qtl")
qtlversion()
library("FactoMineR")
library("dplyr")
library("car")

# Paquetes 
library (lme4)
library (emmeans)
library ("nlmrt")
library ("easynls")
library (tidyverse)
library ("lattice")
library ("latticeExtra")
library (multcompView)
library("dplyr")
library("ggridges")
library("ggjoy")
library("viridis")
library("lmerTest")
library("lubridate")
library (ggcorrplot)
library (sjPlot)
library (sjlabelled)
library (sjmisc)
library (ggplot2)
library (coefplot2)
library (ggpubr)
library(qtl)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")


# Cargar datos como data frame

group.feno.Elite <- read_delim ( file = "./Data/procdata/group.feno.Elite.txt",
                                delim =";", quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                                col_names = TRUE, col_types = NULL, na = "NA")


group.feno.Elite <- mutate (group.feno.Elite,genotype = fct_recode  (genotype,"DM6.8" = "testigo1"))
group.feno.Elite <- mutate (group.feno.Elite,genotype = fct_recode  (genotype,"NA5509RG" = "testigo2"))
group.feno.Elite <- mutate (group.feno.Elite,genotype = fct_recode  (genotype,"LEO1706-07" = "testigo3"))
group.feno.Elite <- mutate (group.feno.Elite,genotype = fct_recode  (genotype,"DM5.9i" = "testigo4"))
group.feno.Elite <- mutate (group.feno.Elite,genotype = fct_recode  (genotype,"NA5009RG" = "testigo5"))

group.feno.Elite <- group.feno.Elite%>%
                    dplyr::mutate (genotypes = paste ("G",genotype,sep="_"))%>%
                    dplyr::select (-genotype)


ids.pheno.elite <- group.feno.Elite %>%
                   dplyr::select(genotypes)


## Geno data
G.soja.data <- read.table ("./Data/rawdata/soja_geno.2.txt",
                      header = T, sep = "\t",
                      dec = ".", na.strings = "-")
dim (G.soja.data)
str (G.soja.data)

G.soja.data <- G.soja.data %>%
               dplyr::arrange (genotypes)

setdiff (ids.pheno.elite$genotypes, G.soja.data$genotypes )

setdiff ( G.soja.data$genotypes , ids.pheno.elite$genotypes)

#ids.pheno.elite <- mutate (ids.pheno.elite,genotypes = fct_recode  (genotypes,"G_LEO_1706_07"= "G_LEO1706-07"))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes,  "G_NA_5909_RG"= "G_NIDERA_A_5909_RG" ))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes,  "G_LEO1706-07"= "G_LEO_1706_07" ))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes,   "G_NA5509RG"="G_NIDERA_A_5509_RG"))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes, "G_DM5.9i"= "G_DON_MARIO_5.9i" ))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes,  "G_DM6.8"= "G_DON_MARIO_6.8_i"))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes, "G_NA5009RG" ="G_NIDERA_A_5009_RG"))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes,  "G_NA_6126_RG"= "G_NIDERA_A_6126_RG" ))
G.soja.data  <- mutate (G.soja.data,genotypes = fct_recode  (genotypes,  "G_SRM_6900" ="G_RM_6900"))

setdiff (ids.pheno.elite$genotypes, G.soja.data$genotypes )


G.soja.elite <- G.soja.data %>%
                    dplyr::filter (genotypes %in% ids.pheno.elite$genotypes )

# Exporto los datos de la poblacion elite para generar el cross

write_delim (G.soja.elite, file="./Data/procdata/G.soja.elite.txt",
            delim = ";")
####
elite.cross <- read.cross (format="csv",
                           dir="./Data/procdata",
                           file="elite.cross.csv", 
                           na.strings= "NA",sep=";",dec=",",
                           genotypes=c("AA","BB"), 
                           alleles=c("A","B"),
                           estimate.map=FALSE, crosstype="dh",
                           convertXdata=FALSE, error.prob=0.0001)

summary(elite.cross)

nind(elite.cross)
nchr(elite.cross)
totmar(elite.cross)
nmar(elite.cross)
nphe(elite.cross)


###### vamos a sacar los misssing data ##########
plotMissing (elite.cross) 

# Eliminar marcadores con mas del 50 % de datos faltantes
n.missing <- nmissing(elite.cross, what="mar")[(nmissing(elite.cross, what="mar"))/sum(summary(elite.cross)$n.ind)> 0.50]  
names.marker <- c (names(n.missing))
length(names.marker)

elite.cross.1 <- drop.markers (elite.cross, names.marker)

plotMissing(elite.cross.1, main="elite.cross.1")
summary(elite.cross.1)
nmar(elite.cross.1)
totmar(elite.cross.1)

# eliminar los individuos con mas de 50 no genotipado
elite.cross.2 <- subset(elite.cross.1, 
                        ind=(ntyped(elite.cross.1) > ((totmar(elite.cross.1) * 50)/100)))

indiv <- subset(elite.cross.1, 
                ind=(ntyped(elite.cross.1)  < ((totmar(elite.cross.1) * 50)/100)))

indiv$pheno$id
plotMissing(elite.cross.2, main="elite.cross.2")

nmar(elite.cross.2)
totmar(elite.cross.2)

summary (elite.cross.2 )

write.cross (elite.cross.2, format="csv",
             filestem="./Data/procdata/elite.cross.2")

# cargo la matriz para el PCA
Elite_dat <-  read_delim ( file ="./Data/procdata/elite.cross.2.csv"  , delim =",",
                           quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                           col_names = TRUE, col_types = NULL, na = "-")

Elite_dat  <- Elite_dat [-c(1,2),]

Elite_dat  [Elite_dat == "AA"] <- "0"
Elite_dat  [Elite_dat  == "BB"] <- "1"
Elite_dat  [Elite_dat  == "AB"] <- "0.5"


write_delim (Elite_dat ,file ="./Data/procdata/Elite_dat.1.txt" ,
             delim = ",", na = "NA")

## reingreso los datos 
Elite_dat.1  <- read_delim ( file ="./Data/procdata/Elite_dat.1.txt", delim =",",
                            quote = "\"", escape_backslash = FALSE, escape_double = TRUE,
                            col_names = TRUE, col_types = NULL, na = "NA")

Elite_dat.1.pca <- PCA (Elite_dat.1, scale.unit = TRUE, ncp = 5, 
                          ind.sup = NULL, quanti.sup = NULL, 
                          quali.sup = 1, row.w = NULL, 
                          col.w = NULL, graph = TRUE, axes = c(1,2))


## contribucion de las DIM 
eig.val.G53.Elite.G <- get_eigenvalue ( Elite_dat.1.pca)
eig.val.G53.Elite.G

## contribucion de las variables a las DIM
# Contributions to the principal components
(contrib <-  Elite_dat.1.pca$var$contrib)


### Plot de la contribucion de la dimensiones 
fviz_eig (Elite_dat.1.pca, addlabels = TRUE, ylim = c(0, 7.5))

row.names (Elite_dat.1.pca)


var.pheno.G53.Elite.G <- get_pca_var (Elite_dat.1.pca)

var.pheno.G53.Elite.G.1 <- var.pheno.G53.Elite.G$contrib [,1] > 0.1

class ( var.pheno.G53.Elite.G)


# Coordenadas de las variables
head (var.pheno.G53.Elite.G$coord)
# Cos2: quality on the factore map
head (var.pheno.G53.Elite.G$cos2)



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
fviz_contrib (Elite_dat.1.pca, choice = "var", axes = 1, top = 19)

# plot Contributions of variables to PC2
fviz_contrib (Elite_dat.1.pca, choice = "var", axes = 2, top = 19)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)

res.mrk.G53.Elite <- kmeans (var.pheno.G53.Elite.G$coord, centers = 3 , nstart = 25)

grp.G53.Elite.G <- as.factor(res.mrk.G53.Elite$cluster)

# Plot la agrupacion de variables Color variables by groups
fviz_pca_var (Elite_dat.1.pca, col.var = grp.G53.Elite.G,
              repel =TRUE, 
              arrowsize=0.75,
              palette = c("red", "navyblue", "darkgreen"),
              legend.title = "Cluster")

##### individuos
ind.G53.Elite.G <- get_pca_ind (Elite_dat.1.pca)
ind.G53.Elite.G

#fviz_pca_ind (G53.Elite.pca.1)


#fviz_pca_ind (G53.Elite.pca.1,geom = c("point"),
#             palette = "jco", repel = TRUE)


# Clustering, auto nb of clusters:
hc.pheno.G53.Elite.G <- HCPC (Elite_dat.1.pca, nb.clust=-1)

clust.pca.pheno.G53.Elite.G <- hc.pheno.G53.Elite.G$data.clust 


# Plot cluster individuos 

fviz_pca_ind (Elite_dat.1.pca,
              geom.ind = "point", # show points only (nbut not "text")
              pointsize = 2,
              col.ind = clust.pca.pheno.G53.Elite.G$clust, # color by groups
              #addEllipses = TRUE,
              palette = c("green", "red", "navyblue", "darkorange", "gray48", "black"),
              ellipse.level=0.95,# Concentration ellipses
              legend.title = "Groups")+ 
  theme_minimal()+
  scale_shape_manual(values=c(19, 14, 15, 16,17,8))


# ordeno para  fusionar y  exportar

clust.pca.G53.Elite.G <-  hc.pheno.G53.Elite.G$data.clust 

head (clust.pca.G53.Elite.G )
# ordeno para  fusionar y  exportar
clust.pca.G53.Elite.G <- clust.pca.G53.Elite.G %>%
                         dplyr::arrange(genotypes) %>%
                         dplyr::select (c(genotypes, clust))

prev.CA.Elite.G <- group.feno.Elite  %>%
                   dplyr::inner_join( clust.pca.G53.Elite.G, by= "genotypes") %>%
                   dplyr::rename (clust.PCA = clust) %>%
                   dplyr::select (c(genotypes, group.tm, group.Gwtm, clust.PCA ,everything()))


write.table (prev.CA.Elite.G, file = "./Data/procdata/prev.CA.Elite.G.txt", 
             append = FALSE, quote = TRUE, 
             sep = ",",
             eol = "\n", 
             na = "NA", 
             dec = ".", row.names = FALSE,
             col.names = TRUE)

######## veo de armar el data frame para el CA ###

prev.CA.Elite.G.1 <- prev.CA.Elite.G %>%
                     dplyr::select (genotypes, group.tm, group.Gwtm, clust.PCA) 

#### para la agrupacion con tm ########3
prev.CA.Elite.G.1$clust.PCA <- as.factor (prev.CA.Elite.G.1$clust.PCA )
prev.CA.Elite.G.1$group.tm <- as.factor  (prev.CA.Elite.G.1$group.tm )

tm.PCApheno <-  prev.CA.Elite.G.1 %>%
                dplyr::group_by (group.tm) %>%
                count(clust.PCA) %>%
                dplyr::ungroup()

CA.Elite.tm.PCApheno  <- tm.PCApheno %>%
                          tidyr::pivot_wider(names_from = clust.PCA, values_from = n)

colnames (CA.Elite.tm.PCApheno) <- c("var.factor", "clust.1", "clust.2", "clust.3","clust.4", "clust.5", "clust.6" )

CA.Elite.tm.PCApheno <- as.data.frame(CA.Elite.tm.PCApheno )

row.names (CA.Elite.tm.PCApheno) <- CA.Elite.tm.PCApheno$var.factor

summary (CA.Elite.tm.PCApheno)

# cluster tm
CA.Elite.tm.PCApheno <- CA.Elite.tm.PCApheno %>%
                       dplyr::select (-var.factor)


res.CA.Elite.tm.PCApheno <- CA (CA.Elite.tm.PCApheno, row.sup = NULL, col.sup = NULL)

plot (res.CA.Elite.tm.PCApheno, title="CA.x")

#### para la agrupacion con GWtm ########3
prev.CA.Elite.G.1$clust.PCA <- as.factor (prev.CA.Elite.G.1$clust.PCA )
prev.CA.Elite.G.1$group.Gwtm <- as.factor  (prev.CA.Elite.G.1$group.Gwtm )

Gwtm.PCApheno <- prev.CA.Elite.G.1 %>%
                 dplyr::group_by (group.Gwtm) %>%
                 count(clust.PCA) %>%
                 dplyr::ungroup()

CA.Elite.Gwtm.PCApheno  <- Gwtm.PCApheno %>%
  tidyr::pivot_wider(names_from = clust.PCA, values_from = n)

colnames (CA.Elite.Gwtm.PCApheno) <- c("var.factor", "clust.1", "clust.2", "clust.3","clust.4", "clust.5", "clust.6" )

CA.Elite.Gwtm.PCApheno <- as.data.frame(CA.Elite.Gwtm.PCApheno )

row.names (CA.Elite.Gwtm.PCApheno) <- CA.Elite.Gwtm.PCApheno$var.factor

summary (CA.Elite.Gwtm.PCApheno)

# cluster Gwtm
CA.Elite.Gwtm.PCApheno <- CA.Elite.Gwtm.PCApheno %>%
                         dplyr::select (-var.factor)

CA.Elite.Gwtm.PCApheno [3,4] <-0
res.CA.Elite.Gwtm.PCApheno <- CA (CA.Elite.Gwtm.PCApheno, row.sup = NULL, col.sup = NULL)

plot (res.CA.Elite.Gwtm.PCApheno, title="CA.x")
