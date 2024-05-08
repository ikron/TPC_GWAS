##Run GAPIT for gr20 with 4 PCs

# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))
libpath <- .libPaths()[1]

#Load sources and libraries
source("/projappl/project_2000350/scripts/Neuro_functions.R") #This works
source("/projappl/project_2000350/gwas/association_scripts.R")
library(dplyr)

##Run association analysis

##Load genotype data

myG <- read.table(file = "/projappl/project_2000350/gwas/data/Neuro_hapmap_natpop_gwas_temp.txt", header = FALSE, stringsAsFactors = FALSE)
geno.names <- myG[1,-c(1:11)]

##Load phenotypic data

splineres <- read.csv(file = "/projappl/project_2000350/gwas/data/splineresults.csv", header = T, sep = ",")

myY <- splineres[splineres$Genotype %in% geno.names,]
myY <- myY[,c(3,2)] #Drop the umax column and reorder column, not needed
colnames(myY) <- c("Taxa", "Topt")

##Run GAPIT
myGAPIT20 <- GAPIT(Y = myY, PCA.total = 4, G = myG, model = "Blink")
