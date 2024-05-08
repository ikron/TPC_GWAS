##Run GAPIT for gr35 with 4 PCs

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

genomeans <- read.csv(file = "/projappl/project_2000350/gwas/data/temp_genomeans_new.csv", header = T, sep = ",")

gr35 <- dplyr::filter(genomeans, Temp == 35)

myY <- gr35[gr35$Genot %in% geno.names,]
myY <- myY[,-c(2,3)] #Drop the temp column, not needed
colnames(myY) <- c("Taxa", "gr35")

##Run GAPIT
myGAPIT35 <- GAPIT(Y = myY, PCA.total = 4, G = myG, model = "Blink")
