## Analysis and association mapping of constant temperatures

##################################################################
### Set up folders to load scripts, load libraries             ###
### Allows for easier loading of scripts if file paths change  ###
##################################################################

#Folder for script files
scriptfolder <- "~/Documents/tutkijatohtori/scripts/" #Change this as necessary
source(paste(scriptfolder, "Neuro_functions.R", sep = ""))

#Folder for association scripts
assofolder <- "~/Documents/tutkijatohtori/association/" #Change this as necessary
source(paste(assofolder, "association_scripts.R", sep = ""))

#Data folder
datafolder <- "./data/" #Change this is necessary
#Results folder
resultsfolder <- "./results/"

#Load libraries
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(brms)
library(coda)
library(latex2exp)
library(data.table)

### * Load and process phenotypic data

aineisto <- read.csv(paste(datafolder, "phenotypes.csv", sep = ""))

##Calculating growth rates
gr <- calc.lin.gr(aineisto, 8:15, 5, tindex)
aineisto <- cbind(aineisto, gr)

### Make the genotypic means file (Needs to be run only once)

genomeans <- summarise(group_by(aineisto, Genot, Temp, Family), meangr = mean(growthrate, na.rm = T))

write.table(genomeans, file = paste(datafolder, "genomeans.csv", sep = ""), sep = ",", quote = F, row.names = F)

### Load genotypic means

genomeans <- read.csv(paste(datafolder, "genomeans.csv", sep = ""), header = T, sep = ",")

### * Summarise the data

##Plot of genotypic means
my.xlab <- expression(paste("Temperature (", degree,"C)"))
p1 <- ggplot(genomeans, aes(y = meangr, x = Temp, group = Genot)) +
    geom_point(alpha = 0.1) +
    geom_line(alpha = 0.1) +
    ylab("Growth rate (mm / h)") +
    scale_y_continuous(limits = c(0, 5.5), breaks = seq(from = 0, to = 5.5, by = 0.5)) +
    xlab(my.xlab)

### * Spline fits

## using splines for each genotype and then estimating optimum temp and broadness
library(splines)

#This function fits a spline for a single genotype and the result
plot.spline <- function(aineisto, genotype)
{
    datamat <- aineisto[aineisto$Genot == genotype,] #Take data for a genotype
    fit <- spline(datamat$Temp, datamat$growthrate, n = 50, method = "natural") #Fit the spline
    fit <- as.data.frame(fit) #change to dataframe

    my.xlab <- expression(paste("Temperature (", degree,"C)"))
    ggplot(fit, aes(x = x, y = y)) +
        geom_line() +
        geom_point(data = datamat, aes(x = Temp, y = growthrate)) +
        ylab("Growth rate (mm / h)") +
        xlab(my.xlab)     
    
}

plot.spline(aineisto, genotype = "XG23")
plot.spline(aineisto, genotype = "X2489")
plot.spline(aineisto, genotype = "P4486")
plot.spline(aineisto, genotype = "XG9")
plot.spline(aineisto, genotype = "X10911")
plot.spline(aineisto, genotype = "X10917")

##Constructing a datalist for splines, Predicted fit for each genotype
##pred.data <- data.frame(
##Some genotypes need to be dropped: G2
splinedata <- filter(aineisto, Genot != "XG2")
splinedata$Genot <- factor(splinedata$Genot)
foo <- list(0)
datalist <- rep(foo, nlevels(splinedata$Genot))
allgenot <- unique(splinedata$Genot)
#splineres <- data.frame(Genotype = allgenot, umax = rep(0, length(allgenot)), T
for(i in 1:nlevels(splinedata$Genot)) {
    current <- filter(splinedata, Genot == allgenot[i])
    print(as.character(allgenot[i]))
    fit <- spline(current$Temp, current$growthrate, n = 100, method = "natural")
    datalist[[i]] <- data.frame(Temp = fit$x, growthrate = fit$y)
}


splineres <- matrix(rep(0,length(datalist)*2), ncol = 2)
splineres <- data.frame(splineres)
colnames(splineres) <- c("umax", "Topt")
##Search for umax and Topt from spline fits
for(i in 1:length(datalist)) {
    splineres[i,1] <- max(datalist[[i]][,2]) #Get umax
    umax.ind <- which(datalist[[i]][,2] == max(datalist[[i]][,2])) #Get index of umax
    splineres[i,2] <- datalist[[i]][umax.ind,1] }

splineres$Genotype <- allgenot

#Save the spline results for GWAS
write.table(splineres, file = "./data/splineresults.csv", sep = ",", row.names = F, quote = F)

### Loading spline results
splineres <- read.csv(paste0(datafolder, "splineresults.csv"), header = T, sep =",")

#Generate covariate tables that include umax
#Load the PCAs for covariates
PCAs <- read.table("./data/GAPIT.Genotype.PCA.csv", header = T, sep = ",")

sortorder <- match(splineres[,3], PCAs[,1])
PCAs <- PCAs[sortorder,]
PCAs$umax <- splineres$umax
PCAs <- filter(PCAs, is.na(taxa) == F) #There was one row (XG9) with missing values

#Write a table that includes PCA covariates and maximum growth rate (to control for elevation)
write.table(PCAs, file = "./data/PCAs_umax.csv", sep = ",", row.names = F, quote = F)

#cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)   


##Figure of the spline fits
my.xlab <- expression(paste("Temperature (", degree,"C)"))
pdf(file = "splines.pdf")
ggplot() +
    lapply(datalist, function(dat) {
        geom_line(data = dat, aes(x = Temp, y = growthrate), alpha = 0.1) } ) +
    ylab("Growth rate (mm / h)") +
    xlab(my.xlab) +
    scale_y_continuous(limits = c(0,6), breaks = seq(0, 5.5, 0.5) ) +
    scale_x_continuous(limits = c(20, 40), expand = c(0,0), breaks = seq(20, 40, 5)) +
    theme(plot.margin = margin(t = 10, r = 15))
dev.off()

##Making a data.frame with genotype and T_opt and umax
splineres <- data.frame(genotype = allgenot, umax = splineres[,1], Topt = splineres[,2])

##Save the spline data
#write.table(splineres, file = "temp_spline_est.csv", sep = ",", row.names = FALSE)
 

### * Association mapping

##Association mapping to run a cluster

### ** Analysis of GWAS results

### *** Load the data

#Folder for association data for a given trait (Change when necessary)
assodata <- "~/Genomics/Neurospora/association/"

#For growth rate at 20 C
#asso.gr20.MLM <- read.csv(paste(assodata, "gr20/GAPIT.Association.GWAS_Results.MLM.gr20.csv", sep = ""))
#asso.gr20.MLMM <- read.csv(paste(assodata, "gr20/GAPIT.Association.GWAS_Results.MLMM.gr20.csv", sep = ""))
asso.gr20.BLINK <- read.csv(paste(assodata, "gr20/GAPIT.Association.GWAS_Results.BLINK.gr20.csv", sep = ""))
colnames(asso.gr20.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

#asso.gr20.MLM$method <- "MLM"
#asso.gr20.MLMM$method <- "MLMM"
#asso.gr20.BLINK$method <- "BLINK"

#gwasdata.gr20 <- rbind(asso.gr20.MLM, asso.gr20.MLMM, asso.gr20.BLINK)
#gwasdata.gr20 <- arrange(gwasdata.gr20, Chr, Pos)

#For growth rate at 25 C
#asso.gr25.MLM <- read.csv(paste(assodata, "gr25/GAPIT.Association.GWAS_Results.MLM.gr25.csv", sep = ""))
#asso.gr25.MLMM <- read.csv(paste(assodata, "gr25/GAPIT.Association.GWAS_Results.MLMM.gr25.csv", sep = ""))
asso.gr25.BLINK <- read.csv(paste(assodata, "gr25/GAPIT.Association.GWAS_Results.BLINK.gr25.csv", sep = ""))
colnames(asso.gr25.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

#asso.gr25.MLM$method <- "MLM"
#asso.gr25.MLMM$method <- "MLMM"
#asso.gr25.BLINK$method <- "BLINK"

#gwasdata.gr25 <- rbind(asso.gr25.MLM, asso.gr25.MLMM, asso.gr25.BLINK)
#gwasdata.gr25 <- arrange(gwasdata.gr25, Chr, Pos)

#For growth rate at 30 C
#asso.gr30.MLM <- read.csv(paste(assodata, "gr30/GAPIT.Association.GWAS_Results.MLM.gr30.csv", sep = ""))
#asso.gr30.MLMM <- read.csv(paste(assodata, "gr30/GAPIT.Association.GWAS_Results.MLMM.gr30.csv", sep = ""))
asso.gr30.BLINK <- read.csv(paste(assodata, "gr30/GAPIT.Association.GWAS_Results.BLINK.gr30.csv", sep = ""))
colnames(asso.gr30.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

#asso.gr30.MLM$method <- "MLM"
#asso.gr30.MLMM$method <- "MLMM"
#asso.gr30.BLINK$method <- "BLINK"

#gwasdata.gr30 <- rbind(asso.gr30.MLM, asso.gr30.MLMM, asso.gr30.BLINK)
#gwasdata.gr30 <- arrange(gwasdata.gr30, Chr, Pos)

#For growth rate at 35 C
#asso.gr35.MLM <- read.csv(paste(assodata, "gr35/GAPIT.Association.GWAS_Results.MLM.gr35.csv", sep = ""))
#asso.gr35.MLMM <- read.csv(paste(assodata, "gr35/GAPIT.Association.GWAS_Results.MLMM.gr35.csv", sep = ""))
asso.gr35.BLINK <- read.csv(paste(assodata, "gr35/GAPIT.Association.GWAS_Results.BLINK.gr35.csv", sep = ""))
colnames(asso.gr35.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

#asso.gr35.MLM$method <- "MLM"
#asso.gr35.MLMM$method <- "MLMM"
#asso.gr35.BLINK$method <- "BLINK"

#gwasdata.gr35 <- rbind(asso.gr35.MLM, asso.gr35.MLMM, asso.gr35.BLINK)
#gwasdata.gr35 <- arrange(gwasdata.gr35, Chr, Pos)

#For growth rate at 37.5 C
#asso.gr375.MLM <- read.csv(paste(assodata, "gr375/GAPIT.Association.GWAS_Results.MLM.gr375.csv", sep = ""))
#asso.gr375.MLMM <- read.csv(paste(assodata, "gr375/GAPIT.Association.GWAS_Results.MLMM.gr375.csv", sep = ""))
asso.gr375.BLINK <- read.csv(paste(assodata, "gr375/GAPIT.Association.GWAS_Results.BLINK.gr375.csv", sep = ""))
colnames(asso.gr375.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

#asso.gr375.MLM$method <- "MLM"
#asso.gr375.MLMM$method <- "MLMM"
#asso.gr375.BLINK$method <- "BLINK"

#gwasdata.gr375 <- rbind(asso.gr375.MLM, asso.gr375.MLMM, asso.gr375.BLINK)
#gwasdata.gr375 <- arrange(gwasdata.gr375, Chr, Pos)

#For growth rate at 40 C
#asso.gr40.MLM <- read.csv(paste(assodata, "gr40/GAPIT.Association.GWAS_Results.MLM.gr40.csv", sep = ""))
#asso.gr40.MLMM <- read.csv(paste(assodata, "gr40/GAPIT.Association.GWAS_Results.MLMM.gr40.csv", sep = ""))
asso.gr40.BLINK <- read.csv(paste(assodata, "gr40/GAPIT.Association.GWAS_Results.BLINK.gr40.csv", sep = ""))
colnames(asso.gr40.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

#asso.gr40.MLM$method <- "MLM"
#asso.gr40.MLMM$method <- "MLMM"
#asso.gr40.BLINK$method <- "BLINK"

#gwasdata.gr40 <- rbind(asso.gr40.MLM, asso.gr40.MLMM, asso.gr40.BLINK)
#gwasdata.gr40 <- arrange(gwasdata.gr40, Chr, Pos)


### *** Looking at significant SNPs and Manhattan plot
### **** Associations in different temperatures
#holm <- adj.p(asso.gr20.MLM[,4], method = "Holm")
bonft <- bonf.threshold(asso.gr20.BLINK[,4], 0.01) #Significance threshold (This is the same for all temps, because of the number of SNPs is the same)

#Significant SNPs
signif.SNPs.gr20 <- filter(asso.gr20.BLINK, P.value < bonft)
hl.gr20 <- signif.SNPs.gr20[,1]

signif.SNPs.gr25 <- filter(asso.gr25.BLINK, P.value < bonft)
hl.gr25 <- signif.SNPs.gr25[,1]

signif.SNPs.gr30 <- filter(asso.gr30.BLINK, P.value < bonft)
hl.gr30 <- signif.SNPs.gr30[,1]

signif.SNPs.gr35 <- filter(asso.gr35.BLINK, P.value < bonft)
hl.gr35 <- signif.SNPs.gr35[,1]

signif.SNPs.gr375 <- filter(asso.gr375.BLINK, P.value < bonft)
hl.gr375 <- signif.SNPs.gr375[,1]

signif.SNPs.gr40 <- filter(asso.gr40.BLINK, P.value < bonft)
hl.gr40 <- signif.SNPs.gr40[,1]

hl.all <- c(as.character(signif.SNPs.gr20[,1]), as.character(signif.SNPs.gr25[,1]), as.character(signif.SNPs.gr30[,1]), as.character(signif.SNPs.gr35[,1]), as.character(signif.SNPs.gr375[,1]), as.character(signif.SNPs.gr40[,1]))

label20 <- expression(paste("Growth rate at 20 ", degree, "C"))
label25 <- expression(paste("Growth rate at 25 ", degree, "C"))
label30 <- expression(paste("Growth rate at 30 ", degree, "C"))
label35 <- expression(paste("Growth rate at 35 ", degree, "C"))
label375 <- expression(paste("Growth rate at 37.5 ", degree, "C"))
label40 <- expression(paste("Growth rate at 40 ", degree, "C"))

colnames(asso.gr20.BLINK)[c(2,3)] <- c("Chromosome", "Position")
colnames(asso.gr25.BLINK)[c(2,3)] <- c("Chromosome", "Position")
colnames(asso.gr30.BLINK)[c(2,3)] <- c("Chromosome", "Position")
colnames(asso.gr35.BLINK)[c(2,3)] <- c("Chromosome", "Position")
colnames(asso.gr375.BLINK)[c(2,3)] <- c("Chromosome", "Position")
colnames(asso.gr40.BLINK)[c(2,3)] <- c("Chromosome", "Position")

gr20.plot <- manhattan.plot(asso.gr20.BLINK, signift = bonft, mylabel = label20, hl = hl.all, myy = 18.5)

gr25.plot <- manhattan.plot(asso.gr25.BLINK, signift = bonft, mylabel = label25, hl = hl.all, myy = 18.5)

gr30.plot <- manhattan.plot(asso.gr30.BLINK, signift = bonft, mylabel = label30, hl = hl.all, myy = 18.5)

gr35.plot <- manhattan.plot(asso.gr35.BLINK, signift = bonft, mylabel = label35, hl = hl.all, myy = 18.5)

gr375.plot <- manhattan.plot(asso.gr375.BLINK, signift = bonft, mylabel = label375, hl = hl.all, myy = 18.5)

gr40.plot <- manhattan.plot(asso.gr40.BLINK, signift = bonft, mylabel = label40, hl = hl.all, myy = 18.5)

##Extract legend and add it to my manhattan plot
#manforlegend <- ggplot(gwasdata.gr20, aes(x = Pos, y = -log10(P.value), shape = method)) +
 # geom_point(alpha = 0.75, size = 2.5) +
 # scale_y_continuous(expand = c(0,0), limits = c(0, 16.5), breaks = c(seq(0,16,2))) +
 # theme(legend.title = element_blank(), legend.position = "top" )
#legend.man <- get_legend(manforlegend)

manplot.alltraits.hl <- plot_grid(gr20.plot, gr35.plot, gr25.plot, gr375.plot, gr30.plot, gr40.plot, align = "v", ncol = 2)

#final.man <- plot_grid(legend.man, manplot.alltraits.hl, ncol = 1, rel_heights = c(0.05, 1))

save_plot("./fig/mymanhattan.png",manplot.alltraits.hl, base_height = 8, base_width = 16)

###Save a summary of significant SNPs
signif.SNPs.gr20$Temperature <- 20
signif.SNPs.gr25$Temperature <- 25
signif.SNPs.gr30$Temperature <- 30
signif.SNPs.gr35$Temperature <- 35
signif.SNPs.gr375$Temperature <- 37.5
signif.SNPs.gr40$Temperature <- 40
signif.SNPs <- rbind(signif.SNPs.gr20, signif.SNPs.gr25, signif.SNPs.gr30, signif.SNPs.gr35, signif.SNPs.gr375, signif.SNPs.gr40)

##Extract all significant SNPs that were significant in one temperature from all temperatures
#Also scale the allelic effects by the mean growth rate in a given temperature (load genomeans from earlier)
tempmeans <- summarize(group_by(genomeans, Temp), avgr = mean(meangr, na.rm = T))

signf.by.temp.gr20 <- filter(asso.gr20.BLINK, SNP %in% signif.SNPs[,1])
signf.by.temp.gr20$Temperature <- rep(20, nrow(signf.by.temp.gr20))
signf.by.temp.gr20$Scaled.eff <- signf.by.temp.gr20$Effect / as.numeric(tempmeans[1,2]) #Scale allelic effect

signf.by.temp.gr25 <- filter(asso.gr25.BLINK, SNP %in% signif.SNPs[,1])
signf.by.temp.gr25$Temperature <- rep(25, nrow(signf.by.temp.gr25))
signf.by.temp.gr25$Scaled.eff <- signf.by.temp.gr25$Effect / as.numeric(tempmeans[2,2])

signf.by.temp.gr30 <- filter(asso.gr30.BLINK, SNP %in% signif.SNPs[,1])
signf.by.temp.gr30$Temperature <- rep(30, nrow(signf.by.temp.gr30))
signf.by.temp.gr30$Scaled.eff <- signf.by.temp.gr30$Effect / as.numeric(tempmeans[3,2])

signf.by.temp.gr35 <- filter(asso.gr35.BLINK, SNP %in% signif.SNPs[,1])
signf.by.temp.gr35$Temperature <- rep(35, nrow(signf.by.temp.gr35))
signf.by.temp.gr35$Scaled.eff <- signf.by.temp.gr35$Effect / as.numeric(tempmeans[4,2])

signf.by.temp.gr375 <- filter(asso.gr375.BLINK, SNP %in% signif.SNPs[,1])
signf.by.temp.gr375$Temperature <- rep(37.5, nrow(signf.by.temp.gr375))
signf.by.temp.gr375$Scaled.eff <- signf.by.temp.gr375$Effect / as.numeric(tempmeans[5,2])

signf.by.temp.gr40 <- filter(asso.gr40.BLINK, SNP %in% signif.SNPs[,1])
signf.by.temp.gr40$Temperature <- rep(40, nrow(signf.by.temp.gr40))
signf.by.temp.gr40$Scaled.eff <- signf.by.temp.gr40$Effect / as.numeric(tempmeans[6,2])

#Combine all together
signf.by.temp <- rbind(signf.by.temp.gr20, signf.by.temp.gr25, signf.by.temp.gr30, signf.by.temp.gr35, signf.by.temp.gr375, signf.by.temp.gr40)

#Change sign of allelic effects to be relative to the major allele (i.e if minor allele decreases growth rate, then allelic effect has a negative sign)
signvec <- c(1, -1, 1, -1, -1, -1, -1, -1, -1, 1, 1)
signf.by.temp$Effect <- signf.by.temp$Effect * rep(signvec, 6)
signf.by.temp$Scaled.eff <- signf.by.temp$Scaled.eff * rep(signvec, 6)

save(signif.SNPs, signf.by.temp, file = "./data/significantSNPs.RData")

#Load significant SNPs
load("./data/significantSNPs.RData")

##Making a table of the significant SNPs
#arrange(signif.SNPs, Temperature, Chr, Pos)

### Make a plot that shows the effects of all significant SNPs in all temperatures

#Make a label variable for the SNPs
signf.by.temp$snplabel <- paste("Chr ", signf.by.temp[,2], ": ", signf.by.temp[,3], sep = "")

templabel <-  expression(paste("Temperature (", degree, "C", ")"))
eff.plot <- ggplot(signf.by.temp, aes(x = Temperature, y = Effect)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, lty = "dashed") +
    xlab(templabel) +    
    facet_wrap( ~ snplabel)

save_plot("./fig/SNP_effect.pdf", eff.plot, base_height = 9, base_width = 12)

scaled.eff.plot <- ggplot(signf.by.temp, aes(x = Temperature, y = Scaled.eff)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, lty = "dashed") +
    #scale_y_continuous(
    ylab("Scaled allelic effect") +
    xlab(templabel) +    
    facet_wrap( ~ snplabel)

save_plot("./fig/SNP_effect_scaled.pdf", scaled.eff.plot, base_height = 9, base_width = 12)


### **** Allelic effects and allele frequency in natural populations

#Load significant SNPs
load("./data/significantSNPs.RData")

#I need to filter those SNPs in data.frame "signf.by.temp" to include only those that were statistically significant in in GWAS

sig.20 <- unique(filter(signif.SNPs, Temperature == 20)$SNP) #Names of SNPs that were significant in 20 C
sig.25 <- unique(filter(signif.SNPs, Temperature == 25)$SNP)
sig.30 <- unique(filter(signif.SNPs, Temperature == 30)$SNP)
sig.35 <- unique(filter(signif.SNPs, Temperature == 35)$SNP)
sig.375 <- unique(filter(signif.SNPs, Temperature == 37.5)$SNP)
sig.40 <- unique(filter(signif.SNPs, Temperature == 40)$SNP)


#To get significant SNPs from a given temperature
allelic.effects.20 <- filter(signf.by.temp, Temperature == 20 & SNP %in% sig.20)
allelic.effects.25 <- filter(signf.by.temp, Temperature == 25 & SNP %in% sig.25)
allelic.effects.30 <- filter(signf.by.temp, Temperature == 30 & SNP %in% sig.30)
allelic.effects.35 <- filter(signf.by.temp, Temperature == 35 & SNP %in% sig.35)
allelic.effects.375 <- filter(signf.by.temp, Temperature == 37.5 & SNP %in% sig.375)
allelic.effects.40 <- filter(signf.by.temp, Temperature == 40 & SNP %in% sig.40)

allelic.effects <- rbind(allelic.effects.20, allelic.effects.25, allelic.effects.30, allelic.effects.35, allelic.effects.375, allelic.effects.40)
allelic.effects$Temperature <- factor(allelic.effects$Temperature)

#Get all unique SNPs
SNPs.for.maf <- unique(allelic.effects$SNP)

#Need to calculate minor allele frequency only among natural pop strains
#Load the genotypes 
myG <- read.table("~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop_gwas_temp.txt", header = FALSE, stringsAsFactors = FALSE)
geno.names <- as.character(myG[1,-c(1:4)])

#Natural pops are first in the file
myG.natpop <- myG[,1:121]

#Get indexes where the interesting SNPs occur
snp.indices <- match(SNPs.for.maf, myG.natpop[,1])
#Select the SNPs I want for calculating minor allele frequency
myG.natpop <- myG.natpop[c(1,snp.indices),] 

#Calculate allele frequencies
af <- new.af(myG.natpop[-1,-c(1:11)])
MAF <- apply(af, MARGIN = 1, min)

data.frame(SNP = myG.natpop[,1][-1], MAF = MAF)

maf.result <- data.frame(SNP = myG.natpop[,1][-1], MAF = MAF)

#Since some SNPs occur multiple times in the allelic effects
for(i in 1:nrow(allelic.effects)) {
    index <- which(allelic.effects$SNP[i] == maf.result$SNP)
    allelic.effects$MAF[i] <- maf.result$MAF[index]
}

#Allelic effects with the MAF from natural populations
save(allelic.effects, file = "./data/allelic_effects.RData")

#Can load allelic effects again here...
load(file = "./data/allelic_effects.RData")

### Make a plot, with the raw data and a log-log plot with fit.

#Then plotting the allelic effects
#Untransformed data
ae1 <- ggplot(allelic.effects, aes(x = MAF, y = abs(Scaled.eff))) +
    geom_point(data = allelic.effects, aes(x = MAF, y = abs(Scaled.eff), colour = Temperature)) +
    #geom_smooth() +  
    ylab("Absolute scaled allelic effect") +
    xlab("Minor allele frequency") +
    scale_y_continuous(breaks = seq(0, 0.4, 0.1)) +
    theme(legend.position = "none")    

#Transformed data with 
ae2 <- ggplot(allelic.effects, aes(x = log(MAF), y = log(abs(Scaled.eff)))) +
    geom_point(data = allelic.effects, aes(x = log(MAF), y = log(abs(Scaled.eff)), colour = Temperature)) +
    geom_smooth(data = allelic.effects, aes(x = log(MAF), y = log(abs(Scaled.eff))), method = "lm") +
    ylab("log(absolute scaled allelic effect)") +
    xlab("log(minor allele frequency)")    

#Combine plots
ae.plot <- plot_grid(ae1, ae2, labels = c("A", "B"), rel_widths = c(1,1.3))

save_plot("./fig/allefplot.pdf", ae.plot, base_height = 4, base_width=10)


fit <- lm(log(abs(Scaled.eff)) ~ log(MAF), data = allelic.effects)
summary(fit)

### **** Get the types of signicant SNPs (which AA change etc.)

##Using Variant Effect Predictor tool for this
#http://fungi.ensembl.org/Neurospora_crassa/Tools/VEP
#First need to format the data accordingly

#Load the genotype data
myG <- read.table("~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop_gwas_temp.txt", header = FALSE, stringsAsFactors = FALSE)

#Make sure that significant SNPs are loaded

variants <- filter(myG, V1 %in% signif.SNPs[,1])[,1:4]
#Need to reformat this to for Variant effect predictor
#Format is:
#chromosome (as roman numerals (ensembl fungi for N. crassa)
#position
#.
#allele1
#allele2
#. . .
chrs <- recode(variants[,3], '1' = "I", '2' = "II", '3' = "III", '4' = "IV", '5' = "V", '6' = "VI", '7' = "VII")
alleles <- str_split_fixed(variants[,2], "", 2)
pos <- variants[,4]

#Formatting the SNPs for variant effect predictor
forVEP <- cbind(chrs, pos, rep(".", nrow(variants)), alleles[,1], alleles[,2], rep(".", nrow(variants)), rep(".", nrow(variants)),rep(".", nrow(variants)))

write.table(forVEP, file = "./data/SNPs_forVEP.txt", quote = FALSE, sep = " ", row.names = F, col.names = F)

### **** Looking at associations genotype-phenotype

#library(data.table)

#load the genotypic means
genomeans <- read.csv(paste(datafolder, "genomeans.csv", sep = ""), header = T, sep = ",")

#Load the genotypes 
myG <- read.table("~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop_gwas_temp.txt", header = FALSE, stringsAsFactors = FALSE)
geno.names <- as.character(myG[1,-c(1:4)])

#Need to filter again for SNP with low minor allele counts
#ma.counts <- calc.alleles(as.matrix(myG[,12:ncol(myG)])) #Calculate minor allele counts (better than freq since totals are not the same for all)
#check <- ma.counts > 5
#check[1] <- T ## Set the first row as true so that genotype names etc are not removed
#myG <- myG[check,] #Filter SNPs with low minor allele counts

#Load significant SNPs
load("./data/significantSNPs.RData")

#Write a function that gets genotypes from SNP
get.genotypes <- function(phenot, SNP) {
    snpname <- SNP[2,1]
    phenot$mysnp <- rep(0, nrow(phenot))
    SNP <- SNP[,-c(1:11)] #Drop unused columns
    genonames <- as.character(SNP[1,]) #Names of genotypes
    for(i in 1:nrow(phenot)) { #Loop over genotypes and get genotypes
        current.genot <- as.character(phenot[i,1])
        if((current.genot %in% genonames) == TRUE) { #Is the current genotype in genonames?
            index <- which(current.genot == genonames) #Check index
            phenot$mysnp[i] <- SNP[2,index] #Store the genotype
        } else { phenot$mysnp[i] <- NA }
    } #Done looping over all genotypes
    colnames(phenot)[5] <- snpname
    phenot[,5] <- replace(phenot[,5], phenot[,5] == "N", NA) #Replace 'N' with NA
    return(phenot)
}
    
#Load genome annotations and do some formatting
annotations <- read.table("~/Genomics/Neurospora/reference/neurospora_crassa_or74a_12_transcripts.gtf", header = F, sep = "\t")
colnames(annotations) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
annotations$Chr <- str_replace(as.character(annotations$seqname), "Supercontig_12.", "")
temp <- strsplit(as.character(annotations$attribute), ";", fixed = T)
temp <- unname(sapply(temp, '[[', 1))
annotations$GeneID <- str_replace(temp, "gene_id ", "")

#Load the PCAs for covariates
PCAs <- read.table("./data/GAPIT.Genotype.PCA.csv", header = T, sep = ",")

sortorder <- match(koe[,1], PCAs[,1])
PCAs[sortorder,]
cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)       

### Making a plot to view possible segregation of SNPs
### For snp649355 
SNP <- rbind(myG[1,], filter(myG, V1 == "snp649355"))
phenot <- filter(genomeans, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)
#Make the factor for faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "D"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family D", nrow(filter(koe, Family == "D"))), rep("Family G", nrow(filter(koe, Family == "G")))))
#fin <- filter(fin, is.na(fin[,5]) == FALSE)
#modelAll <- lm(data = filter(fin, fam == "All"), meangr ~ PC1 + PC2 + PC3 + PC4)
#modelO <- lm(data = filter(fin, fam == "Natural populations"), meangr ~  PC1 + PC2 + PC3 + PC4)
#modelO2 <- lm(data = filter(fin, fam == "Natural populations"), meangr ~ snp649355)        
#modelD <- lm(data = filter(fin, fam == "Family D"), meangr ~  snp649355 + PC1 + PC2 + PC3 + PC4)
#modelG <- lm(data = filter(fin, fam == "Family G"), meangr ~  snp649355 + PC1 + PC2 + PC3 + PC4)        
#fin$predgr <- c(predict(modelAll), predict(modelO), predict(modelD), predict(modelG))        
        
ssnp1 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp649355)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 1:7214188") + 
    ylab("Growth rate (mm / h) at 20 °C") +
    facet_wrap( . ~ fam, ncol = 4)

 ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = predgr, x = snp649355)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 1:7214188") + 
    ylab("Predicted growth rate (mm / h) at 20 °C") +
    facet_wrap( . ~ fam, ncol = 4)        

hist(filter(fin, fam == "All")$meangr)
hist(predict(modelAll))
hist(fitted(modelAll))        
        
#How about significance and LD?
testzoom <- plot.GWAS.res(asso.gr20.BLINK, annotations, 1, 7210000, 7220000) #NCU00720 (tca-17) at 30 C        

## For snp1598356
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1598356"))
phenot <- filter(genomeans, Temp == 25)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "D"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family D", nrow(filter(koe, Family == "D")))))

ssnp2 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp1598356)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 3:3119653") + 
    ylab("Growth rate (mm / h) at 25 °C") +
    facet_wrap( . ~ fam, ncol = 3)

## For snp2713600
SNP <- rbind(myG[1,], filter(myG, V1 == "snp2713600"))
phenot <- filter(genomeans, Temp == 25)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#Make the factor for faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "D"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family D", nrow(filter(koe, Family == "D"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp3 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp2713600)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 5:2747660") + 
    ylab("Growth rate (mm / h) at 25 °C") +
    facet_wrap( . ~ fam, ncol = 4)

## For snp3414038
SNP <- rbind(myG[1,], filter(myG, V1 == "snp3414038"))
phenot <- filter(genomeans, Temp == 25)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp4 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp3414038)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 6:3669489") + 
    ylab("Growth rate (mm / h) at 25 °C") +
    facet_wrap( . ~ fam, ncol = 3)


## For snp1574452
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1574452"))
phenot <- filter(genomeans, Temp == 30)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O")))))

ssnp5 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp1574452)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 3:2881603") + 
    ylab("Growth rate (mm / h) at 30 °C") +
    facet_wrap( . ~ fam, ncol = 2)

## For  snp3407080
SNP <- rbind(myG[1,], filter(myG, V1 == "snp3407080"))
phenot <- filter(genomeans, Temp == 30)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp6 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp3407080)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 6:3600716") + 
    ylab("Growth rate (mm / h) at 30 °C") +
    facet_wrap( . ~ fam, ncol = 3)

## For snp1519029
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1519029"))
phenot <- filter(genomeans, Temp == 35)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp7 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp1519029)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 3:2179382") + 
    ylab("Growth rate (mm / h) at 35 °C") +
    facet_wrap( . ~ fam, ncol = 3)

## For snp3407080
SNP <- rbind(myG[1,], filter(myG, V1 == "snp3407080"))
phenot <- filter(genomeans, Temp == 35)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp8 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp3407080)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 6:3600716") + 
    ylab("Growth rate (mm / h) at 35 °C") +
    facet_wrap( . ~ fam, ncol = 3)

## For snp1991930
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1991930"))
phenot <- filter(genomeans, Temp == 37.5)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting       
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "C"), filter(koe, Family == "D"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family C", nrow(filter(koe, Family == "C"))), rep("Family D", nrow(filter(koe, Family == "D"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp9 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp1991930)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 4:1399273") + 
    ylab("Growth rate (mm / h) at 37.5 °C") +
    facet_wrap( . ~ fam, ncol = 5)

## For snp3406826
SNP <- rbind(myG[1,], filter(myG, V1 == "snp3406826"))
phenot <- filter(genomeans, Temp == 37.5)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))

ssnp10 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp3406826)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 6:3599026") + 
    ylab("Growth rate (mm / h) at 37.5 °C") +
    facet_wrap( . ~ fam, ncol = 3)

##For snp1232635
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1232635"))
phenot <- filter(genomeans, Temp == 40)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "B"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family B", nrow(filter(koe, Family == "B")))))

ssnp11 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp1232635)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 2:3743231") + 
    ylab("Growth rate (mm / h) at 40 °C") +
    facet_wrap( . ~ fam, ncol = 3)

##For snp2709206
SNP <- rbind(myG[1,], filter(myG, V1 == "snp2709206"))
phenot <- filter(genomeans, Temp == 40)
koe <- get.genotypes(phenot, SNP)
koe[,5] <- factor(koe[,5])
#faceting        
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "B"), filter(koe, Family == "D"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family B", nrow(filter(koe, Family == "B"))), rep("Family D", nrow(filter(koe, Family == "D"))), rep("Family G", nrow(filter(koe, Family == "G")))))
        
ssnp12 <- ggplot(filter(fin, is.na(fin[,5]) == FALSE), aes(y = meangr, x = snp2709206)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    xlab("Chr 5:2708042") + 
    ylab("Growth rate (mm / h) at 40 °C") +
    facet_wrap( . ~ fam, ncol = 5)

##Making a final plot
test <- plot_grid(ssnp1, ssnp2, ssnp3, ssnp4, ssnp5, ssnp6, ssnp7, ssnp8, ssnp9, ssnp10, ssnp11, ssnp12, nrow = 12)
save_plot(filename = "./fig/phenoplot.pdf", test, base_height = 36)

#Separate plots for segragation figures
save_plot(filename = "./fig/seg20.pdf", ssnp1, base_height = 3.71, base_width = 3*4)

row25.1 <- plot_grid(ssnp2, NULL, rel_widths = c(3,1))
row25.2 <- plot_grid(ssnp3)
row25.3 <- plot_grid(ssnp4, NULL, rel_widths = c(3,1))
test25 <- plot_grid(row25.1, row25.2, row25.3, nrow = 3)

save_plot(filename = "./fig/seg25.pdf", test25, base_height = 3.71*3, base_width = 3*4)

row30.1 <- plot_grid(ssnp5, NULL, NULL, rel_widths = c(2,1,1), ncol = 3)
row30.2 <- plot_grid(ssnp6, NULL, rel_widths = c(3,1))
test30 <- plot_grid(row30.1, row30.2, nrow = 2)
save_plot(filename = "./fig/seg30.pdf", test30, base_height = 3.71*2, base_width = 3*4)

row35.1 <- plot_grid(ssnp7, NULL, rel_widths = c(3,1))
row35.2 <- plot_grid(ssnp8, NULL, rel_widths = c(3,1))
test35 <- plot_grid(row35.1, row35.2, nrow = 2)
save_plot(filename = "./fig/seg35.pdf", test35, base_height = 3.71*2, base_width = 3*4)

row37.1 <- plot_grid(ssnp9)
row37.2 <- plot_grid(ssnp10, NULL, NULL, rel_widths = c(3,1,1), ncol = 3)
test37 <- plot_grid(row37.1, row37.2, nrow = 2)
save_plot(filename = "./fig/seg37.pdf", test37, base_height = 3.71*2, base_width = 3*4)

row40.1 <- plot_grid(ssnp11, NULL, NULL, rel_widths = c(3,1,1), ncol = 3)
row40.2 <- plot_grid(ssnp12)
test40 <- plot_grid(row40.1, row40.2, nrow = 2)
save_plot(filename = "./fig/seg40.pdf", test40, base_height = 3.71*2, base_width = 3*4)

#Make plots by temperature or something...

### **** Check marker spacing and LD

## Marker spacing
## Can load data from association results.

asso.gr20.BLINK <- read.csv(paste(assodata, "gr20/GAPIT.Association.GWAS_Results.BLINK.gr20.csv", sep = ""))

chr1 <- filter(asso.gr20.BLINK, Chr == 1)
chr2 <- filter(asso.gr20.BLINK, Chr == 2)
chr3 <- filter(asso.gr20.BLINK, Chr == 3)
chr4 <- filter(asso.gr20.BLINK, Chr == 4)
chr5 <- filter(asso.gr20.BLINK, Chr == 5)
chr6 <- filter(asso.gr20.BLINK, Chr == 6)
chr7 <- filter(asso.gr20.BLINK, Chr == 7)

#Calculate median distance between two consecutive SNPs in the dataset. Does that fall between LD?
#Distance between two consecutive SNPs (in bp)
SNP.dist <- function(data) {
    resvec <- rep(0, nrow(data)-1)
    for(i in 1:(nrow(data) - 1)) { resvec[i] <- data$Pos[i+1] - data$Pos[i] }
    return(resvec)
}
        

dist.chr1 <- SNP.dist(chr1)
dist.chr2 <- SNP.dist(chr2)
dist.chr3 <- SNP.dist(chr3)
dist.chr4 <- SNP.dist(chr4)
dist.chr5 <- SNP.dist(chr5)
dist.chr6 <- SNP.dist(chr6)
dist.chr7 <- SNP.dist(chr7)

#Combine all distances together
SNP.distances <- c(dist.chr1, dist.chr2, dist.chr3, dist.chr4, dist.chr5, dist.chr6, dist.chr7)
chrs <- c(rep(1, length(dist.chr1)), rep(2, length(dist.chr2)), rep(3, length(dist.chr3)), rep(4, length(dist.chr4)), rep(5, length(dist.chr5)), rep(6, length(dist.chr6)), rep(7, length(dist.chr7)) )

summary(SNP.distances)
quantile(SNP.distances, probs = c(0.025, 0.5, 0.975))
#95% of consecutive SNP distances are within 139 bp

quantile(SNP.distances, probs = c(0.025, 0.5, 0.99))
#99% of consecutive distances are within 215 bp

#Make a dataframe for plotting
distances <- data.frame(dist = SNP.distances, chr = chrs)

inset <- ggplot(filter(distances, dist > 1000), aes(x = dist/1000)) +
    geom_histogram(fill = "white", col = "black") +
    scale_y_continuous(expand = c(0,0)) +
    #scale_x_continuous(limits = c(-1,1000)) +
    xlab("Distance between SNPs (kb)") +
    ylab("Count")


disthist <- ggplot(filter(distances, dist < 300), aes(x = dist)) +
    geom_histogram(fill = "white", col = "black") +
    scale_y_continuous(expand = c(0,0)) +
    #scale_x_continuous(limits = c(-1,1000)) +
    xlab("Distance between SNPs (bp)") +
    ylab("Count")

histplot <- ggdraw(disthist) +
    draw_plot(inset, 0.45, 0.45, 0.5, 0.5)


### LD along the chromosomes

################ Testing locally ###################

##First get SNPs only for natpop
#Load the genotypes 
myG <- read.table("~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop_all_filtered.txt", header = FALSE, stringsAsFactors = FALSE)
geno.names <- as.character(myG[1,-c(1:4)])

#Natural pops are first in the file
myG.natpop <- myG[,1:129]



#Need to run chromosome by chromosome
myG.chr1 <- filter(myG.natpop, V3 == 1)
myG.chr2 <- filter(myG.natpop, V3 == 2)
myG.chr3 <- filter(myG.natpop, V3 == 3)
myG.chr4 <- filter(myG.natpop, V3 == 4)
myG.chr5 <- filter(myG.natpop, V3 == 5)
myG.chr6 <- filter(myG.natpop, V3 == 6)
myG.chr7 <- filter(myG.natpop, V3 == 7)

### Testing ###
small <- myG.natpop[300:400,]

#Minor alleles
calc.alleles <- function(aineisto) {
    results <- rep(0, dim(aineisto)[1])
    for(i in 1:dim(aineisto)[1]) {
        count <- table(as.character(aineisto[i,]), exclude = "N")
        results[i] <- min(count)
    }
    return(results)
}

test4 <- calc.alleles(as.matrix(small[,12:129]))

check <- test4 > 3
small <- small[check,]

test4 <- calc.alleles(as.matrix(myG.chr3[,12:129]))
check <- test4 > 3
myG.chr3 <- myG.chr3[check,]

#Calculate pairwise LD for all SNPs
#Assumes that only SNP genotypes are present in the data
test <- calc.pw.LD(small[12:129])

#This function calculates pairwise physical distance of two SNPs along a chromosome
#input is a vector physical coordinates along the chromosome
pw.phys.dist <- function(aineisto) {
        aineisto <- as.numeric(aineisto) #Make sure distances are numeric
        #Number of SNPs
        nsnp <- length(aineisto)
        if(nsnp < 2) {stop(cat("Need more than 1 SNP to calculate LD!")) }

        #Initialize results matrix
        res.mat <- matrix(rep(0, nsnp*nsnp), ncol = nsnp)

        #Calculate distance
        for(i in 1:(nsnp-1)) {
            for(j in (i+1):nsnp) {
                res.mat[i,j] <- aineisto[j] - aineisto[i] #Distance between SNPs
            }
        }

        return(res.mat)
    }


calculate.LD <- function(P1, P2) {
    #P1 <- aineisto[i,] #Genotypes of first locus
    #P2 <- aineisto[j,] #Genotypes of second locus
    P12 <- paste(P1, P2, sep = "") #Haplotypes of both
    index <- !grepl("N", P12) #Those cases with complete haplotype data
    ntotal <- sum(index)      #Count haplotypes with no missing data
    n1 <- sort(unique(as.character(P1)))[1] #These need to be sorted
    n2 <- sort(unique(as.character(P2)))[1]
    label <- paste(n1, n2, sep = "") #Get the string for haplotype 11
    nhaplo <- sum(grepl(label, P12)) #Count how many haplotype 11

    #Haplotype frequency
    hapfreq <- nhaplo / ntotal

    #Calculate allele frequencies for those inds with complete haplotype data
    P1 <- P1[index]
    P2 <- P2[index]
    af.p1 <- new.af(P1)
    af.p2 <- new.af(P2)

    D <- hapfreq - af.p1[1]*af.p2[1] #Calculate D
    r2 <- (D^2)/(af.p1[1]*(1-af.p1[1])*af.p2[1]*(1-af.p2[1])) #Calculate LD
    #Store LD
    if((is.na(r2) == FALSE) & (r2 < Inf)) { result <- r2 } else {result <- NA }
    return(result)
}

#This function samples pairs of SNPs at random (input should be a chromosome)
#and calculates LD and physical distance between the SNPs
#npairs = number of SNP pairs
sample.SNPs.LD <- function(aineisto, npairs) {
    nsnp <- nrow(aineisto) #Number of SNPs in the dataset

    #initialize results
    res.mat <- matrix(rep(0, npairs*2), ncol = 2)
    colnames(res.mat) <- c("dist", "LD")

    #Sample SNP indices
    snp.indices <- matrix(sample(1:nsnp, size = npairs*2, replace = FALSE), ncol = 2)

    for(i in 1:nrow(snp.indices)) {
        #Take a pair on SNPs
        snp1 <- aineisto[snp.indices[i,1],]
        snp2 <- aineisto[snp.indices[i,2],]
        
        #Calculate physical distance
        res.mat[i,1] <- abs(as.numeric(snp1[4]) - as.numeric(snp2[4]))

        #Calculate LD
        res.mat[i,2] <- calculate.LD(snp1[12:129], snp2[12:129])
    }
    return(res.mat)
}
         


test2 <- pw.phys.dist(small[,4])

test.data <- data.frame(r2 = test[upper.tri(test)], dist = test2[upper.tri(test2)])

test[upper.tri(test)]

test2[upper.tri(test2)]

ggplot(test.data, aes(x = dist, y = r2)) +
    geom_point()

koe <- sample.SNPs.LD(myG.chr3, 10000)
ggplot(data.frame(koe), aes(x = dist, y = LD)) +
    geom_point()

#Sampling SNPs

#Select ten SNPsets spaced evenly along the chromosome
snpsets <- round(seq(100, nrow(myG.chr3)-101, length.out = 10), 0)

#Sample SNP sets from chromosome, leave space at the chr end
#snpsets <- sample(1:(nrow(myG.chr3)-100), size = 10, replace = FALSE)
foo <- list(0)
reslist <- rep(foo, 10)

for(i in 1:10) { #Loop over the 10 SNP sets
    sampledsnp <- snpsets[i]
    myG.SNP <- myG.chr3[sampledsnp:(sampledsnp+99),] #Take 100 SNPs

    LD <- calc.pw.LD(myG.SNP[12:129]) #Calculate pairwise LD for each SNP
    distbp <- pw.phys.dist(myG.SNP[,4]) #Pairwise physical distances for each SNP

    #Make a dataframe of the results
    resmat <- data.frame(r2 = LD[upper.tri(LD)], dist = distbp[upper.tri(distbp)])

    #Store the results in the list
    reslist[[i]] <- resmat
}

#Combine the results of different SNP sets
finalresults <- do.call("rbind", reslist)


    
### Okay, sample SNP and calculate LD within 10 kb

#100 SNPs is about 6 kb
#100 * 100 SNPs is 10000 pairwise comparisons
#Sample 10 sets from each chromosome -> 100 000 pairwise comparisons
#Across the seven chromosomes -> 700 000 pairwise comparisons



########################################################

### Load LD calculation results done on the cluster
### This is data for course grained LD over long distances


load("./data/pw_LD_chr1.RData")
LD.chr1 <- LD.data

load("./data/pw_LD_chr2.RData")
LD.chr2 <- LD.data

load("./data/pw_LD_chr3.RData")
LD.chr3 <- LD.data

load("./data/pw_LD_chr4.RData")
LD.chr4 <- LD.data

load("./data/pw_LD_chr5.RData")
LD.chr5 <- LD.data

load("./data/pw_LD_chr6.RData")
LD.chr6 <- LD.data

load("./data/pw_LD_chr7.RData")
LD.chr7 <- LD.data

LD.data.all <- rbind(LD.chr1, LD.chr2, LD.chr3, LD.chr4, LD.chr5, LD.chr6, LD.chr7)

LD.data.all$bin <- cut(LD.data.all$dist/1000000, breaks = seq(0,10,0.01)) #Make bins

#Calculate mean and lower and upper quantiles for LD across bins
bindata <- summarise(group_by(LD.data.all, bin), meanLD = mean(r2, na.rm = T), upper = quantile(r2, probs = 0.975, na.rm = T), lower = quantile(r2, probs = 0.025, na.rm = T))

#Make a numeric bin vector
##extract first element from list (TODO)
bin1 <- strsplit(gsub("\\(", "", as.character(bindata$bin)), split = ",")
bin1 <- as.numeric(unname(sapply(bin1, '[[', 1)))
bindata$numbin <- bin1

fitLD <- nls(r2 ~ r0*exp(k*dist), data = LD.data.all, start = list(r0 = 1, k = -0.001))
predLD <- predict(fitLD, list(dist = seq(0, 9668405, length.out = 10000)))
predictionLD <- data.frame(r2 = predLD, dist = seq(0, 9668405, length.out = 10000))

### Making a plot
LD.longrange <- ggplot(LD.data.all, aes(x = dist/1000000, y = r2)) +
    geom_hex() +
    scale_fill_viridis_c() +
    #geom_smooth(col = "red") +
    #geom_line(data = bindata, aes(x = numbin, y = meanLD), col = "red") +
    geom_line(data = predictionLD, aes(y = r2, x = dist / 1000000), col = "red") +
    ylab(TeX("$r^2$")) +
    xlab("Distance (Mb)")

save_plot("./fig/LD_longrange.pdf", LD.longrange)

ggplot(LD.data.all, aes(x = bin, y = r2)) +
    geom_boxplot()


### Data for LD on a more fine scale

load("./data/pw_LD_fine_chr1.RData")
LDfine.chr1 <- finalresults

load("./data/pw_LD_fine_chr2.RData")
LDfine.chr2 <- finalresults

load("./data/pw_LD_fine_chr3.RData")
LDfine.chr3 <- finalresults

load("./data/pw_LD_fine_chr4.RData")
LDfine.chr4 <- finalresults

load("./data/pw_LD_fine_chr5.RData")
LDfine.chr5 <- finalresults

load("./data/pw_LD_fine_chr6.RData")
LDfine.chr6 <- finalresults

load("./data/pw_LD_fine_chr7.RData")
LDfine.chr7 <- finalresults

LD.fine.all <- rbind(LDfine.chr1, LDfine.chr2, LDfine.chr3, LDfine.chr4, LDfine.chr5, LDfine.chr6, LDfine.chr7)

#Exponential decay
fit1 <- nls(r2 ~ r0*exp(k*dist), data = LD.fine.all, start = list(r0 = 1, k = -0.001))
pred1 <- predict(fit1, list(dist = seq(0,23000,1)))
prediction1 <- data.frame(r2 = pred1, dist = 1:23001)

#Logistic decay
#fit2 <- nls(r2 ~ c / (1 + a*exp(-b*dist)), data = LD.fine.all, start = list(c =8.761141e-01, a = 9.119125e+00, b = -5.091930e-06), trace = T, nls.control(printEval = T, minFactor = 1/60024))

#fit2 <- nls(r2 ~ SSlogis(dist, Asym, xmid, scal), data = LD.fine.all)
### Model does not fit the data very well...

### Making a plot
LD.fine <- ggplot(LD.fine.all, aes(x = dist, y = r2)) +
    geom_hex() +
    scale_fill_viridis_c() +
    geom_line(data = prediction1, lwd = 1.5, col = "red") +
    #geom_line(data = bindata, aes(x = numbin, y = meanLD), col = "red") +
    ylab(TeX("$r^2$")) +
    xlab("Distance (bp)")


##Making a plot with

full.LD.plot <- plot_grid(LD.fine, histplot, ncol = 2, labels = c("A", "B"))

save_plot("./fig/fullLD.pdf", full.LD.plot, base_height=3.71, base_width = 3.71*1.618*2)

### **** Association for Topt and extremes conditioned on max growth rate

### ***** Association and manhattan

#Folder for association data for a given trait (Change when necessary)
assodata <- "~/Genomics/Neurospora/association/"

#For optimum temperature, growth at 20 C conditioned on maximum growth rate, and growth at 40 C conditioned on maximum growth rate
asso.Topt.BLINK <- read.csv(paste(assodata, "Topt/GAPIT.Association.GWAS_Results.BLINK.Topt.csv", sep = ""))
colnames(asso.Topt.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

asso.gr20.umax.BLINK <- read.csv(paste(assodata, "gr20_umax/GAPIT.Association.GWAS_Results.BLINK.gr20.csv", sep = ""))
colnames(asso.gr20.umax.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed

asso.gr40.umax.BLINK <- read.csv(paste(assodata, "gr40_umax/GAPIT.Association.GWAS_Results.BLINK.gr40.csv", sep = ""))
colnames(asso.gr40.umax.BLINK) <- c("SNP", "Chromosome", "Position", "P.value", "nobs", "Effect", "H.B.P.Value", "MAF") #Colnames are wrong in GAPIT output and they need to be fixed


#holm <- adj.p(asso.gr20.MLM[,4], method = "Holm")
bonft <- bonf.threshold(asso.Topt.BLINK[,4], 0.01) #Significance threshold (This is the same for all temps, because of the number of SNPs is the same)

#Significant SNPs
signif.SNPs.Topt <- filter(asso.Topt.BLINK, P.value < bonft)
hl.Topt <- signif.SNPs.Topt[,1]

signif.SNPs.gr20.umax <- filter(asso.gr20.umax.BLINK, P.value < bonft)
hl.gr20.umax <- signif.SNPs.gr20.umax[,1]

signif.SNPs.gr40.umax <- filter(asso.gr40.umax.BLINK, P.value < bonft)
hl.gr40.umax <- signif.SNPs.gr40.umax[,1]

#Labels
labelTopt <- expression(paste("Optimum temperature (", degree, "C)"))
labelgr20.umax <- expression(paste("Growth at 20", degree, "C", " with maximum growth rate as covariate"))
labelgr40.umax <- expression(paste("Growth at 40", degree, "C", " with maximum growth rate as covariate"))

#Highlight significant SNPs
hl.tpc <- c(as.character(signif.SNPs.Topt[,1]), as.character(signif.SNPs.gr20.umax[,1]), as.character(signif.SNPs.gr40.umax[,1]))


#Manhattan plots
Topt.plot <- manhattan.plot(asso.Topt.BLINK, signift = bonft, mylabel = labelTopt, hl = hl.tpc, myy = 22.5)

gr20.umax.plot <- manhattan.plot(asso.gr20.umax.BLINK, signift = bonft, mylabel = labelgr20.umax, hl = hl.tpc, myy = 22.5)

gr40.umax.plot <- manhattan.plot(asso.gr40.umax.BLINK, signift = bonft, mylabel = labelgr40.umax, hl = hl.tpc, myy = 22.5)

#Making the final manhattan plot

manplot.tpc.hl <- plot_grid(Topt.plot, gr20.umax.plot, gr40.umax.plot, align = "v", ncol = 1)

#final.man <- plot_grid(legend.man, manplot.alltraits.hl, ncol = 1, rel_heights = c(0.05, 1))

save_plot("./fig/TPC_manhattan.png",manplot.tpc.hl, base_height = 8, base_width = 8)

#Save significant TPC SNPs

signif.SNPs.TPC <- rbind(signif.SNPs.Topt, signif.SNPs.gr20.umax, signif.SNPs.gr40.umax)
signif.SNPs.TPC$trait <- c("Topt", rep("gr20.umax", 7), rep("gr40.umax", 3))

save(signif.SNPs.TPC, file = "./data/significantSNPs_TPC.RData")



### *****  Information about TPC SNPs

##Using Variant Effect Predictor tool for this
#http://fungi.ensembl.org/Neurospora_crassa/Tools/VEP
#First need to format the data accordingly

#load("./data/significantSNPs_TPC.RData")

#Load the genotype data
myG <- read.table("~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop_gwas_temp.txt", header = FALSE, stringsAsFactors = FALSE)

#Make sure that significant SNPs are loaded

variants <- filter(myG, V1 %in% signif.SNPs.TPC[,1])[,1:4]
#Need to reformat this to for Variant effect predictor
#Format is:
#chromosome (as roman numerals (ensembl fungi for N. crassa)
#position
#.
#allele1
#allele2
#. . .
chrs <- recode(variants[,3], '1' = "I", '2' = "II", '3' = "III", '4' = "IV", '5' = "V", '6' = "VI", '7' = "VII")
alleles <- str_split_fixed(variants[,2], "", 2)
pos <- variants[,4]

#Formatting the SNPs for variant effect predictor
forVEP <- cbind(chrs, pos, rep(".", nrow(variants)), alleles[,1], alleles[,2], rep(".", nrow(variants)), rep(".", nrow(variants)),rep(".", nrow(variants)))

write.table(forVEP, file = "./data/SNPs_forVEP_TPC.txt", quote = FALSE, sep = " ", row.names = F, col.names = F)

### ***** Allelic effects of TPC shape SNPs

load(file = "./data/significantSNPs_TPC.RData")

splinedata <- read.csv("./data/splineresults.csv", header = T, sep = ",")
genomeans <- read.csv(paste(datafolder, "genomeans.csv", sep = ""), header = T, sep = ",")

#Combine phenotypic data
splinedata <- arrange(splinedata, Genotype)
geno.gr20 <- arrange(filter(genomeans, Temp == 20), Genot)
geno.gr20 <- filter(geno.gr20, Genot != "XG2")
gr20.data <- cbind(geno.gr20, splinedata)
geno.gr40 <- arrange(filter(genomeans, Temp == 40), Genot)
geno.gr40 <- filter(geno.gr40, Genot != "XG2")
gr40.data <- cbind(geno.gr40, splinedata)

signif.by.trait.gr20.umax <- filter(asso.gr20.umax.BLINK, SNP %in% signif.SNPs.TPC[,1])
signif.by.trait.gr20.umax$trait <- "gr20.umax"
signif.by.trait.gr20.umax$Scaled.eff <- signif.by.trait.gr20.umax$Effect / mean(gr20.data$meangr)

signif.by.trait.gr40.umax <- filter(asso.gr40.umax.BLINK, SNP %in% signif.SNPs.TPC[,1])
signif.by.trait.gr40.umax$trait <- "gr40.umax"
signif.by.trait.gr40.umax$Scaled.eff <- signif.by.trait.gr40.umax$Effect / mean(gr40.data$meangr)

signif.by.trait.Topt <- filter(asso.Topt.BLINK, SNP %in% signif.SNPs.TPC[,1])
signif.by.trait.Topt$trait <- "Topt"
signif.by.trait.Topt$Scaled.eff <- signif.by.trait.Topt$Effect / mean(splinedata$Topt)

#Combine all together
signif.by.trait <- rbind(signif.by.trait.gr20.umax, signif.by.trait.gr40.umax, signif.by.trait.Topt)

#Make a label variable for the SNPs
signif.by.trait$snplabel <- paste("Chr ", signif.by.trait[,2], ": ", signif.by.trait[,3], sep = "")

scaled.eff.trait.plot <- ggplot(signif.by.trait, aes(x = trait, y = Scaled.eff)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, lty = "dashed") +
    #scale_y_continuous(
    ylab("Scaled allelic effect") +
    xlab("Trait") +    
    facet_wrap( ~ snplabel)

#Change sign of allelic effects to be relative to the major allele (i.e if minor allele decreases growth rate, then allelic effect has a negative sign)
#signvec <- c(1, -1, 1, -1, -1, -1, -1, -1, -1, 1, 1)
#signf.by.temp$Effect <- signf.by.temp$Effect * rep(signvec, 6)
#signf.by.temp$Scaled.eff <- signf.by.temp$Scaled.eff * rep(signvec, 6)

#save(signif.SNPs, signf.by.temp, file = "./data/significantSNPs.RData")

### Checking the direction of allelic effects
#Load the genotypes 
myG <- read.table("~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop_gwas_temp.txt", header = FALSE, stringsAsFactors = FALSE)
geno.names <- as.character(myG[1,-c(1:4)])

#Write a function that gets genotypes from SNP
get.genotypes <- function(phenot, SNP) {
    snpname <- SNP[2,1]
    phenot$mysnp <- rep(0, nrow(phenot))
    SNP <- SNP[,-c(1:11)] #Drop unused columns
    genonames <- as.character(SNP[1,]) #Names of genotypes
    for(i in 1:nrow(phenot)) { #Loop over genotypes and get genotypes
        current.genot <- as.character(phenot[i,1])
        if((current.genot %in% genonames) == TRUE) { #Is the current genotype in genonames?
            index <- which(current.genot == genonames) #Check index
            phenot$mysnp[i] <- SNP[2,index] #Store the genotype
        } else { phenot$mysnp[i] <- NA }
    } #Done looping over all genotypes
    colind <- ncol(phenot)
    colnames(phenot)[colind] <- snpname
    phenot[,colind] <- replace(phenot[,colind], phenot[,colind] == "N", NA) #Replace 'N' with NA
    return(phenot)
}

#Load the PCAs for covariates
PCAs <- read.table("./data/GAPIT.Genotype.PCA.csv", header = T, sep = ",")

### Making a plot to view possible segregation of SNPs
### For snp169331 
SNP <- rbind(myG[1,], filter(myG, V1 == "snp169331"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp169331 = rep(koe$snp169331, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp169331) == F)

##Need to check which families there is segregation

model.all <- lm(data = koe, meangr ~ snp169331 + umax + PC1 + PC2 + PC3 + PC4 + snp169331*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp169331 + umax + PC1 + PC2 + PC3 + PC4 + snp169331*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.C <- lm(data = filter(koe, Family == "C"), meangr ~ snp169331 + umax + PC1 + PC2 + PC3 + PC4 + snp169331*Temp)
summary(model.C)
pred.C <- predict(model.C)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp169331 + umax + PC1 + PC2 + PC3 + PC4 + snp169331*Temp)
summary(model.G)
pred.G <- predict(model.G)

#Make the factor for faceting
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "C"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))),rep("Family C", nrow(filter(koe, Family == "C"))), rep("Family G", nrow(filter(koe, Family == "G")))))
fin$predicted <- c(pred.all, pred.nat, pred.C, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp169331, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

my.xlab <- expression(paste("Temperature (", degree,"C)"))
chr1.1 <- ggplot(fin, aes(y = predicted, x = snp169331, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp169331, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp169331, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red"), name = my.xlab) +
    xlab("Chr 1: 1 638 766") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 4) +
    theme(legend.position = "none")    

#Get the legend for later
test <- chr1.1 + theme(legend.position = "top", legend.justification = "center")
legenda <- get_legend(test)

#####################
#snp189210
SNP <- rbind(myG[1,], filter(myG, V1 == "snp189210"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp189210 = rep(koe$snp189210, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp189210) == F)

model.all <- lm(data = koe, meangr ~ snp189210 + umax + PC1 + PC2 + PC3 + PC4 + snp189210*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp189210 + umax + PC1 + PC2 + PC3 + PC4 + snp189210*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp189210 + umax + PC1 + PC2 + PC3 + PC4 + snp189210*Temp)
summary(model.G)
pred.G <- predict(model.G)

#Facetig
fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))
fin$predicted <- c(pred.all, pred.nat, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp189210, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr1.2 <- ggplot(fin, aes(y = predicted, x = snp189210, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp189210, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp189210, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 1: 1 934 979") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")     


##################
#snp564926
SNP <- rbind(myG[1,], filter(myG, V1 == "snp564926"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp564926 = rep(koe$snp564926, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp564926) == F)

model.all <- lm(data = koe, meangr ~ snp564926 + umax + PC1 + PC2 + PC3 + PC4 + snp564926*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp564926 + umax + PC1 + PC2 + PC3 + PC4 + snp564926*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp564926 + umax + PC1 + PC2 + PC3 + PC4 + snp564926*Temp)
summary(model.G)
pred.G <- predict(model.G)


fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))
fin$predicted <- c(pred.all, pred.nat, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp564926, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr1.3 <- ggplot(fin, aes(y = predicted, x = snp564926, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp564926, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp564926, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 1: 6 331 236") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")     

#Interaction is significant, but the main effect is not... weak evidence...

###################################

###################################
#snp1228644
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1228644"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp1228644 = rep(koe$snp1228644, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp1228644) == F)

model.all <- lm(data = koe, meangr ~ snp1228644 + umax + PC1 + PC2 + PC3 + PC4 + snp1228644*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp1228644 + umax + PC1 + PC2 + PC3 + PC4 + snp1228644*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.B <- lm(data = filter(koe, Family == "B"), meangr ~ snp1228644 + umax + PC1 + PC2 + PC3 + PC4 + snp1228644*Temp)
summary(model.B)
pred.B <- predict(model.B)

fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "B"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family B", nrow(filter(koe, Family == "B")))))
fin$predicted <- c(pred.all, pred.nat, pred.B)

pred.means <- summarise(group_by(fin, Temp, snp1228644, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr2.1 <- ggplot(fin, aes(y = predicted, x = snp1228644, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp1228644, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp1228644, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 2: 3 718 388") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")    

 #There is no trade-off in natural populations

#############################

#snp1543447
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1543447"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp1543447 = rep(koe$snp1543447, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp1543447) == F)

model.all <- lm(data = koe, meangr ~ snp1543447 + umax + PC1 + PC2 + PC3 + PC4 + snp1543447*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp1543447 + umax + PC1 + PC2 + PC3 + PC4 + snp1543447*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.C <- lm(data = filter(koe, Family == "C"), meangr ~ snp1543447 + umax + PC1 + PC2 + PC3 + PC4 + snp1543447*Temp)
summary(model.C)
pred.C <- predict(model.C)

model.D <- lm(data = filter(koe, Family == "D"), meangr ~ snp1543447 + umax + PC1 + PC2 + PC3 + PC4 + snp1543447*Temp)
summary(model.D)
pred.D <- predict(model.D)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp1543447 + umax + PC1 + PC2 + PC3 + PC4 + snp1543447*Temp)
summary(model.G)
pred.G <- predict(model.G)

fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "C"), filter(koe, Family == "D"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family C", nrow(filter(koe, Family == "C"))), rep("Family D", nrow(filter(koe, Family == "D"))), rep("Family G", nrow(filter(koe, Family == "G")))  ))
fin$predicted <- c(pred.all, pred.nat, pred.C, pred.D, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp1543447, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr3.1 <- ggplot(fin, aes(y = predicted, x = snp1543447, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp1543447, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp1543447, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 3: 2 484 109") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 5) +
    theme(legend.position = "none")    

#There is no trade-off in natural populations

###########################################

############
#snp1648785
SNP <- rbind(myG[1,], filter(myG, V1 == "snp1648785"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp1648785 = rep(koe$snp1648785, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp1648785) == F)

model.all <- lm(data = koe, meangr ~ snp1648785 + umax + PC1 + PC2 + PC3 + PC4 + snp1648785*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp1648785 + umax + PC1 + PC2 + PC3 + PC4 + snp1648785*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp1648785 + umax + PC1 + PC2 + PC3 + PC4 + snp1648785*Temp)
summary(model.G)
pred.G <- predict(model.G)


fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G")))))
fin$predicted <- c(pred.all, pred.nat, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp1648785, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr3.2 <- ggplot(fin, aes(y = predicted, x = snp1648785, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp1648785, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp1648785, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 3: 3 616 832") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")    
 #There is no trade-off in natural populations

####################################################


###########
#snp2683552
SNP <- rbind(myG[1,], filter(myG, V1 == "snp2683552"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp2683552 = rep(koe$snp2683552, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp2683552) == F)

model.all <- lm(data = koe, meangr ~ snp2683552 + umax + PC1 + PC2 + PC3 + PC4 + snp2683552*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp2683552 + umax + PC1 + PC2 + PC3 + PC4 + snp2683552*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.B <- lm(data = filter(koe, Family == "B"), meangr ~ snp2683552 + umax + PC1 + PC2 + PC3 + PC4 + snp2683552*Temp)
summary(model.B)
pred.B <- predict(model.B)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp2683552 + umax + PC1 + PC2 + PC3 + PC4 + snp2683552*Temp)
summary(model.G)
pred.G <- predict(model.G)


fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "B"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family B", nrow(filter(koe, Family == "B"))), rep("Family G", nrow(filter(koe, Family == "G"))) ))
fin$predicted <- c(pred.all, pred.nat, pred.B, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp2683552, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr5.1 <- ggplot(fin, aes(y = predicted, x = snp2683552, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp2683552, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp2683552, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 5: 2 439 327") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 4) +
    theme(legend.position = "none")    

 #There is no trade-off in natural populations

#########################################33

############
#snp2715914
SNP <- rbind(myG[1,], filter(myG, V1 == "snp2715914"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp2715914 = rep(koe$snp2715914, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp2715914) == F)

model.all <- lm(data = koe, meangr ~ snp2715914 + umax + PC1 + PC2 + PC3 + PC4 + snp2715914*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp2715914 + umax + PC1 + PC2 + PC3 + PC4 + snp2715914*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.G <- lm(data = filter(koe, Family == "G"), meangr ~ snp2715914 + umax + PC1 + PC2 + PC3 + PC4 + snp2715914*Temp)
summary(model.G)
pred.G <- predict(model.G)

fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "G"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family G", nrow(filter(koe, Family == "G"))) ))
fin$predicted <- c(pred.all, pred.nat, pred.G)

pred.means <- summarise(group_by(fin, Temp, snp2715914, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr5.2 <- ggplot(fin, aes(y = predicted, x = snp2715914, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp2715914, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp2715914, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 5: 2 767 945") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")    


###########################################

##
#snp3555731
SNP <- rbind(myG[1,], filter(myG, V1 == "snp3555731"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp3555731 = rep(koe$snp3555731, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp3555731) == F)

model.all <- lm(data = koe, meangr ~ snp3555731 + umax + PC1 + PC2 + PC3 + PC4 + snp3555731*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp3555731 + umax + PC1 + PC2 + PC3 + PC4 + snp3555731*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

model.D <- lm(data = filter(koe, Family == "D"), meangr ~ snp3555731 + umax + PC1 + PC2 + PC3 + PC4 + snp3555731*Temp)
summary(model.D)
pred.D <- predict(model.D)

fin <- rbind(koe, filter(koe, Family == "O"), filter(koe, Family == "D"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))), rep("Family D", nrow(filter(koe, Family == "D"))) ))
fin$predicted <- c(pred.all, pred.nat, pred.D)

pred.means <- summarise(group_by(fin, Temp, snp3555731, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr7.1 <- ggplot(fin, aes(y = predicted, x = snp3555731, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp3555731, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp3555731, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 7: 720 256") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")    



#############################3

###
#snp3833472
SNP <- rbind(myG[1,], filter(myG, V1 == "snp3833472"))
phenot <- filter(gr20.data, Temp == 20)
koe <- get.genotypes(phenot, SNP)
koe$gr40 <- gr40.data[,4]
koe[,8] <- factor(koe[,8])
koe <- cbind(koe, PCAs[match(koe[,1], PCAs[,1]),-1]) #One liner to get PCAs (sorting is correct)

koe <- data.frame(Genot = rep(koe$Genot,2), meangr = c(koe$meangr, koe$gr40), Temp = factor(c(rep(20, nrow(koe)), rep(40, nrow(koe)))), snp3833472 = rep(koe$snp3833472, 2), umax = rep(koe$umax, 2), PC1 = rep(koe$PC1,2), PC2 = rep(koe$PC2, 2), PC3 = rep(koe$PC3, 2), PC4 = rep(koe$PC4), Family = rep(koe$Family,2) )
koe <- filter(koe, is.na(snp3833472) == F)

model.all <- lm(data = koe, meangr ~ snp3833472 + umax + PC1 + PC2 + PC3 + PC4 + snp3833472*Temp)
summary(model.all)
pred.all <- predict(model.all)

model.nat <- lm(data = filter(koe, Family == "O"), meangr ~ snp3833472 + umax + PC1 + PC2 + PC3 + PC4 + snp3833472*Temp)
summary(model.nat)
pred.nat <- predict(model.nat)

#Check which families are polymorphic
#koe[,c(1,4,10)]

fin <- rbind(koe, filter(koe, Family == "O"))
fin$fam <- factor(c(rep("All", nrow(koe)), rep("Natural populations", nrow(filter(koe, Family == "O"))) ))
fin$predicted <- c(pred.all, pred.nat)

pred.means <- summarise(group_by(fin, Temp, snp3833472, fam), meanpred = mean(predicted), npred = n(), sdpred = sd(predicted, na.rm = T), lower = mean(predicted) - 1.96*(sd(predicted)/sqrt(nobs(predicted))), upper = mean(predicted) + 1.96*(sd(predicted)/sqrt(nobs(predicted))) )

chr7.2 <- ggplot(fin, aes(y = predicted, x = snp3833472, fill = Temp)) +
    geom_point(position = position_jitterdodge(dodge.width=0.75), alpha = 0.25, shape = 21) +
    geom_line(data=pred.means, aes(y = meanpred, x = snp3833472, group = Temp), position = position_dodge(width=0.75)) +    
    geom_pointrange(data=pred.means, aes(y = meanpred, x = snp3833472, fill = Temp, ymin = lower, ymax = upper), position = position_dodge(width=0.75), shape = 21, size = 0.75) +
    scale_fill_manual(values = c("blue", "red")) +
    xlab("Chr 7: 3 358 427") +
    ylab("Predicted growth") +
    facet_wrap(. ~ fam, ncol = 3) +
    theme(legend.position = "none")    


###########################################

#Making the final figures
row1 <- plot_grid(chr1.2, NULL, rel_widths = c(3,1))
row2 <- plot_grid(chr1.3, NULL, rel_widths = c(3,1))

row3 <- plot_grid(chr7.1, NULL, rel_widths = c(3,1))
#row4 <- plot_grid(chr7.2, NULL, rel_widths = c(2,1), nrow=1)

#testi <- plot_grid(legenda, chr1.1, row1, row2, nrow = 4, rel_heights = c(0.1,1,1,1))

#save_plot(filename = "./fig/TPC_tradeoff.pdf", testi, base_height = 8.1, base_width = 8)

#testi2 <- plot_grid(legenda, chr7.1, row4, nrow = 3, rel_heights = c(0.1,1,1,1))

#save_plot(filename = "./fig/TPC_tradeoff2.pdf", testi2, base_height = 5+(1/3)+0.1, base_width = 6)


row4 <- plot_grid(chr7.2, NULL, NULL, rel_widths = c(2,1,1), nrow=1)
test3 <- plot_grid(legenda, chr1.1, row1, row2, row3, row4, nrow = 6, rel_heights = c(0.1,1,1,1,1,1))

save_plot(filename = "./fig/TPC_tradeoff3.pdf", test3, base_height = 12.1, base_width = 8)


rivi1 <- plot_grid(chr2.1, NULL, NULL, rel_widths = c(3,1,1), nrow = 1)
rivi2 <- chr3.1
rivi3 <- plot_grid(chr3.2, NULL, NULL, rel_widths = c(3,1,1), nrow = 1)
rivi4 <- plot_grid(chr5.1, NULL, rel_widths = c(4,1), nrow = 1)
rivi5 <- plot_grid(chr5.2, NULL, NULL, rel_widths = c(3,1,1), nrow = 1)

kuva3 <- plot_grid(legenda, rivi1, rivi2, rivi3, rivi4, rivi5, nrow = 6, rel_heights = c(0.1,1,1,1,1,1))

save_plot(filename = "./fig/TPC_tradeoff_non1.pdf", kuva3, base_height = 12.1, base_width = 10)

### * Illumina seq reports

load(paste0(datafolder, "seq_report.RData"))

### * PCA variance explained figure

PCAvar <- read.csv(paste0(datafolder, "GAPIT.Genotype.PCA_eigenvalues.csv"), header = T, dec = ".")

#Calculate percentage variance explained
PCAvar$percentage <- (PCAvar$x / sum(PCAvar$x)) * 100
PCAvar$component <- factor(1:nrow(PCAvar))
PCAvar <- PCAvar[1:10,] #No need to show all the tiny components
PCAvar$component <- factor(PCAvar$component)

PCAplot <- ggplot(PCAvar, aes(y = percentage, x = component)) +
    geom_bar(stat = "identity", colour = "black", fill = "grey") +
    ylab("Variance explained (%)") +
    xlab("Principle component") +
    scale_y_continuous(expand = c(0,0), limits = c(0,12))


save_plot(filename = "./fig/PCAplot.pdf", PCAplot)
