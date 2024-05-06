#This file contains scripts used in the association / power analysis of the E. coli project

#Load libraries required by GAPIT
library(multtest)
library(gplots)
##Need to install dependency snpStats first...
#BiocManager::install(c("snpStats"))
#install.packages("ape")
#install.packages("EMMREML")
library(LDheatmap) #install.packages("/home/ililkron/Documents/tutkijatohtori/association/LDheatmap_0.99-5.tar.gz", repos=NULL, type="source")
library(genetics)
library(compiler) #this library is already installed in R
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("gapit_functions_20160216.txt")
source("http://zzlab.net/GAPIT/emma.txt")



#This function calculates allele frequencies
#Have to deal with loci with multiple alleles somehow...(still to do) now using lapply and just picking the smallest freq
#Sometimes I want to return a list and sometimes a matrix
#The allele.freq function now takes an argument should the output be a list or a matrix
#output takes arguments "list" or "matrix", defaults to "list"
allele.freq <- function(aineisto, output = "list")
{
	counts <- apply(aineisto, MARGIN = 2, table)
        if(output == "list") {
	if(is.list(counts) == TRUE) { freqs <- lapply(counts, prop.table) } else { freqs <- as.list(as.data.frame(apply(counts, MARGIN = 2, prop.table))) }
    }
        if(output == "matrix") {
            if(is.list(counts) == TRUE) { freqs <- data.frame(lapply(counts, prop.table)) } else { freqs <- apply(counts, MARGIN = 2, prop.table) }
        }
        

	return(freqs)
}

#This function calculates minor allele frequencies from the output of the previous function
minor.allele.freq <- function(allele.freq)
{
	minor.af <- lapply(allele.freq, min)
	minor.af <- unlist(minor.af)
	return(minor.af)
}

##Just to make sure I'm not breaking any stuff make new allele frequency calc. functions
##This function assumes that individuals are on columns and loci are in rows
##Should work with the hapmap data format (remember to drop all but genotype columns)
new.af <- function(aineisto)
    {
        counts <- apply(aineisto, MARGIN = 1, table, exclude = "N")
        if(length(counts) == 1) {freqs <- 1 } else { #If monomorphic set frequency to 1
        freqs <- t(apply(counts, MARGIN = 2, prop.table)) }
        return(freqs)
    }

##Function for calculating pairwise LD between SNPs
##How to deal with missing data? -> Only complete cases are used. That is only those individuals with complete haplotype data are used in allele frequency calculations
calc.pw.LD <- function(aineisto)
    {
        #Number of SNPs
        nsnp <- nrow(aineisto)
        if(nsnp < 2) {stop(cat("Need more than 1 SNP to calculate LD!")) }

        #Initialize results matrix
        resmat <- matrix(rep(0, nsnp*nsnp), ncol = nsnp)

        #Calculate LD
        for(i in 1:(nsnp-1)) {
            for(j in (i+1):nsnp) {
                P1 <- aineisto[i,] #Genotypes of first locus
                P2 <- aineisto[j,] #Genotypes of second locus
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
                if((is.na(r2) == FALSE) & (r2 < Inf)) { resmat[i,j] <- r2 } else {resmat[i,j] <- NA }
            }
        }

        return(resmat)
    }
            
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


#This function counts the number of missing data for each locus
count.NA <- function(aineisto)
{
	counts <- apply(aineisto, MARGIN = 2, is.na)
	nas <- apply(counts, MARGIN = 2, sum)
	return(nas)
}

#Minor alleles
calc.alleles <- function(aineisto) {
    results <- rep(0, dim(aineisto)[1])
    for(i in 1:dim(aineisto)[1]) {
        count <- table(as.character(aineisto[i,]), exclude = "N")
        results[i] <- min(count)
    }
    return(results)
}


##Function used in the analysis of haplotypes



#Some function needed for power analysis

#Simulate QTNs
#I simulate the QTN genotypes since it is more straight forward to do that than to deal with missing data
### STOP THE PRESS!!!, of course we need to use the actual SNPs for the QTNs ###
#This needs to be modified, just have to deal with the missing data somehow...

#note that SNPs have have their colnames already in place before this is run
simulate.QTN.genotype <- function(SNP.data, nloc, maf, na.threshold)
{
	#Need to sample SNPs with a given maf
	allele.freqs <- allele.freq(SNP.data)
	minor.af <- minor.allele.freq(allele.freqs)

	#Some SNPs have quite a few NA's
	na.counts <- count.NA(SNP.data)
	
	#Select SNPs
	nind <- nrow(SNP.data) #Number of individuals
	na.prop <- na.counts / nind

	#Accepting only SNPs with "na.threshold % or more data
	na.ok <- na.prop <= na.threshold

	#Correct minor.af
	maf.ok <- minor.af >= (maf - 0.05) & minor.af <= (maf + 0.05)

	#Correct minor allele frequency and not too much missing data
	snp.selection <- na.ok == T & maf.ok == T
	#Check that there are enough SNPs that fulfill the criteria
	if(sum(snp.selection) < nloc) stop(cat("Not enough SNPs that fill criteria. Adjust selection parameters.", "\n"))

	candidate.snps <- SNP.data[,snp.selection] #These are the SNPs that fill the criteria
	sample.index <- sample(x = 1:ncol(candidate.snps), size = nloc, replace = FALSE)
	
	final.snps <- candidate.snps[,sample.index]
	
	#Store QTN colnames now!
	write.csv(colnames(final.snps), file = "causal_SNPs.txt", quote = FALSE, row.names = FALSE)
	
	#Converting the SNPs to numeric format for phenotype simulation
	res.mat <- matrix(rep(0, nloc*nind), ncol = nloc)

	for(i in 1:nloc)
	{
		res.mat[,i] <- as.numeric(final.snps[,i]) #Convert to numeric, minor allele not considered
		
	}

	return(res.mat)
}

#This function simulates phenotypes, needs as input the number of different genotypes, number of times each genotype is replicated,
# standard deviation among genotypes, standard deviation due to environmental effect (error), and grand mean of phenotype
simulate.phenotypes <- function(ngenot, nrep, sd.genot, sd.env, grand.mean)
{
	#Prepare data structure
	genotype.vec <- gl(ngenot, nrep) #vector of genotypes
	y <- rep(0, ngenot*nrep) #vector of phenotypes

	#Simulate phenotypic data
	genotypic.effect <- rnorm(n = ngenot, mean = 0, sd = sd.genot)
	environmental.effect <- rnorm(n = ngenot*nrep, mean = 0, sd = sd.env)

	for(i in 1:(ngenot*nrep))
	{
		y[i] <- grand.mean + genotypic.effect[genotype.vec[i]] + environmental.effect[i]
	}

	#Assemble the data
	aineisto <- as.data.frame(cbind(y, genotype.vec))
	aineisto[,2] <- factor(aineisto[,2])
	colnames(aineisto) <- c("Z", "Genotype")

	return(aineisto)
}


#This function simulates phenotypes but with QTNs
#I assume that each QTN has (nearly) identical allelic effect (allele frequencies differ a bit)
simulate.phenotypes.QTN <- function(ngenot, nrep, sd.genot, sd.env, grand.mean, SNP.data, nloc, maf, na.threshold)
{
	#Prepare data structure
	genotype.vec <- gl(ngenot, nrep) #vector of genotypes
	y <- rep(0, ngenot*nrep) #vector of phenotypes
	QTN.vec <- simulate.QTN.genotype(SNP.data = SNP.data, nloc = nloc, maf = maf, na.threshold = na.threshold) #matrix of QTN genotypes
	
	#Allele frequencies for QTNs
	QTN.freq <- allele.freq(QTN.vec, output = "matrix")

	### Simulate the phenotypic data ###
	
	var.genot <- sd.genot^2 #Genotypic variance
	#Calculating the allelic effects for each locus
	allelic.effects <- rep(0, nloc)
	#effects of minor allele = a, effects of major allele = 0
	#I assume that each locus makes an equal contribution to the genetic variance
	#sigma^2.g = sum 2p1p2a^2, thus a_i = +- sqrt( (sigma^2.g / nloc) / (2*maf*(1-maf)) )
	for(i in 1:nloc)
	{
		allelic.effects[i] <- sqrt( (var.genot/nloc) / (QTN.freq[1,i]*QTN.freq[2,i]) ) * sample(c(1,-1), size = 1)
	}	
	
	#Calculating genotypic effects by summing over all allelic effects
	genotypic.effect <- rep(0, ngenot)	
	for(i in 1:ngenot)
	{
		genotypic.effect[i] <- sum(allelic.effects[(QTN.vec[i,] == "1")], na.rm = TRUE)
	}
	
	#Generate environmental effects
	environmental.effect <- rnorm(n = ngenot*nrep, mean = 0, sd = sd.env)
	
	#Calculate phenotypes
	for(i in 1:(ngenot*nrep))
	{
		y[i] <- grand.mean + genotypic.effect[genotype.vec[i]] + environmental.effect[i]
	}

	#Assemble the data
	aineisto <- as.data.frame(cbind(y, genotype.vec))
	aineisto[,2] <- factor(aineisto[,2])
	colnames(aineisto) <- c("Z", "Genotype")

	#Could change this to return also the QTN identities since QTN.vec is no longer needed
	return(list(aineisto, QTN.vec)) #Returns a list with phenotype data in element [[1]] and genotype data in element [[2]]

}

#This function converts the SNP data into numeric format for GAPIT
#Afterward it maybe easier if minor allele is converted to 2 (but not necessary)
SNP.2.numeric <- function(SNP.data)
{
	nSNP <- ncol(SNP.data) #Number of columns
	nind <- nrow(SNP.data) #Number of individuals

	res.mat <- matrix(rep(0, nSNP*nind), ncol = nSNP)

	for(i in 1:nSNP)
	{
		converted <- as.numeric(SNP.data[,i]) #Convert to numeric, minor allele not considered
		res.mat[,i] <- replace(converted, converted == 1, 0) #Replace 1 with 0 and store
	}

	return(res.mat)

}

#Need to convert the SNP data into hapmap format, note: this function assumes that SNPs are already in rows, inds at cols
hapmap.format <- function(SNP.data)
{
	nind <- ncol(SNP.data) #number of individuals
	#Hapmap convention for missing genotype data is "N" for single bit and "NN" for double bit
	for(i in 1:nind)
	{
		na.index <- is.na(SNP.data[,i])
		SNP.data[,i] <- replace(SNP.data[,i], na.index, "N")
	}

	return(SNP.data)
}

#SNP positions
#This functions generates some random SNP positions for test purposes
#E.coli genome about 4.5 Mb = 4500000
generate.SNP.positions <- function(nSNP, chromosome.length)
{
	positions <- sample.int(n = chromosome.length, size = nSNP, replace = FALSE)
	return(positions)
}


#Check for triallelic SNPs
is.triallelic <- function(SNP.data)
{
	nSNP <- ncol(SNP.data)	
	res.mat <- rep(0, nSNP)
	
	for(i in 1:nSNP)
	{
		res.mat[i] <- nlevels(SNP.data[,i])
	}
	res.mat <- res.mat == 3
	return(res.mat)
}

#Check for triallelic SNPs in data where SNPs are not factors and in hapmap format
is.triallelic2 <- function(SNP.data)
{
	nSNP <- nrow(SNP.data)
	res.mat <- rep(0, nSNP)

	for(i in 1:nSNP)
	{
		res.mat[i] <- sum(names(table(unlist(SNP.data[i,]))) != "N") > 2
	}

	return(res.mat)
}

#This function flags all SNPs that are not biallelic (not polymorphic or are triallelic)
is.not.biallelic <- function(SNP.data)
{
	nSNP <- nrow(SNP.data)
	res.mat <- rep(0, nSNP)

	for(i in 1:nSNP)
	{
		res.mat[i] <- sum(names(table(unlist(SNP.data[i,]))) != "N") != 2
	}

	return(res.mat)
}


	
#This function shows the genotypes of the individuals for a certain SNP
#Input is the myG matrix
#Chromosome is in column 3, and position in column 4

get.genotypes <- function(Gdata, chr, pos)
    {
        res <- Gdata[Gdata[,3] == chr & Gdata[,4] == pos,] #Get correct row
        colnames(res) <- unlist(lapply(Gdata[1,], as.character))
        return(res)
    }


##This function reads the genotype names from hapmap and phenotype data,
##Compares them and makes genotype and phenotype variables with matching genotypes
##Assumes that there are two columns ("Genotype" and "Phenotype")
##Can drop extra columns with drop.col argument
##Phenotype column can also be renamed with phenotype.name argument
format.GAPIT.data.phenotype <- function(geno.names, pheno.names, phenodata, drop.col = NULL, phenotype.name = "trait1") {

    myY <- phenodata[pheno.names %in% geno.names,]
    if(is.null(drop.col) == FALSE) { myY <- myY[,-drop.col] } 
    colnames(myY) <- c("Taxa", phenotype.name)

    return(myY)
}

##This function reads the genotype names from hapmap and phenotype data,
##Compares them and makes a genotype dataset with matching genotypes
##Function assumes that data is in hapmap format with 11 columns before individuals start
format.GAPIT.data.genotype <- function(geno.names, pheno.names, genotypes) {

    genofilt <- geno.names %in% pheno.names #Get the columns that match
    myG <- genotypes[,c(rep(T,11),genofilt)] #Filter genotypes

    return(myG)
}


#This function calculates the bonferroni threshold
bonf.threshold <- function(pvalues, threshold) {
    #Bonferroni threshold
    #adjusted = pvalue * n
    #pvalue = adjusted / n
    bonft <- threshold / length(pvalues)  #Bonferroni significance threshold
    return(bonft)
}

#This function calculates the FDR significance threshold
#Assumes that pvalues are in ascending order
fdr.threshold <- function(fdr.pvalues, pvalues, threshold) {
    ##False discovery rate
    #set p-values in ascending order (from smallest to largest)
    #Benjamini-Hochberg procedure
    #Adjusted p-value is calculated as
    #FDRadjusted = (pvalue*m)/k, m = number of tests, k = rank of p-value (1 = smallest value, and so on)
    #BH threshold -> find the last significant test at 0.05 level -> threhold is between that p-value and the next nonsignificant one
    firstno <- which(fdr.pvalues > threshold)[1]
    lastyes <- firstno - 1
    #
    FDR.threshold <- (pvalues[firstno] + pvalues[lastyes])/2
    return(FDR.threshold)
}


##This makes a manhattan plot using ggplot2 commands
manhattan.plot <- function(gwasdata, signift, suggest, mylabel = NULL, hl = NULL, myy = 30.5) {
    ##Making a custom Manhattan plot using ggplot2
    ##Making cumulative base pair position
    nCHR <- length(unique(gwasdata$Chromosome))
    gwasdata$BPcum <- NA
    s <- 0
    nbp <- c()
    for (i in 1:nCHR) {
        nbp[i] <- max(gwasdata[gwasdata$Chromosome == i,]$Position)
        gwasdata[gwasdata$Chromosome == i,"BPcum"] <- gwasdata[gwasdata$Chromosome == i,"Position"] + s
        s <- s + nbp[i]
    }

    axis.set <- group_by(gwasdata, Chromosome)
    axis.set <- summarize(axis.set, center = (max(BPcum) + min(BPcum))/2)
    ylim <- abs(floor(log10(min(gwasdata$P.value)))) + 2
    bonft <- bonf.threshold(gwasdata$P.value, threshold = 0.01) #Calc bonferroni threshold
    fdrt <- fdr.threshold(gwasdata$FDR_Adjusted_P.values, gwasdata$P.value, threshold = 0.05) #Calc FDR threshold
    ylabel <- expression(paste(-"log"[10], "(p)", sep = ""))
    hl.bpcum <- NULL
    if(is.null(hl) == F) {
    #Lines to highlight candidates
    hl.ind <- match(hl, gwasdata[,1])
    hl.bpcum <- gwasdata[hl.ind, "BPcum"] }

ggplot(gwasdata, aes(x = BPcum, y = -log10(P.value), 
  color = as.factor(Chromosome))) + #, #size = -log10(P.value))) +
  geom_vline(xintercept = hl.bpcum, color = "grey40", linetype = "dotted") +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(bonft), color = "grey40") +
  geom_hline(yintercept = -log10(fdrt), color = "grey40", linetype = "dashed") +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, myy), breaks = c(seq(0,myy-0.5,4))) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  #scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = ylabel, subtitle = mylabel) +
  #ggtitle(label = mylabel) +     
  #theme_minimal() +
  theme( 
    legend.position = "none",
  #  panel.border = element_blank(),
  #  panel.grid.major.x = element_blank(),
  #  panel.grid.minor.x = element_blank(),
  #  axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )
    
}

##This is probably best done so that raise the 0-level and plot the gene models beneath
##Need some sort of automatic functio to plot the gene models

#This function returns midpoints and widths
box.midpoint <- function(dat) { dat[1] + (dat[2]-dat[1])/2 }
box.width <- function(dat) { dat[2]-dat[1] }  

##Function for making nice plots of associated SNPs for shorter regions
##Need to add genome wide significance and suggestive lines (to do!)
##Need to add functionality for plotting direction of transcription (i.e. strand) (to do!)
##Need to to able to deal with overlapping genes (to do!)

plot.GWAS.res <- function(GWASres, annodata, chr, start.pos, end.pos) {
    #First filter the date 
    smalldata <- filter(GWASres, Chromosome == chr, Position >= start.pos, Position <= end.pos)
    smallanno <- filter(annodata, Chr == chr, start >= start.pos, end <= end.pos)

    #Making a dataframe of all genes in this interval
    genes <- unique(smallanno$GeneID)
    genelist <- data.frame(matrix(rep(0, length(genes)*6), ncol = 6))
    colnames(genelist) <- c("start", "end", "mid", "ID", "strand", "ycoord")

    #Loop over all genes and build dataframe 
    for(i in 1:length(genes)) {
        tempdata <- filter(smallanno, GeneID == genes[i])
        if(tempdata$strand[1] == "+") {
        genelist[i,1:3] <- c(min(tempdata$start), max(tempdata$end), box.midpoint(c(min(tempdata$start), max(tempdata$end))))
        genelist$ycoord[i] <- -1 } else {  genelist[i,1:3] <- c(max(tempdata$end), min(tempdata$start), box.midpoint(c(min(tempdata$start), max(tempdata$end))))
        genelist$ycoord[i] <- -2 }                                
        genelist[i,4] <- unique(as.character(tempdata$GeneID))      
    }

    
    #Make y-axis label
    my.ylab <- expression(paste(-"log"[10], "(p)", sep = ""))

    ##start plotting, note that need to use two separate lapply commands for gene models and labels
    ggplot() +
        ##First plot gene models
        lapply(split(genelist, 1:nrow(genelist)), function(dat) {   
         geom_segment(data = dat, aes(x = start, xend = end, y = ycoord, yend = ycoord), size = 3, lineend = "round", linejoin = "mitre", arrow = arrow(length = unit(0.005, "inches"))) } ) +
        ##Then plot gene labels    
        lapply(split(genelist, 1:nrow(genelist)), function(dat) {    
        geom_text(data = dat, aes(x = mid, y = ycoord - 0.5, label = ID)) } ) +
        ##Plot association p-values    
        geom_point(data = smalldata, aes(x = Position, y = -log10(P.value))) +
        scale_y_continuous(expand = c(0,0), limits = c(-3, 12), breaks = c(0,seq(1:12))) +
        ylab(my.ylab) +
        xlab("Position (bp)")    

}
