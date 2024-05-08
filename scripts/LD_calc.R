### * Load libraries
library(dplyr)

### * Functions to calculate LD

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
### This function needs to be changed to have smaller numbers of SNPs at a time
### Figure a distance to where LD is very small, i.e 1000 SNPs downstream or so...
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

#This function calculates the minor allele count (better to use this than freq, since total is not the same for all)
calc.alleles <- function(aineisto) {
    results <- rep(0, dim(aineisto)[1])
    for(i in 1:dim(aineisto)[1]) {
        count <- table(as.character(aineisto[i,]), exclude = "N")
        results[i] <- min(count)
    }
    return(results)
}

#This function calculates LD between two SNPs
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

print("Functions loaded")

### * Doing the calculations by chromosome

##Usage: Rscript LD_calc.R <chr>

##--MAIN--##
args <- commandArgs(trailingOnly=TRUE)

cur.chr <- as.integer(args[1]) #Chromosome to analyze

##Load the SNP data
myG <- read.table("/projappl/project_2000350/gwas/data/Neuro_hapmap_natpop_all_filtered.txt", header = FALSE, stringsAsFactors = FALSE)

#Natural pops are first in the file
myG.natpop <- myG[,1:129]
rm(myG) #Delete this to free up memory

myG.chr <- filter(myG.natpop[-1,], V3 == cur.chr) #Take only the current chromosome, drop first row (names)

ma.counts <- calc.alleles(as.matrix(myG.chr[,12:129])) #Calculate minor allele counts (better than freq since totals are not the same for all)
check <- ma.counts > 3
myG.chr <- myG.chr[check,] #Filter SNPs with low minor allele counts

#pw.LD <- calc.pw.LD(myG.chr[,12:129]) #Calculate pairwise LD between all SNPs on a chromosome, and keep only genotype columns)

#pw.dist <- pw.phys.dist(myG.chr[,4]) #Calculate pairwise physical distances (bp) for all SNPs, position is column 4

LD <- sample.SNPs.LD(myG.chr, 10000) #Sample 10 000 SNP pairs and calculate LD between each

#Make the final dataframe
LD.data <- data.frame( r2 = LD[,2], dist = LD[,1], chr = rep(cur.chr, nrow(LD)) )

filename <- paste0("pw_LD_chr", cur.chr, ".RData") #Setting up filenames

save(LD.data, file = paste0("/projappl/project_2000350/gwas/results/", filename)) #Save the data


