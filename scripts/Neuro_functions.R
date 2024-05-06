#These are functions required for the analysis of Neurospora growth rate data
#Load this file before analyzing the data

## Functions needed ##

#Function to drop zeroes from raw data
#The raw data includes zeroes. Obviously for estimateing growth rate this is not what we want.
#However, these will be useful when estimating  the "lag phase" or time until germination. Note that this is a very crude measure of time until germination. Mainly intended for screening.

#Function to plot growth rate of a single replicate.
#This is very important for visually looking the data

#Many of the growth experiments will have different time points ( points on x-axis will be different)
#There has to be a function that associates the correct time points with the correct measurements
#Will be probably the easiest to add an index to all measurements that indicates which time points are used

#Function to calculate linear growth rates for all
#Should the first data point be included or not? Have to make a decision about this.

#Dealing with dataframes etc.

#Wrapper for the Bayesian analysis

######

#library(splines)
#library(car) #Functions for calculating type III sums of squares ;)
#library(plyr)
#library(rjags)
#library(ggplot2)
#library(lme4)
#library(xtable)
#theme_set(theme_bw(base_size=16))
#vplayout <- function(x,y) {viewport(layout.pos.row = x, layout.pos.col = y) } #Need this function for plotting

#Alternative palette
# The palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2")

#theme_set(theme_bw()) #Setting the BW ggplot theme
#theme_update(axis.text = element_text(size = 14), axis.title = element_text(size=14), legend.text = element_text(size = 12), legend.title = element_text(size = 12))

#library(grid)
#library(gridBase)
#library(scales)
#library(cowplot)

### Time indexes ###
#In the data there is a variable called tindex, this variable defines at which time points the measurements were taken
#There is some variation in the index among replicates because, sometimes it was not possible to measure different replicates at the same time

#The time indexes (indexii?) are

tindex.1 <- c(0, 8, 24, 30, 48, 56, 72, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.2 <- c(0, 8, 24, 32, 48, 55, 72, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.3 <- c(0, 8, 24, 32, 48, 56, 72, 79, 96, 104, 120, 128, NA, NA, NA)
tindex.4 <- c(0, 7, 24, 32, 48, 56, 72, 80, 96, 104, 120, 128, 144, 152, NA)
tindex.5 <- c(0, 8, 24, 32, 48, 56, 72, 80, 96, 104, 120, 128, 144, 152, NA)
tindex.6 <- c(0, 8, 24, 32, 49, 56, 72, 80, 96, 104, 120, 127, 144, 152, NA)
tindex.7 <- c(0, 8, 24, 32, 48, 56, 72, 79, 96, 104, 120, 128, 144, 152, NA)
tindex.8 <- c(0, 8, 25, 32, 49, 56, 72, 80, 96, 103, 120, 128, 144, 152, NA)
tindex.9 <- c(0, 8, 24, 31, 48, 56, 75, 96, 120, 144, NA, NA, NA, NA, NA)
tindex.10 <- c(0, 8, 24, 32, 49, 56, 73, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.11 <- c(0, 8, 24, 32, 48, 56, 72, 95, 120, 128, 144, 152, NA, NA, NA)
tindex.12 <- c(0, 8, 25, 32, 48, 56, 72, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.13 <- c(0, 6, 23, 30, 47, 55, 71, 79, 96, 104, 120, 128, NA, NA, NA)
tindex.14 <- c(0, 8, 24, 32, 48, 56, 72, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.15 <- c(0, 8, 24, 32, 48, 57, 72, 81, 96, 104, 120, 128, NA, NA, NA)
tindex.16 <- c(0, 10, 24, 32, 48, 55, 72, 80, 98, 104, 120, 128, NA, NA, NA)
tindex.17 <- c(0, 8, 22, 31, 47, 55, 71, 79, 95, 103, 120, 128, NA, NA, NA)
tindex.18 <- c(0, 8, 24, 32, 48, 56, 73, 81, 96, 104, 120, 128, NA, NA, NA)
tindex.19 <- c(0, 23, 46, 72, 80, 96, 120, 128, 145, 152, NA, NA, NA, NA, NA)
tindex.20 <- c(0, 8, 24, 32, 48, 56, 96, 104, 120, 128, NA, NA, NA, NA, NA)
tindex.21 <- c(0, 9, 24, 32, 48, 56, 72, 80, 96, 104, NA, NA, NA, NA, NA)
tindex.22 <- c(0, 7, 24, 32, 48, 56, 72, 80, 96, 104, NA, NA, NA, NA, NA)
tindex.23 <- c(0, 6, 22, 30, 48, 56, 71, 79, 95, 103, NA, NA, NA, NA, NA)
tindex.24 <- c(0, 8, 24, 32, 46, 54, 71, 79, 95, 105, NA, NA, NA, NA, NA)
tindex.25 <- c(0, 8, 24, 32, 48, 56, 72, 80, 96, 104, NA, NA, NA, NA, NA)
tindex.26 <- c(0, 8, 25, 49, 57, 73, 80, 96, 104, 120, NA, NA, NA, NA, NA)
tindex.27 <- c(0, 8, 24, 32, 48, 56, 73, 81, 97, 105, NA, NA, NA, NA, NA)
tindex.28 <- c(0, 8, 24, 32, 48, 56, 72, 80, 97, 104, NA, NA, NA, NA, NA)
tindex.29 <- c(0, 7, 23, 31, 47, 55, 71, 81, 95, 104, 119, 128, NA, NA, NA)
tindex.30 <- c(0, 8, 24, 32, 49, 56, 74, 80, 97, 105, 121, 129, NA, NA, NA)
tindex.31 <- c(0, 8, 24, 33, 48, 56, 71, 80, 94, 104, 120, 128, NA, NA, NA)
tindex.32 <- c(0, 8, 24, 32, 48, 56, 71, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.33 <- c(0, 6, 22, 30, 46, 54, 70, 78, 94, 102, 118, 142, NA, NA, NA)
tindex.34 <- c(0, 8, 24, 32, 48, 56, 72, 80, 97, 104, 120, 128, NA, NA, NA)
tindex.35 <- c(0, 8, 24, 32, 48, 56, 72, 81, 96, 104, 120, 128, NA, NA, NA)
tindex.36 <- c(0, 8, 24, 32, 48, 55, 72, 79, 96, 103, 120, 128, NA, NA, NA)
tindex.37 <- c(0, 4, 21, 29, 46, 52, 69, 77, 96, 102, 120, 128, NA, NA, NA)
tindex.38 <- c(0, 5, 24, 30, 48, 54, 72, 78, 96, 102, 120, 128, NA, NA, NA)
tindex.39 <- c(0, 8, 24, 34, 48, 56, 72, 80, 96, 104, 120, 128, NA, NA, NA)
tindex.40 <- c(0, 7, 23, 25, 27, 29, 31, 48, 55, 71, 79, 95, NA, NA, NA)
tindex.41 <- c(0, 7, 23, 31, 48, 49, 50, 51, 52, 53, 55, 71, NA, NA, NA)
tindex.42 <- c(0, 17, 18, 19, 20, 21, 22, 24, 40, 48, 64, 72, NA, NA, NA)
tindex.43 <- c(0, 16, 24, 41, 48, 64, 72, 90, 91, 92, 93, 94, 95, 96, 114)
tindex.44 <- c(0, 7, 24, 31, 48, 55, 72, 79, 105, 128, NA, NA, NA, NA, NA)
tindex.45 <- c(0, 7, 24, 31, 49, 55, 73, 79, 130, 138, NA, NA, NA, NA, NA)
tindex.46 <- c(0, 7, 25, 31, 48, 55, 72, 79, 97,  103, NA, NA, NA, NA, NA)
tindex.47 <- c(0, 9, 27, 33, 51, 57, 73, 79, 96, 104, NA, NA, NA, NA, NA)
tindex.48 <- c(0, 6, 24, 31, 48, 55, 73, 79, 96, 104, NA, NA, NA, NA, NA)
tindex.49 <- c(0, 7, 25, 32, 47, 54, 72, 79, 96, 104, NA, NA, NA, NA, NA)
tindex.50 <- c(0, 8, 25, 31, 48, 53, 72, 79, 96, 104, NA, NA, NA, NA, NA)
tindex.51 <- c(0, 8, 24, 25, 26, 27, 28, 29, 30, 31, 32, 48, 56, NA, NA)
tindex.52 <- c(0, 8, 25, 26, 27, 28, 29, 30, 31, 32, 33, 48, 56, NA, NA) 
tindex.53 <- c(0, 7, 23, 24, 25, 26, 27, 28, 29, 30, 31, 47, 55, NA, NA)
tindex.54 <- c(0, 8, 24, 25, 26, 27, 28, 29, 30, 31, 32, 50, 57, NA, NA)
tindex.55 <- c(0, 8, 24, 32, 48, 56, 72, 80, NA, NA, NA, NA, NA, NA, NA)
tindex.56 <- c(0, 6, 22, 30, 46, 54, 70, 78, NA, NA, NA, NA, NA, NA, NA)
tindex.57 <- c(0, 8, 26, 32, 50, 57, 72, 80, NA, NA, NA, NA, NA, NA, NA)
tindex.58 <- c(0, 7, 24, 31, 48, 55, 72, 79, NA, NA, NA, NA, NA, NA, NA)
tindex.59 <- c(0, 6, 24, 30, 48, 54, 72, 78, NA, NA, NA, NA, NA, NA, NA)
tindex.60 <- c(0, 7, 24, 30, 48, 54, 72, 78, NA, NA, NA, NA, NA, NA, NA)


#Since some time indexes are of different lengths  added NA's to the ends to make them all the same length


tindex <- rbind(tindex.1, tindex.2, tindex.3, tindex.4, tindex.5, tindex.6, tindex.7, tindex.8, tindex.9, tindex.10, tindex.11, tindex.12, tindex.13, tindex.14, tindex.15, tindex.16, tindex.17, tindex.18, tindex.19, tindex.20, tindex.21, tindex.22, tindex.23, tindex.24, tindex.25, tindex.26, tindex.27, tindex.28, tindex.29, tindex.30, tindex.31, tindex.32, tindex.33, tindex.34, tindex.35, tindex.36, tindex.37, tindex.38, tindex.39, tindex.40, tindex.41, tindex.42, tindex.43, tindex.44, tindex.45, tindex.46, tindex.47, tindex.48, tindex.49, tindex.50, tindex.51, tindex.52, tindex.53, tindex.54, tindex.55, tindex.56, tindex.57, tindex.58, tindex.59, tindex.60)

#temp.tindex <- rbind(tindex.1, tindex.2, tindex.3)
#row.names(temp.tindex) <- c("1", "2", "3")
#temp.tindex2 <- rbind(tindex.4, tindex.5, tindex.6)
#row.names(temp.tindex2) <- c("4", "5", "6")

##################

### Plotting one replicate ###

#The idea is to plot the growth rate of one replicate to view the data

#Need to have a funtion, that takes a row index as input
#We need to give the variable where the data is
#And an index of columns that contain growth rate data (as these may change since some experiments may have more or less data
#For generelity I should also give the column that contains tindex
#And of course the tindex itself

plot.gr <- function(datamat, row.index, datarange, tindex.col, tindex)
{
	aineisto <- datamat[row.index,] #Take the row we want
	grdata <- aineisto[,datarange] #Take the columns that contain measurements
	
	#Okey here is a problem, since we will take the wrong time index. Not all of them are the same size -> need different variables -> Numbers not correct. Need to use the tindexes in a different way
	timepoints <- tindex[as.numeric(aineisto[tindex.col]),] #Take the correct time index
	timepoints <- timepoints[!is.na(timepoints)] #Remove NA's 
	
	#Calculate the growth rate, using linear regression
	model <- lm(unlist(grdata) ~ timepoints)

	plot(timepoints, grdata, type = "p", xlab = "Time (h)", ylab = "Distance (mm)", las = 1)
	abline(coef(model)) #Add the 
	
}

#Usage example, plots row 3
#plot.gr(aineisto, 3, 7:18, 4, temp.tindex)


########### Some quality control functions ###################

#####

#This functions checks if  there are observations in the growth rate data that are likely to be human errors
#Look for data points that are smaller than the previous one or larger than the next one
qc.neuro <- function(datamat, datarange)
{
	nrows <- nrow(datamat)
	
	#Loop over samples
	for(i in 1:nrows)
	{
		aineisto <- datamat[i,datarange] #Take row and those columns that contain measurements
		aineisto <- aineisto[!is.na(aineisto)] #Remove NA's
		nobs <- length(aineisto) #Number of measurements
		#print(i) #debug
		
		#Need to check some conditions first
		check1 <- all(aineisto == 0) #All observations 0
		check2 <- nobs == 0 #All observations were NA
		check3 <- nobs < 3 #Only two or less observations
		conditions <- c(check1, check2, check3)
				
		if(any(conditions) == FALSE) { #Check that samples passes tests
		#Loop over all measurements (from the second measurement to the last
			for(j in 2:nobs)
			{
				if( (aineisto[j] > aineisto[j-1]) == FALSE) { cat("Problem in sample", i, "check the data", "\n") }
				if(j < nobs) {
					if( (aineisto[j] < aineisto[j+1]) == FALSE) { cat("Problem in sample", i, "check the data", "\n") }
				}
			}
		}
		
	}
	
}

#############################################################

#####
#This function calculates linear growth rate for each sample in the data matrix
#Store growth rate, and R^2 values for troubleshooting
calc.lin.gr <- function(datamat, datarange, tindex.col, tindex, int = FALSE)
{
	#Make result matrix first
	nrows <- dim(datamat)[1]

        if(int == FALSE) {
	res.mat <- matrix(rep(0, nrows*2), ncol = 2)
	colnames(res.mat) <- c("growthrate", "R2")
        }

        if(int == TRUE) {
            res.mat <- matrix(rep(0, nrows*3), ncol = 3)
            colnames(res.mat) <- c("growthrate", "intercept", "R2")
        }
	
	#Loop over all samples
	for(i in 1:nrows)
	{
		aineisto <- datamat[i,] #Take a row
		grdata <- datamat[i,datarange] #Take the columns that contain measurements for row i
		timepoints <- tindex[as.numeric(datamat[i,tindex.col]),] #Take the correct time index
		timepoints <- timepoints[!is.na(timepoints)] #Remove NA's 
		
		#Need to perform NA check here
		if(all(is.na(grdata)) == TRUE) { res.mat[i,1] <- NA; res.mat[i,2] <- NA } else {

		#Calculate growth rate
		model <- lm(unlist(grdata) ~ timepoints)
		
		#Store the data
                if(int == FALSE) {
		res.mat[i,1] <- coef(model)[2] #Store linear growth rate
		res.mat[i,2] <- summary(model)$adj.r.squared #Store R²
		}

                if(int == TRUE) {
                    res.mat[i,1] <- coef(model)[2] #Store linear growth rate
                    res.mat[i,2] <- coef(model)[1] #Store the intercept
                    res.mat[i,3] <- summary(model)$adj.r.squared #Store R²
                }
              }
        }

	return(res.mat)
}

###########################

		
#Usage example
#testi <- calc.lin.gr(aineisto, 7:18, 4, temp.tindex)

#pdf("tarkastelua.pdf")
#plot(dim2$Level, dim2$growthrate, type = "b", col = "red", ylim = c(1,5), xlab = "Temperature", ylab = "Growth rate (mm / h)")
#points(cont$Level, cont$growthrate, type = "b", col = "black")
#points(set7$Level, set7$growthrate, type = "b", col = "blue")
#points(set2$Level, set2$growthrate, type = "b", col = "hotpink")
#points(hda2$Level, hda2$growthrate, type = "b", col = "yellow")
#points(nst2$Level, nst2$growthrate, type = "b", col = "green")
#points(dcl2$Level, dcl2$growthrate, type = "b", col = "cyan")
#dev.off()

#Calculating growth rates between each time point

#Need a function that calculates growth rate between each time point
#Function assumes that distances are a matrix (drop all other columns)
#and times and a vector
calc.each.gr <- function(distances, times) {
    #How many data point intervals there are?
    no.interval <- ncol(distances) - 1
    no.sample <- nrow(distances)
    times <- times[is.na(times) == FALSE] #Dropping time points with no data

    #Initialize the results matrix
    res.mat <- matrix(rep(0, no.interval*no.sample), ncol = no.interval)

    for(i in 1:no.interval) {
        res.mat[,i] <- (distances[,i+1] - distances[,i]) / (times[i+1] - times[i])
    }

    return(res.mat)
}




#Some functions with splines

#This function fits a spline for a single genotype and the result
plot.spline <- function(aineisto, genotype)
{
    datamat <- aineisto[aineisto$Genotype == genotype,] #Take data for a genotype
    fit <- spline(datamat$Level, datamat$growthrate, n = 50, method = "natural") #Fit the spline
    
    plot(datamat$Level, datamat$growthrate, las = 1)
    lines(fit$x, fit$y, col = "blue", lwd = 2)
}

#This function fits a spline for each genotype and calculates the optimum and other values
#Takes only one environment at a time
calc.spline <- function(aineisto)
{
    ngenot <- nlevels(aineisto$Genotype) #Number of genotypes
    genotypes <- levels(aineisto$Genotype) #genotype vector
    res.mat <- matrix(rep(0, ngenot*4), ncol = 4) #Make the results matrix

    for(i in 1:ngenot)
        {
            genodata <- aineisto[aineisto$Genotype == genotypes[i],] #Take a genotype
            sfunktio <- splinefun(genodata$Level, genodata$growthrate, method = "natural")

            #now need to calculate the optimum
            #Take min and max of the environmental variable
            min.env <- min(aineisto$Level)
            max.env <- max(aineisto$Level)

            optimum <- optimise(sfunktio, c(min.env,max.env), maximum = TRUE, deriv = 0)
            min.env.gr <- sfunktio(min.env, deriv = 0)
            max.env.gr <- sfunktio(max.env, deriv = 0)
            res.mat[i,1:4] <- c(optimum$maximum, optimum$objective, min.env.gr, max.env.gr) #Store the results
        }

    #Clean up the results
    res.mat <- data.frame(genotypes, res.mat)
    colnames(res.mat) <- c("Genotype", "opt.env", "opt.gr", "min.env.gr", "max.env.gr")
    return(res.mat)
}

### * Some wrapper functions for Bayesian analysis

#This function calculates contrasts from posterior distributions
#Takes a contrast matrix as input (e.g. mutant - control), assumes that columns are named (genotype names)
posterior.contrasts <- function(contrast.mat)
    {
        ncols <- ncol(contrast.mat)
        contrast.res <- array(0, c(ncols,4))
        for(i in 1:ncols) {
            contrast.res[i,c(1,3)] <- quantile(contrast.mat[,i], probs = c(0.025, 0.975))
            contrast.res[i,2] <- mean(contrast.mat[,i])
            contrast.res[i,4] <- sum(contrast.mat[,i] > 0) / (nrow(contrast.mat))
        }

        contrast.res <- data.frame(factor(colnames(contrast.mat)), contrast.res)
        colnames(contrast.res) <- c("Genotype", "lower", "mean", "upper", "P(dif >= 0)")
        
        return(contrast.res)
    }

#This function extracts standardised posterior values and returns them to the original scale, expects one factor
#Takes as input the mcmc results, a0 column index, number of levels of a, original mean, original sd, number of chains, names of a
#Remember to change a0 column index and number of levels of a when using the function!
extract.stand.rescale <- function(jags.results, a0.col.ind, nlevel.a, orig.mean, orig.sd, n.chains, xnames)
    {
        #Extract a values from mcmc results
        a0Sample <- NULL
        for(i in 1:n.chains) { a0Sample <- c(a0Sample, jags.results[[i]][,a0.col.ind]) } #extract values from all chains
        chainlength <- length(a0Sample) #Total length of all chains
        one.chain <- chainlength/n.chains #Length of a single chain
        aSample <- array(0, dim = c(chainlength, nlevel.a))
        for(i in 1:n.chains) { aSample[ (1+((i-1)*one.chain)): (i*one.chain), 1:nlevel.a ] <- jags.results[[i]][,1:nlevel.a] }
        #Convert to zero-centred b values
        mSample <- array(0, dim = c(chainlength, nlevel.a))
        for(i in 1:chainlength) { mSample[i,] <- a0Sample[i] + aSample[i,] }
        b0Sample <- apply(mSample, 1, mean)
        bSample <- mSample - matrix(rep(b0Sample, nlevel.a), ncol = nlevel.a)
        #Convert from standardised b values to the original scale b values:
        b0Sample <- b0Sample * orig.sd + orig.mean
        bSample <- bSample * orig.sd + orig.mean
        colnames(bSample) <- xnames

        return(list(b0Sample, bSample))
    }

#This is a more general function of posterior contrasts
HPD <- function(contrast.mat)
    {
        ncols <- ncol(contrast.mat)
        contrast.res <- array(0, c(ncols,4))
        for(i in 1:ncols) {
            contrast.res[i,c(1,3)] <- quantile(contrast.mat[,i], probs = c(0.025, 0.975))
            contrast.res[i,2] <- mean(contrast.mat[,i])
            contrast.res[i,4] <- sum(contrast.mat[,i] > 0) / (nrow(contrast.mat))
        }

        contrast.res <- data.frame(contrast.res)
        colnames(contrast.res) <- c("lower", "mean", "upper", "P(dif >= 0)")
        
        return(contrast.res)
    }


### Functions to deal with fluctuating environment
#Mainly manipulation of time strings

#This function converts the time records obtained from growth chamber logs to running time
convert.time <- function(time.col, date.col, orig.type = "american") {
    #Make the concatenated time string
    fulltime <- paste(time.col, " ", date.col, sep = "")
    #Then convert it into a time format
    #If type = "american" this is month/day/year (weird format)
    if(orig.type == "american") {
    timeconv <- strptime(fulltime, format = "%H:%M:%S %m/%d/%Y")
    }

    if(orig.type == "euro") {
    timeconv <- strptime(fulltime, format = "%H.%M.%S %d.%m.%Y")
    }

    #Time in minutes starting from 0
    timerun <- as.numeric(timeconv - timeconv[1])/60 # divided by 60 as returns seconds

    return(timerun)
}


#This function calculates genotype by environment interactions for all mutants
#To test for differences in reaction norm shape againts control
#Note 20.07.2017: I've modified this function to be more general. now it needs column numbers for growthrate column, genotype column, and environment column for input
#This change may affect some older scripts, update to new version as needed
pairwise.ANOVA.GxE <- function(datamat, mutants, round = TRUE, control.name = "cont4200", genotype.col.no, environ.col.no, phenotype.col.no)
    {
        pairwise.res <- matrix(rep(0, length(mutants)*7), ncol = 7)
        colnames(pairwise.res) <- c("Genotype", "F-value Genotype", "p-value Genotype", "p-cor Genotype", "F-value GxE", "p-value GxE", "p-cor GxE")
        #pairwise.res[,1] <- mutants
        #Take only control and mutant
        for(i in 1:length(mutants)) {
            datamat2 <- datamat[datamat[,genotype.col.no] == control.name | datamat[,genotype.col.no] == mutants[i],]
            datamat2[,genotype.col.no] <- factor(datamat2[,genotype.col.no]) #Clean genotype factor
            result <- anova(lm(datamat2[,phenotype.col.no] ~ datamat2[,genotype.col.no] + datamat2[,environ.col.no] + datamat2[,genotype.col.no]:datamat2[,environ.col.no]))
            pairwise.res[i,5] <- result[3,4] #Store F-value for GxE
            pairwise.res[i,6] <- result[3,5] #Store p-value for GxE
            pairwise.res[i,2] <- result[1,4] #Store F-value for Genotypic effect
            pairwise.res[i,3] <- result[1,5] #Store p-value for Genotypic effect
        }

        #Need to apply multiple testing correction, otherwise reviewers are sad...  :´(
        #Bonferroni-Holm
        pairwise.res[,7] <- adj.p(pairwise.res[,6], method = "Holm")
        pairwise.res[,4] <- adj.p(pairwise.res[,3], method = "Holm")
        pairwise.res <- data.frame(pairwise.res)
        pairwise.res[,1] <- mutants

        if(round == TRUE) { #Round values if needed for tables etc.
            pairwise.res[,2] <- sprintf("%.3f", pairwise.res[,2]) #Round F-values to 3 digits
            pairwise.res[,5] <- sprintf("%.3f", pairwise.res[,5])

            pairwise.res[,3] <- sprintf("%.3G", pairwise.res[,3]) #Round p-values to 3 digits but keep exp. notation
            pairwise.res[,4] <- sprintf("%.3G", pairwise.res[,4])
            pairwise.res[,6] <- sprintf("%.3G", pairwise.res[,6])
            pairwise.res[,7] <- sprintf("%.3G", pairwise.res[,7])
        }

        #Performing some formatting for the 

        return(pairwise.res)
    }

        
#This function calculates adjusted p-values given a vector of p.values
adj.p <- function(pvalues, method = c("Bonferroni", "Holm"))
{
	n <- length(pvalues) #number of observations
		
	if(method == "Bonferroni") {
		adjusted <- pvalues*n #Bonferroni method Rp(r)
		for(i in 1:n) { if(adjusted[i] > 1) adjusted[i] <- 1 } #Setting p-values larger than 1 to 1
	}
	
	if(method == "Holm") {
		o <- order(pvalues) #order of pvalues
		pvalues.sorted <- sort(pvalues) #Sorting
		holm.p <- rep(0,n)
		for(i in 1:n) { 
			holm.p[i] <- (n - i + 1)*pvalues.sorted[i]  #Method of Holm
			if(i >= 2) {
				if(holm.p[i] < holm.p[i-1]) holm.p[i] <- holm.p[i-1] } #Correct sequential Bonferroni values to Holm values, i. e. no greater p-value may appear earlier in sequence
		}
		for(j in 1:n) { if(holm.p[j] > 1) holm.p[j] <- 1 }
		sortmat <- cbind(holm.p,o) #Back sorting adjusted p-values
		o2 <- order(sortmat[,2])
		adjusted <- sortmat[order(sortmat[,2]),1]
	}
	
	return(adjusted)
}

#This function converts spectophometric absorbance measured at 600 nm into conidial numbers
#Using a standard curve
abs2conidia <- function(data) { exp((data+13.8268)/0.810368) - 2.70826*10^7 }

#This function plots the reaction norm data
plot.geno <- function(aineisto, geno) {
    koe <- filter(aineisto, Genotype == geno)
    ggplot(koe, aes(x = Temp, y = growthrate)) +
       geom_point()
}

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
#summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
#                      conf.interval=.95, .drop=TRUE) {
#    require(plyr)
#
#    # New version of length which can handle NA's: if na.rm==T, don't count them
#    length2 <- function (x, na.rm=FALSE) {
#        if (na.rm) sum(!is.na(x))
#        else       length(x)
#    }
#
#    # This does the summary. For each group's data frame, return a vector with
#    # N, mean, and sd
#    datac <- ddply(data, groupvars, .drop=.drop,
#      .fun = function(xx, col) {
#        c(N    = length2(xx[[col]], na.rm=na.rm),
#          mean = mean   (xx[[col]], na.rm=na.rm),
#          sd   = sd     (xx[[col]], na.rm=na.rm)
#        )
#      },
#      measurevar
#    )
#
#    # Rename the "mean" column    
#    datac <- rename(datac, c("mean" = measurevar))
#
#    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
#
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
#    datac$ci <- datac$se * ciMult
#
#   return(datac)
#}


##This function standardizes the data. Input is the data to be standardized
data.std <- function(x) {
    ( x - mean(x) )/sd(x)
}


##This function is for sampling G and E matrices from a multivariate posterior
##Used in ms: "Quantitative genetics of temperature tolerance curves of Neurospora crassa".
sample.matrices.from.posterior <- function(mvpost) {

    rowind <- sample(1:nrow(mvpost), size = 1) #Sample row index
    
##Sample the G-matrix
G <- matrix(rep(0, 6*6), ncol = 6)
diag(G) <- mvpost[rowind,7:12]^2 #Sample genetic variance to diagonal
G[1,2] <- mvpost[rowind,7]*mvpost[rowind,8]*mvpost[rowind,13]
G[1,3] <- mvpost[rowind,7]*mvpost[rowind,9]*mvpost[rowind,14]
G[1,4] <- mvpost[rowind,7]*mvpost[rowind,10]*mvpost[rowind,16]
G[1,5] <- mvpost[rowind,7]*mvpost[rowind,11]*mvpost[rowind,19]
G[1,6] <- mvpost[rowind,7]*mvpost[rowind,12]*mvpost[rowind,23]
G[2,3] <- mvpost[rowind,8]*mvpost[rowind,9]*mvpost[rowind,15]
G[2,4] <- mvpost[rowind,8]*mvpost[rowind,10]*mvpost[rowind,17]
G[2,5] <- mvpost[rowind,8]*mvpost[rowind,11]*mvpost[rowind,20]
G[2,6] <- mvpost[rowind,8]*mvpost[rowind,12]*mvpost[rowind,24]
G[3,4] <- mvpost[rowind,9]*mvpost[rowind,10]*mvpost[rowind,18]
G[3,5] <- mvpost[rowind,9]*mvpost[rowind,11]*mvpost[rowind,21]
G[3,6] <- mvpost[rowind,9]*mvpost[rowind,12]*mvpost[rowind,25]
G[4,5] <- mvpost[rowind,10]*mvpost[rowind,11]*mvpost[rowind,22]
G[4,6] <- mvpost[rowind,10]*mvpost[rowind,12]*mvpost[rowind,26]
G[5,6] <- mvpost[rowind,11]*mvpost[rowind,12]*mvpost[rowind,27]

#Take upper triangle and replace the lower triangle
G[2:6,1] <- G[1,2:6]
G[3:6,2] <- G[2,3:6]
G[4:6,3] <- G[3,4:6]
G[5:6,4] <- G[4,5:6]
G[6,5] <- G[5,6]

##The environmental variance-covarianve matrix E
E <- matrix(rep(0, 6*6), ncol = 6)
diag(E) <- mvpost[rowind,43:48]^2 #Sample environmental variance to diagonal
E[1,2] <- mvpost[rowind,43]*mvpost[rowind,44]*mvpost[rowind,28]
E[1,3] <- mvpost[rowind,43]*mvpost[rowind,45]*mvpost[rowind,29]
E[1,4] <- mvpost[rowind,43]*mvpost[rowind,46]*mvpost[rowind,31]
E[1,5] <- mvpost[rowind,43]*mvpost[rowind,47]*mvpost[rowind,34]
E[1,6] <- mvpost[rowind,43]*mvpost[rowind,48]*mvpost[rowind,38]
E[2,3] <- mvpost[rowind,44]*mvpost[rowind,45]*mvpost[rowind,30]
E[2,4] <- mvpost[rowind,44]*mvpost[rowind,46]*mvpost[rowind,32]
E[2,5] <- mvpost[rowind,44]*mvpost[rowind,47]*mvpost[rowind,35]
E[2,6] <- mvpost[rowind,44]*mvpost[rowind,48]*mvpost[rowind,39]
E[3,4] <- mvpost[rowind,45]*mvpost[rowind,46]*mvpost[rowind,33]
E[3,5] <- mvpost[rowind,45]*mvpost[rowind,47]*mvpost[rowind,36]
E[3,6] <- mvpost[rowind,45]*mvpost[rowind,48]*mvpost[rowind,40]
E[4,5] <- mvpost[rowind,46]*mvpost[rowind,47]*mvpost[rowind,37]
E[4,6] <- mvpost[rowind,46]*mvpost[rowind,48]*mvpost[rowind,41]
E[5,6] <- mvpost[rowind,47]*mvpost[rowind,48]*mvpost[rowind,42]

#Take upper triangle and replace the lower triangle
E[2:6,1] <- E[1,2:6]
E[3:6,2] <- E[2,3:6]
E[4:6,3] <- E[3,4:6]
E[5:6,4] <- E[4,5:6]
E[6,5] <- E[5,6]

#I also want the mean phenotypes for possible standardizations
z <- matrix(mvpost[1,1:6], ncol = 1)

return(list(G = G, E = E, z = z))

}

abssum <- function(x) {sum(abs(x))}

##This function simulates multivariate selection response
##Uncertainty in G-matrix is taken into account by sampling new matrices from the posterior distribution
quant.gen.sim <- function(mvpost, S, M, output = "phenotypes", method = "sel.dif") {

replicates <- 1000
res.mat <- matrix(rep(0,replicates*6), ncol = 6)
gradients <- matrix(rep(0,replicates*6), ncol = 6)

##Calculating response to selection using selection differentials
if(method == "sel.dif") {
    ##Repeat many simulations
    for(i in 1:replicates) {  ##To do
        sampledmats <- sample.matrices.from.posterior(mvpost)
        G <- sampledmats$G
        E <- sampledmats$E
        z <- sampledmats$z
        #
        P <- G + E

        R <- G%*%solve(P)%*%S #Calculate response to selection
        Beta <- solve(P)%*%S #Calculate realized selection gradient
        #
        res.mat[i,] <- R #Store the predicted responses
        gradients[i,] <- Beta #Store the realized selection gradient
        #
    }
}

if(method == "Beta") {
    ##Repeat many simulations
    for(i in 1:replicates) {  ##To do
        sampledmats <- sample.matrices.from.posterior(mvpost)
        G <- sampledmats$G
        E <- sampledmats$E
        z <- sampledmats$z
        #
        P <- G + E

        R <- G%*%S #Calculate response to selection as if S is Beta
        res.mat[i,] <- R #Store the predicted responses
    }
}

if(output == "phenotypes") {
gen <- c(1,3,5) #Number of generations of selection
pred.pheno <- NULL
for(j in 1:length(gen)) {
    #
    new.pheno <- res.mat
    for(i in 1:replicates) { new.pheno[i,] <- new.pheno[i,]*gen[j] + M }
                                        #Mnew <- M + R*gen #Predicted population mean after selection
    predictions <- apply(new.pheno, 2, quantile, probs = c(0.025, 0.165, 0.5, 0.835, 0.975))
    #
    pred.pheno <- rbind(pred.pheno, data.frame(Temp = Temps, gr = predictions[3,], grmin = predictions[1,], grmax = predictions[5,], grmidl = predictions[2,], grmidh = predictions[4,], generations = rep(gen[j],6)))
}

return(pred.pheno)

}

if(output == "response") {

    responses <- apply(res.mat, 1, abssum)
    responses.quant <- quantile(responses, probs = c(0.025, 0.5, 0.975))
    selresponses <- data.frame(Rmin = responses.quant[1], Rmean = responses.quant[2], Rmax = responses.quant[3])

    return(selresponses) }

if(output == "gradient") {

    gradient.quant <- apply(gradients, 2, quantile, probs = c(0.025, 0.5, 0.975))
    gradient.quant <- t(gradient.quant)
    gradient.res <- data.frame(gradient.quant)
    gradient.res$temperature <- c(20, 25, 30, 35, 37.5, 40)
    colnames(gradient.res) <- c("Bmin", "Bmed", "Bhigh", "temperature")

    return(gradient.res) }

}

##This function calculates the angle between response to selection and selection gradient
##Does many replicate simulations
quant.gen.angle <- function(mvpost, S) {
    #
    replicates <- 1000
    res.mat <- matrix(rep(0, replicates*1), ncol = 1)
    #
    ##Repeat many simulations
    for(i in 1:replicates) {
    #
        ##Examining selection gradients
        sampledmats <- sample.matrices.from.posterior(mvpost)
        G <- sampledmats$G
        E <- sampledmats$E
        #                                #
        P <- G + E
        R <- G%*%solve(P)%*%S #Calculate response to selection
        Beta <- solve(P)%*%S #Calculate selection gradient
        #
        ##Calculate angle between R and Beta
        costheta <- (t(Beta)%*%R) / (sqrt(t(R)%*%R)*sqrt(t(Beta)%*%Beta))
        theta <- acos(costheta)*(180/pi)
        #
        res.mat[i,] <- theta #Store the angle between response to selection and selection gradient
    }
#    #
#
return(quantile(res.mat, probs = c(0.025, 0.5, 0.975)))
#
}


##Additional functions to calculate evolvabilities, conditional evolvabilities, and autonomies according to Hansen & Hoyle 2008

#Calculate evolvability given G-matrix and selection gradient
evolvability <- function(G, Beta) {
   evol <- (t(Beta)%*%G%*%Beta) / norm_vec(Beta)^2
   return(evol)
}

#Calculate respondability given G-matrix and selection gradient
respondability <- function(G, Beta) {
    res <- sqrt( t(Beta)%*%G%*%G%*%Beta ) / norm_vec(Beta)
    return(res)
}

#Calculate conditional evolvability given G-matrix and selection gradient
cond.evol <- function(G, Beta) {
    unit.Beta <- Beta / norm_vec(Beta) #unit vector in the direction of Beta
    condevol <- (t(unit.Beta)%*%solve(G)%*%unit.Beta)^(-1)
    return(condevol)
}

#Function to calculate the norm (length) of a vector
norm_vec <- function(x) sqrt(sum(x^2))

#Conditional evolvabilities per trait
con.trait <- function(G) {
    inv.G <- solve(G)
    return(1/(diag(inv.G)))
}

#Autonomies per trait
autom.trait <- function(G) {
    inv.G <- solve(G)
    a <- (diag(inv.G) * diag(G))^(-1)
    return(a)
}

#Calculating expectations for evolvability statistics over a uniform distribution of selection gradients in the entire phenotypic space

#Calculating variance of eigenvalues (includes sample size correction)
var.eig <- function(eigval) {
    k <- length(eigval) #Number of eigenvalues
    eigvar <- sum(eigval^2)/k - (sum(eigval)/k)^2
    return(eigvar)
}

Har <- function(x) { 1/mean(1/x) } #Harmonic mean
Ie <- function(x) { var.eig(x) / (mean(x))^2 } #Variance for eigenvalues

#Expected evolvability
evol.all <- function(G) {
    eigval <- eigen(G)$values
    return(mean(eigval))
}

#Conditional evolvability
conde.all <- function(G) {
    eigval <- eigen(G)$values
    k <- length(eigval) #Number of eigenvalues
    meanc <- Har(eigval)*(1 + (2*Ie(1/eigval))/(k+2))
    return(meanc)
}

res.all <- function(G) {
    eigval <- eigen(G)$values
    k <- length(eigval)
    meanr <- sqrt(mean(eigval^2))*(1 - (Ie(eigval^2))/( 4*(k+2) ) )
    return(meanr)
}

autom.all <- function(G) {
    eigval <- eigen(G)$values
    k <- length(eigval)
    meana <- Har(eigval)/mean(eigval) * (1 + 2*( (Ie(eigval) + Ie(1/eigval) - 1 + Har(eigval)/mean(eigval) + 2*Ie(eigval)*Ie(1/eigval)/(k+2) ) / (k + 2) ))
    return(meana)
}

##This function is the same as sample.matrices.from.posterior but instead using all posterior samples
matrices.from.posterior <- function(mvpost, rowind) {

    #rowind <- sample(1:nrow(mvpost), size = 1) #Sample row index
    
##Sample the G-matrix
G <- matrix(rep(0, 6*6), ncol = 6)
diag(G) <- mvpost[rowind,7:12]^2 #Sample genetic variance to diagonal
G[1,2] <- mvpost[rowind,7]*mvpost[rowind,8]*mvpost[rowind,13]
G[1,3] <- mvpost[rowind,7]*mvpost[rowind,9]*mvpost[rowind,14]
G[1,4] <- mvpost[rowind,7]*mvpost[rowind,10]*mvpost[rowind,16]
G[1,5] <- mvpost[rowind,7]*mvpost[rowind,11]*mvpost[rowind,19]
G[1,6] <- mvpost[rowind,7]*mvpost[rowind,12]*mvpost[rowind,23]
G[2,3] <- mvpost[rowind,8]*mvpost[rowind,9]*mvpost[rowind,15]
G[2,4] <- mvpost[rowind,8]*mvpost[rowind,10]*mvpost[rowind,17]
G[2,5] <- mvpost[rowind,8]*mvpost[rowind,11]*mvpost[rowind,20]
G[2,6] <- mvpost[rowind,8]*mvpost[rowind,12]*mvpost[rowind,24]
G[3,4] <- mvpost[rowind,9]*mvpost[rowind,10]*mvpost[rowind,18]
G[3,5] <- mvpost[rowind,9]*mvpost[rowind,11]*mvpost[rowind,21]
G[3,6] <- mvpost[rowind,9]*mvpost[rowind,12]*mvpost[rowind,25]
G[4,5] <- mvpost[rowind,10]*mvpost[rowind,11]*mvpost[rowind,22]
G[4,6] <- mvpost[rowind,10]*mvpost[rowind,12]*mvpost[rowind,26]
G[5,6] <- mvpost[rowind,11]*mvpost[rowind,12]*mvpost[rowind,27]

#Take upper triangle and replace the lower triangle
G[2:6,1] <- G[1,2:6]
G[3:6,2] <- G[2,3:6]
G[4:6,3] <- G[3,4:6]
G[5:6,4] <- G[4,5:6]
G[6,5] <- G[5,6]

##The environmental variance-covariance matrix E
E <- matrix(rep(0, 6*6), ncol = 6)
diag(E) <- mvpost[rowind,43:48]^2 #Sample environmental variance to diagonal
E[1,2] <- mvpost[rowind,43]*mvpost[rowind,44]*mvpost[rowind,28]
E[1,3] <- mvpost[rowind,43]*mvpost[rowind,45]*mvpost[rowind,29]
E[1,4] <- mvpost[rowind,43]*mvpost[rowind,46]*mvpost[rowind,31]
E[1,5] <- mvpost[rowind,43]*mvpost[rowind,47]*mvpost[rowind,34]
E[1,6] <- mvpost[rowind,43]*mvpost[rowind,48]*mvpost[rowind,38]
E[2,3] <- mvpost[rowind,44]*mvpost[rowind,45]*mvpost[rowind,30]
E[2,4] <- mvpost[rowind,44]*mvpost[rowind,46]*mvpost[rowind,32]
E[2,5] <- mvpost[rowind,44]*mvpost[rowind,47]*mvpost[rowind,35]
E[2,6] <- mvpost[rowind,44]*mvpost[rowind,48]*mvpost[rowind,39]
E[3,4] <- mvpost[rowind,45]*mvpost[rowind,46]*mvpost[rowind,33]
E[3,5] <- mvpost[rowind,45]*mvpost[rowind,47]*mvpost[rowind,36]
E[3,6] <- mvpost[rowind,45]*mvpost[rowind,48]*mvpost[rowind,40]
E[4,5] <- mvpost[rowind,46]*mvpost[rowind,47]*mvpost[rowind,37]
E[4,6] <- mvpost[rowind,46]*mvpost[rowind,48]*mvpost[rowind,41]
E[5,6] <- mvpost[rowind,47]*mvpost[rowind,48]*mvpost[rowind,42]

#Take upper triangle and replace the lower triangle
E[2:6,1] <- E[1,2:6]
E[3:6,2] <- E[2,3:6]
E[4:6,3] <- E[3,4:6]
E[5:6,4] <- E[4,5:6]
E[6,5] <- E[5,6]

#I also want the mean phenotypes for possible standardizations
z <- matrix(mvpost[1,1:6], ncol = 1)

return(list(G = G, E = E, z = z))

}

##Hansen's evolvability calculations
calc.hansen <- function(mvpost, standardization = "mean") {
    iterations <- nrow(mvpost) #Number of iterations

    ##Setup results matrices
    #Single traits
    trait.results <- matrix(rep(0, iterations*6*3), ncol = 6*3)
    colnames(trait.results) <- c(paste(rep("evol.gr", 6), c(20, 25, 30, 35, 375, 40), sep = ""), paste(rep("condevol.gr", 6), c(20, 25, 30, 35, 375, 40), sep = ""), paste(rep("autom.gr", 6), c(20, 25, 30, 35, 375, 40), sep = ""))

    #Expectation over entire phenotypic space
    exp.results <- matrix(rep(0, iterations*4), ncol = 4)
    colnames(exp.results) <- c("mean.e", "mean.r", "mean.c", "mean.a")

    #Loop over all posterior samples
    for(i in 1:iterations) {

        matrices <- matrices.from.posterior(mvpost, rowind = i) #Get matrices
        G <- matrices$G
        E <- matrices$E
        z <- matrices$z

        #Check for standardization
        if(standardization == "mean") {
            G <- G / (z %*% t(z)) #Element-wise division of G
        }

        ##Do calculations
        #Per trait measures
        trait.results[i,1:6] <- diag(G) #Evolvabilities per trait are just G diagonals
        trait.results[i,7:12] <- con.trait(G) #Conditional evolvabilities per trait
        trait.results[i,13:18] <- autom.trait(G) #Autonomies per trait

        #Expectations
        exp.results[i,1] <- evol.all(G) #Evolvability
        exp.results[i,2] <- res.all(G) #Respondability
        exp.results[i,3] <- conde.all(G) #Conditional evolvability
        exp.results[i,4] <- autom.all(G) #Autonomy

    } #Done looping over all iterations

    return(list(trait.results = trait.results, exp.results = exp.results))
}

##Hansen's evolvability calculations, for a given selection gradient
#Takes selection differentials as input and calculates Beta from this
calc.hansen.beta <- function(mvpost, S, standardization = "mean", method = "sel.dif") {
    #iterations <- nrow(mvpost) #Number of iterations
    replicates <- 1000 #Sample from the posterior

    #Setup results matrices
    resmat <- matrix(rep(0, replicates*8), ncol = 8)
    colnames(resmat) <- c("Beta.e", "Beta.c", "Beta.r", "Beta.a", "mean.e", "mean.c", "mean.r", "mean.a")

    if(standardization == "mean") {
        for(i in 1:replicates) {  ##To do
            sampledmats <- sample.matrices.from.posterior(mvpost)
            G <- sampledmats$G
            E <- sampledmats$E
            z <- sampledmats$z
                                        #
            P <- G + E

            if(method == "sel.dif") {
            Beta <- solve(P)%*%S #Calculate realized selection gradient
            }
            if(method == "Beta") {
            Beta <- S #Use selection gradient directly
            }

            #Mean standardization
            G <- G / (z %*% t(z)) #Element-wise division of G
            Beta <- z * Beta #Element-wise multiplication of Beta

            #Evolvability calculations
            resmat[i,1] <- evolvability(G, Beta) #Evolvability
            resmat[i,2] <- cond.evol(G, Beta) #Conditional evolvability
            resmat[i,3] <- respondability(G, Beta) #Respondability
            resmat[i,4] <- resmat[i,2]/resmat[i,1] #Autonomy
            resmat[i,5] <- evol.all(G) #Expected evolvability
            resmat[i,6] <- conde.all(G) #Expected conditional evolvability
            resmat[i,7] <- res.all(G) #Expected respondability
            resmat[i,8] <- autom.all(G) #Expected autonomy
        }
    }

    return(resmat)  
}
