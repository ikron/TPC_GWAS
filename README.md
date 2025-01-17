# Genome-wide association for thermal performance curves in _Neurospora crassa_

This repository describes phenotypic and data and analysis scripts related to the manuscript: [Genome-wide association for loci influencing thermal performance curves in Neurospora crassa](https://www.biorxiv.org/content/10.1101/2024.04.29.591604v1)

Note that there is another repository for the nested association mapping population itself: [Neurospora_NAM_population](https://github.com/ikron/Neurospora_NAM_population) and the genotype file is available at [zenodo](https://zenodo.org/records/11120317)

## Phenotypic data

The phenotypic data used in this manuscript is from  [Moghadam et al. 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/evo.14016) and consist of growth rate measurements of different strains of the filamentous fungus _Neurospora crassa_. This data is combined with SNP data for all of these strains to perform association mapping. The data is also available from [Dryad](https://doi.org/10.5061/dryad.pk0p2ngk9). However, the most recent version of the data is folder /data/

In folder /data/ the data files are:
1. phenotypes.csv
   - A csv file that contains the phenotypic measurements for each growth rate assay. The columns are: Block = experimental block, Temp = temperature, Family = To which family in the NAM population a strain belongs to (natural strains have 'O'), Genotype = name of the strain, tindex = index for time points, GCham = growth chamber, Replicate, replicate ID, t1 to 8 = values of many mm the mycelium has grown (the first time point should have a value of zero), Genot = name of the strain that corresponds to name in the genotype file. This data is used to calculate growth rates can generate the 'genomeans.csv' file.
2. genomeans.csv
   - A csv file that contains mean growth rates for all genotypes. Used as input data in GWAS. The columns are: Genot = name of the strain that corresponds to name in the genotype file, Temp = temperature, Family = To which family in the NAM population a strain belongs to (natural strains have 'O'), meanger = mean growth rate of the strain
3. splineresults.csv
   - A csv file that contains TPC parameters estimated from spline fits for each genotype. Used as input data in GWAS. The columns are: umax = maximum growth rate over the whole TPC estimated from spline fit, Topt = optimum temperature for the strain or the temperatire where maximum growth rate occurs, Genotype = name of the strain that corresponds to name in the genotype file
4. PCAs_umax.csv
   - A csv file that contains covariate information each strain. Used as input data in GWAS. The columns are: taxa = name of the strain that corresponds to name in the genotype file, PC1 to PC4 = value for each strain on the principle component axes that have been calculated from genotype data for population structure correction, umax = maximum growth rate estimated from spline fits as above
5. GAPIT.Genotype.PCA.csv
   - A csv file that contains covariate information each strain. Used as input data in GWAS. The columns are: taxa = name of the strain that corresponds to name in the genotype file, PC1 to PC4 = value for each strain on the principle component axes that have been calculated from genotype data for population structure correction. This is file is identical to PCAs_umax.csv but the umax columns is missing.
6. GAPIT.Genotype.PCA_eigenvalues.csv
   - A csv file that contains eigenvalues for component from the PCA performed on genotype data (this data is used to calculate % variance explained for each component)

## Analysis scripts

Scripts used for data analysis

In folder /scripts/ the script files are:
1. Neuro_functions.R
   - File that contains R scripts that are used in calculating growth rates, the are called from tempasso.R
2. association_scripts.R
   - File that contains R scripts that are used performing GWAS, these are called from tempasso.R
3. tempasso.R
   - File that contains the R commands used for data analysis. Calculating growth rates, doing summary statistics, and analyzing the GWAS results.

Some of the analyses such as LD calculations and GWAS itself have been performed on a national computer cluster the *.sh are script files that control the execuation of the R scripts on the cluster. In order to make the R scripts work in your environment, you would have to change input and output folders etc., but the commands can be looked at anyway.

4. LD_calc.R and LD_calc.sh
   - R script to sample 10 000 SNPs from each chromosome and calculates LD between them (chr by chr)
5. LD_calc_fine.R and LD_calc_fine.sh
   - R script to sample sets of SNPs evenly along the chromosomes and calculate LD between the sampled SNPs in each set
6. The script files to perform GWAS are of the form: asso_gr20_new.R where 'gr20' is the environment (in this case growth rate at 20 degrees), and cluster control file is the form Rjob_gwas_gr20_new.sh. The 'gr20_umax' and 'gr40_umax' are for performing GWAS when using maximum growth rate as covariate

## Intermediate .RData files

In folder /RData/ there are intermediate files used in some of the analyses. Such as the SNPs sampled for LD calculations, significant SNPs etc. So that these can be viewed. See the R script tempasso.R where these are generated.

