# GV-paritioning-tmp
partitioning heritability and genetic covariance by direct and indirect paths


### Input
1. Dosage files for sibling genotypes and parental genotypes
The format follows PLINK dosage files, including columns of rs number, allele 1, allele 2, and dosages for each individual. Do not need to add column names. In our implementation we are using single dosage for each individual. Sibling dosages and parental dosages should be stacked into one dosage file. Chromosomes seperated dosage files are preferred.
2. map files
map files are used to perform SNP-level jackknife analysis. Format follows PLINK .map files. Every block corresponds to one map file.
3. Overall .fam file
An overall fam file is needed to include all participants (including siblings and parents). 

### Output
Estimation of all parameters and their variance estimations.

## Step1


### Usage



## Step 2

### Usage

