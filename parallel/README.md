This is for parallel computing version of original pipeline.

## Step 1 - produce point estimations for all parameters

Input: GRM file
Output: point estimates for the whole dataset

## Step 2 - individual-level jackknife

Input: GRM file
Output: variance of ind-jackknife for all parameters


## Step 3 - snp-level jackknife

Input: GRM file, fam file path, dosage file path, map file path
Output: variance of snp-jackknife for all parameters
