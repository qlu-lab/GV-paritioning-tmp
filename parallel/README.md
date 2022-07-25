This is for parallel computing version of original pipeline.

## Step 1 - produce point estimations for all parameters

Input: GRM file

Output: point estimates for the whole dataset

Usage: 
```{r}
Rscript 1_pointest.r [path to phenotype 1] [path to phenotype 2] [grm_path] [grm_id_path] [dout] 
```


## Step 2 - individual-level jackknife

Input: GRM file

Output: variance of ind-jackknife for all parameters

Usage:
```{r}
Rscript 2_indjack.r [path to phenotype 1] [path to phenotype 2] [total number of blocks] [ith block] [grm_path] [grm_id_path] [dout]
```

## Step 3 - snp-level jackknife

Input: GRM file, fam file path, dosage file path, map file path

Output: variance of snp-jackknife for all parameters

Usage:

```{r}
Rscript 3_snpjack.r [pheno_y1] [pheno_y2] [CHR] [snpblock] [k_snp] 
```
where 

pheno_y1: file path for the first phenotype

phenotype_y2: file path for the second phenotype

CHR: chromosome number {1...22}

snpblock: ith block for chromosome CHR

k_snp: jth block for all blocks
