# GV-paritioning
partitioning heritability and genetic covariance by direct and indirect paths

### Input
1. fn_y1
2. fn_y2
3. chr_all_grm_zs.rel
4. chr_all_grm_zs.rel.id
5. fam_path
6. map_path
7. dosage_path
8. dout
Path of output files
9. snpjack_map
10. line1_rmed_chr
11. plink_path
Path of the plink software(not the folder)
12. K_ind
13. K_snp
14. N_grm
15. M_grm


1. Dosage files for sibling genotypes and parental genotypes
The format follows PLINK dosage files, including columns of rs number, allele 1, allele 2, and dosages for each individual. Do not need to add column names. In our implementation we are using single dosage for each individual. Sibling dosages and parental dosages should be stacked into one dosage file. Chromosomes seperated dosage files are preferred.

2. map files
map files are used to perform SNP-level jackknife analysis. Format follows PLINK .map files. Every block corresponds to one map file. The set of map files used in our paper is available in our repository.

3. Overall .fam file
An overall fam file is needed to include all participants (including siblings and parents).

### Output
Estimation of all parameters and their variance estimations.

## Step1
Generate GRM matrices for the whole dataset.
This can be done using PLINK. If there are enough RAM and storage space in your machine, calculating GRM for each chromosome will save running time in step2.

## Step 2
Run GV-partitioning model.

### Usage
Rscript ~./pipeline_GCsibling.r [fn_y1] [fn_y2] [chr_all_grm_zs.rel] [chr_all_grm_zs.rel.id] [fam_path] [map_path] [dosage_path] [dout] [snpjack_map] [line1_rmed_chr] [plink_path] [K_ind] [K_snp] [N_grm] [M_grm]

If parallel system such as condor is available, we highly recommend to run jobs in parallel. This will substatially save time. A parallel version of our codes is in 'parallel' folder in this repository.

## Citations
If anyone utilized this pipeline in the paper, please cite
[our paper]
