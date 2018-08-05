#!/bin/bash -ve


#$ -l h_rt=7:59:00

#$ -V

#$ -j y

#$ -l mem=4G 

#$ -l rmem=4G
#####################################################################
####### Genomewide
#####################################################################

################# relatedness matrix - HEIDI ########################
# size
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr genomewide --out results/size/size --maxmem 10 

# shape
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr genomewide --out results/shape/shape --maxmem 10

#####################################################################
##### Partition boundaries
#####################################################################

######## region1 : are in the inversions

##### Size
# LG 6
#/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr 6 --out results/size/size --maxmem 10 --from lg6_ord1_scaf7009:4390 --to lg6_ord21_scaf277724:791 --id 1

# LG 14
#/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr 14 --out results/size/size --maxmem 10 --from lg14_ord1_scaf190620:7774 --to lg14_ord16_scaf7422:3513 --id 1

# LG 17
#/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr 17 --out results/size/size --maxmem 10 --from lg17_ord57_scaf74985:227 --to lg17_ord71_scaf1755:13977 --id 1

########## shape
# LG 6
#/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr 6 --out results/shape/shape --maxmem 10 --from lg6_ord1_scaf7009:4390 --to lg6_ord21_scaf277724:791 --id 1

# LG 14
#/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr 14 --out results/shape/shape --maxmem 10 --from lg14_ord1_scaf190620:7774 --to lg14_ord16_scaf7422:3513 --id 1

# LG 17
#/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr 17 --out results/shape/shape --maxmem 10 --from lg17_ord57_scaf74985:227 --to lg17_ord71_scaf1755:13977 --id 1

######## region2 : are rest of the chromosomes

##### Size
# LG 6
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr 6 --out results/size/size --maxmem 10 --from lg6_ord22_scaf77331:5393 --to lg6_ord74_scaf53922:28164 --id 2

# LG 14
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr 14 --out results/size/size --maxmem 10 --from lg14_ord17_scaf82138:18603 --to lg14_ord43_scaf285120:343 --id 2

# LG 17
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/size --covar pheno_file/covar/covar_size --chr 17 --out results/size/size --maxmem 10 --from lg17_ord1_scaf67673:25754 --to lg17_ord56_scaf48522:5218 --id 2

########## shape
# LG 6
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr 6 --out results/shape/shape --maxmem 10 --from lg6_ord22_scaf77331:5393 --to lg6_ord74_scaf53922:28164 --id 2

# LG 14
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr 14 --out results/shape/shape --maxmem 10 --from lg14_ord17_scaf82138:18603 --to lg14_ord43_scaf285120:343 --id 2

# LG 17
/data/bop14pc/heidi/heidi --bfile geno_file/ADM.impute --pheno pheno_file/pheno/shape --covar pheno_file/covar/covar_shape --chr 17 --out results/shape/shape --maxmem 10 --from lg17_ord1_scaf67673:25754 --to lg17_ord56_scaf48522:5218 --id 2