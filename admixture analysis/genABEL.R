args = commandArgs(trailingOnly=TRUE)
library(GenABEL)

#### load and explore data
dat <- load.gwaa.data(phenofile=args[1], genofile=args[2], force=TRUE, makemap=FALSE, sort=TRUE, id = "id")

#### QC
# QC -#no control for hardy weinberg here (option p.value=0), control for maf of given value, defaults for SNP and individual call rates 
qc1 <- check.marker(dat, p.level=0, maf = 0.01) 
data1 <- dat[qc1$idok, qc1$snpok] # data set containing only individuals and snps passing qc1

# prep data for GWAS: calculating kinship matrix, replacing diagonal for egscore funtion purposes, see ?egscore 
data1.gkin <- ibs(data1[, data1\@gtdata\@chromosome != "X"], weight="freq") 
diag(data1.gkin) <- hom(data1[,autosomal(data1)])$Var

##########################################################
## single_SNP GWA with correction for population structure#
## using the method of Price et al. (using PCs of the genomic kinship matrix) here implemented in the egscore function (you just have to provide the kinship matrix)
##############################################################################

#GWA
data1.pca <- egscore(args[3], data1, kin=data1.gkin, propPs=0.95, clambda=FALSE, naxes=4 , times = 1000)
#data1.pca <- egscore(args[3]~size, data=data1, kin=data1.gkin, propPs=0.95, clambda=FALSE, naxes=4 , times = 1000) # trait = shape
#data1.pca <- egscore(args[3]~sex, data=data1, kin=data1.gkin, propPs=0.95, clambda=FALSE, naxes=4 , times = 1000) # trait = size


write.table(results, file=args[4], row.names=FALSE, col.names=TRUE, dec=".", na="NA", eol="\n", sep="\t")
lamb <- estlambda(results[,"P1df"])
write.table(lamb[[1]] , "lambda_with_correction_for_pop_str", row.names = F)
