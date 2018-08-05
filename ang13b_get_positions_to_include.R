rm (list=ls())

# get positions to include in analyses of hybrid zone data


## set maximum distance from mapped position for including SNP
thresh = 1000



################################################################################################################
##### INPUT ####################################################################################################

# get LG information
map = read.table("LinkageGroup_20170405.txt", header=T, stringsAsFactors=F)
spl = data.frame(do.call('rbind', strsplit(as.character(map$Contigs), ".", fixed=T)), stringsAsFactors=F)
map$contig = spl$X1
map$pos = as.numeric(spl$X2)
map = map[map$LG != 0, ] # keep only markers placed on map

# find cases where a contig is associated with multiple LGs, and remove those contigs
per_contig = aggregate(map$LG, by=list(map$contig), function(x) length(unique(x)))
rem = per_contig[per_contig$x>1, "Group.1"]
map = map[(map$contig %in% rem)==F, ] # keep only contigs associated with a single LG

# get all ANG SNPs
SNPs = read.table("ANG13b_all_snps_unique.txt", stringsAsFactors=F, header=F)
names(SNPs) = c("contig", "pos")



################################################################################################################
##### IDENTIFY ALL SNPS ON MAP OR CLOSE TO MAPPED SNPS #########################################################

SNPs = SNPs[SNPs$contig %in% map$contig, ] # first exclude SNPs in contigs not on map at all

# for each SNP, get closest distance to SNP on map
for (number1 in seq(1, length(SNPs$pos))){
   contig = SNPs[number1, "contig"]
   pos = SNPs[number1, "pos"]
   focal = map[map$contig == contig, ] # all SNPs in this contig that are on the map
   SNPs[number1, "mindist"] = min(abs(pos - focal$pos)) # get smallest distance of THIS SNP from mapped SNP
   print (number1)
}

to_use = SNPs[SNPs$mindist<=1000 & (is.na(SNPs$mindist)==F), ]

write.table(to_use[, c("contig", "pos")], "ANG13b_positions_to_keep1.txt", append=F, quote=F, row.names=F, col.names=F)
