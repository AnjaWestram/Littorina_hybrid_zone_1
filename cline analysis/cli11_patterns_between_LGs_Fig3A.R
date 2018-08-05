rm (list=ls())

# look at distribution of nnSNPs across linkage groups


library("Hmisc", lib.loc="~/R/win-library/3.3") # contains function to get weighted variance



################################################################################################################
##### GET CLINE FITTING RESULTS ################################################################################

thr = 35.69  # var.ex threshold used to distinguish between neutral and non-neutral clines
means = read.table(paste("CLI04_general_means", thr, ".txt", sep=""), header=T, stringsAsFactors=F)
means = means[(is.na(means$av)==F), ] # remove unmapped SNPs



################################################################################################################
##### PROPORTIONS OF nnSNPs PER LG #############################################################################

par(mar=c(4.1,6.3,4.1,0.2)) # plotting parameters

# summarise number and proportion of non-neutral SNPs per LG
per_lg_sel= aggregate(means$sel, by=list(means$LG), function(x) length(x[x==T])) # number of non-neutral SNPs per LG
per_lg_total = aggregate(means$sel, by=list(means$LG), function(x) length(x)) # total number of SNPs per LG
per_lg = merge(per_lg_sel, per_lg_total, by="Group.1")
names(per_lg) = c("LG", "sel", "total")
per_lg$rel = per_lg$sel / per_lg$total # proportion of nnSNPs per LG

# add proportion of non-neutral SNPs that are NOT in a putative inversion
per_lg_out = aggregate(means$inv[means$sel==T], by=list(means$LG[means$sel==T]), function(x) length(x[x==F]))
names(per_lg_out) = c("LG", "selOut")
per_lg = merge(per_lg, per_lg_out, by="LG", all=T)
per_lg$relOut = per_lg$selOut/per_lg$total # proportion of total SNPs that are non-neutral and NOT in putative inversion

# plot proportions of nnSNPs per LG (Westram et al. 2018 Fig. 3A)
pl = barplot(per_lg$rel, yaxt='n', ylim=c(0,0.15),
        col=c("white","white","white","white","white",rgb(0, 158, 115, maxColorValue = 255),
              "white","white","white","white","white","white","white",
              rgb(240, 228, 66, maxColorValue = 255),
              "white","white",
              rgb(204, 121, 167, maxColorValue = 255)),
        ylab="Number of nnSNPs /\ntotal number of SNPs", xlab="Linkage group",cex.lab=1.4)
barplot(per_lg$relOut, yaxt='n', ylim=c(0,0.14), col="white", add=T) # proportion that is non-neutral but outside inversion
axis(1, pl, labels=rep("", 17), cex.axis=0.7)
axis(2, seq(0,0.15,0.05), seq(0,0.15,0.05), las=2)
mtext(c("","","","","",6,"","","","","","","",14,"","",17), side=3, at=pl, cex=1)



################################################################################################################
##### PERMUTATION TEST #########################################################################################

# prepare output
out = as.data.frame(matrix(nrow=0, ncol=4))
names(out) = c("repl", "mean_rel", "var", "sum_sel")

means1 = means # data that will be permuted

# first loop is with original data; rest permuted (1000 permutations)
for (number1 in seq(1,1001)){
  
  # summarise per LG
  per_lg_sel   = aggregate(means1$sel, by=list(means1$LG), function(x) length(x[x==T]))
  per_lg_total = aggregate(means1$sel, by=list(means1$LG), function(x) length(x))
  per_lg = merge(per_lg_sel, per_lg_total, by="Group.1")
  names(per_lg) = c("LG", "sel", "total")
  per_lg$rel = per_lg$sel / per_lg$total
  
  # output 
  out[number1, "repl"] = number1
  out[number1, "mean_rel"] = weighted.mean(per_lg$rel, per_lg$total)
  out[number1, "var"] = wtd.var(per_lg$rel, per_lg$total)
  out[number1, "sum_sel"] = sum(per_lg$sel) # number of non-neutral SNPs across all LG (should always be the same)
  
  # permute SNPs across map positions
  means1$LG = sample(means$LG, length(means$LG), replace=F)
  
  print(number1)
}
  

# get significance cutoff - 95% quantile of permuted data
num_repl = length(out$var)-1
thresh = sort(out$var[out$repl!=1])[0.95*num_repl]
sig = out$var[out$repl==1]>thresh
print(sig)