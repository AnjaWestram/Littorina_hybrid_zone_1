rm (list=ls())

# look at distribution of nnSNPs on genetic map


library("Hmisc", lib.loc="~/R/win-library/3.3")


################################################################################################################
##### INPUT & OUTPUT SETTINGS ##################################################################################

# set whether to be done per map position or 1cM interval
cM = F

thr = 35.69 # var.ex threshold used to distinguish between neutral and non-neutral clines
means = read.table(paste("CLI04_general_means", thr, ".txt", sep=""), header=T, stringsAsFactors=F)



################################################################################################################
##### BETWEEN-CONTIG CLUSTERING ################################################################################

# filter input
means = means[is.na(means$sel)==F, ]
means = means[is.na(means$av)==F, ] # only SNPs with a map position

# if looking at the 10cM scale, make 10cM bins for map positions
if (cM == T){
  means$av = ceiling((means$av+0.00000001)/10)*10
}

for (LG in seq(1,17)){#
  
  print(LG)
  
  # get data for focal LG
  focal = means[means$LG==LG,]
  
  # remove map pos with < 4 SNPs
  a = aggregate(focal$av, by=list(focal$av), function(x) length(x))
  exclude = a[a$x<4, "Group.1"]
  focal = focal[(focal$av %in% exclude)==F, ]
  
  # prepare output
  out = as.data.frame(matrix(nrow=0, ncol=4))
  names(out) = c("repl", "mean_rel", "var","sum_sel")
  per_mappos_out = c()
  
  # first loop is with original data; rest permuted
  # set number of permutations here
  for (number1 in seq(1,1000)){
    
    # get weighted average and variance of proportion of sel SNPs per map positions
    per_mappos_sel   = aggregate(focal$sel, by=list(focal$av),
                               function(x) length(x[x==T]))
    per_mappos_total = aggregate(focal$sel, by=list(focal$av),
                               function(x) length(x))
    per_mappos = merge(per_mappos_sel, per_mappos_total, by="Group.1")
    names(per_mappos) = c("av", "sel", "total")
    per_mappos$rel = per_mappos$sel / per_mappos$total # proportion of nnSNPs for each map position / bin
    per_mappos$repl = number1
    
    # output 
    out[number1, "repl"] = number1
    out[number1, "mean_rel"] = weighted.mean(per_mappos$rel, per_mappos$total)
    out[number1, "var"] = wtd.var(per_mappos$rel, per_mappos$total)
    out[number1, "sum_sel"] = sum(per_mappos$sel) # number of nnSNPs at this map position / in this bin
    
    # output that contains each map position
    per_mappos_out = rbind(per_mappos_out, per_mappos)
    
    
    # permute contigs across map positions
    perm = unique(focal[, c("Contig", "av")])
    perm$av_perm = sample(perm$av, length(perm$av), replace=F)
    focal = merge(focal, perm, by=c("Contig", "av"), all=T)
    focal$av = focal$av_perm
    focal$av_perm = NULL
  }
  
  
  # get significance cutoff - 95% quantile of permuted data
  num_repl = length(out$var)-1
  thresh = sort(out$var[out$repl!=1])[0.95*num_repl]
  sig = out$var[out$repl==1]>thresh
  print(sig)
  
  
  # plot proportion of nnSNPs per map position / bin
  plot(per_mappos_out$av[per_mappos_out$repl==1],
       per_mappos_out$sel[per_mappos_out$repl==1]/(per_mappos_out$total[per_mappos_out$repl==1]),
       col=rgb(0.5,0.5,0.5,0.5), cex=1,
       main=paste(LG), ylim=c(0,1), cex.lab=1.2, xlab="Map position",
       ylab="Prop. non-neutral SNPs", yaxt='n')
  axis(2, at=seq(0,1,0.2), seq(0,1,0.2),las=1)
  
  # colour if map position / bin is in nnBlock
  if (LG==6){
    inv = per_mappos_out[per_mappos_out$av<29.5,]
    points(inv$av[inv$repl==1],
           inv$sel[inv$repl==1]/(inv$total[inv$repl==1]),
           col=rgb(0, 158, 115, maxColorValue = 255), cex=1, pch=16)
  }
  if (LG==14){
    inv = per_mappos_out[per_mappos_out$av<12.5,]
    points(inv$av[inv$repl==1],
           inv$sel[inv$repl==1]/(inv$total[inv$repl==1]),
           col=rgb(240, 228, 66, maxColorValue = 255), cex=1, pch=16)
  }
  if (LG==17){
    inv = per_mappos_out[per_mappos_out$av>46,]
    points(inv$av[inv$repl==1],
           inv$sel[inv$repl==1]/(inv$total[inv$repl==1]),
           col=rgb(204, 121, 167, maxColorValue = 255), cex=1, pch=16)
  }
}
