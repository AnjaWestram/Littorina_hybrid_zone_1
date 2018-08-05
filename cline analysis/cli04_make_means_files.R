rm (list=ls())

# sum up results of cline analysis over jittered replicates for each SNP,
# add a few statistics and clean up



################################################################################################################
##### INPUT ####################################################################################################

# decide whether to output for general analyses (include jittered replicates that are tailed)
# or for LD analysis (keep only jittered replicates that are standard (non-tailed) clines)
outp = "general" # LD OR general

thresh = 35.69 # threshold to define non-neutral SNPs



################################################################################################################
##### INPUT & OVERVIEW #########################################################################################

results = read.table(paste("CLI03_clines.txt", sep=""), header=T, stringsAsFactors=F)
results$Wave[results$Wave=="a1"] = 1
results$Wave[results$Wave=="a2"] = 2
results$Wave = as.numeric(results$Wave)

results$cp = paste(results$Contig, results$Position, sep="_")

((((is.numeric(results$Var.Ex)==is.numeric(results$Width))==is.numeric(results$p_diff))==
  is.numeric(results$p_crab))==is.numeric(results$Centre))==T

# total number of SNPs:
length(unique(results$cp))
# number of Dup>HWE:
length(results[results$Type=="Dup>HWE", "cp"])
# number of SL:
length(results[results$Type=="SL", "cp"])
# number of C_Peak:
length(results[results$Type=="C_Peak", "cp"])



################################################################################################################
##### FILTERING OF RESULTS ########################################################################

# keep only Cline and NC SNPs
results_filt = results[(results$Type %in% c("Dup>HWE", "SL", "C_Peak"))==F, ]
unique(results_filt$Type)

# make table which for each SNP says how many of the jittered replicates are of each type
tab = as.data.frame.matrix(table(results_filt[, c("cp", "Type")]))

# categorise each SNP depending on whether one category is most common
# include tailed cline replicates for general downstream analyses but not for LD calculation
tab[(tab$NC + tab$"NC>0.1") >= 8, "cat"] = "NC"
if (outp == "LD"){
  tab[(tab$Cline) >= 8, "cat"] = "Cline" 
}
if (outp == "general"){
  tab[(tab$Cline + tab$rt_Cline + tab$lt_Cline + tab$bt_Cline) >= 8, "cat"] = "Cline"
}
length(tab[(is.na(tab$cat))==T, "cat"]) # unclear loci

table(tab$cat[is.na(tab$cat)==F])



##### MAKE SEPARATE DATA SETS FOR CLINAL AND NC SNPs ###########################################################
##### AVERAGE ACROSS JITTERED REPLICATES #######################################################################

# only clinal SNPs, and for those only the clinal jittered replicates
Clines_cp = row.names(tab[tab$cat == "Cline" & (is.na(tab$cat)==F), ])
Clines = results[results$cp %in% Clines_cp, ]
if (outp == "LD"){
  Clines = Clines[Clines$Type == "Cline", ]
}
if (outp == "general"){
  Clines = Clines[Clines$Type %in% c("Cline", "rt_Cline", "lt_Cline", "bt_Cline"), ]
}


# add some statistics
Clines$slope = Clines$p_diff/Clines$Width
Clines$p_wave = Clines$p_crab+Clines$p_diff
Clines$Hs = (2*Clines$p_crab*(1-Clines$p_crab) + 2*Clines$p_wave*(1-Clines$p_wave))/2
Clines$Ht = 2 * ((Clines$p_crab + Clines$p_wave)/2) * (((1-Clines$p_crab) + (1-Clines$p_wave))/2)
Clines$Fst = (Clines$Ht - Clines$Hs) / Clines$Ht

# get average over replicate jittered cline fittings
Clines_means = aggregate(Clines, by=list(Clines$cp), FUN = function(x) mean(x, na.rm=T))
Clines_means$cp=NULL
names(Clines_means)[1] = "cp"
spl = data.frame(do.call('rbind', strsplit(as.character(Clines_means$cp), "_", fixed=T)), stringsAsFactors=F)
Clines_means$Contig = spl$X1
Clines_means$Position = as.numeric(spl$X2)
Clines_means$Type = "Cline"


# only NC SNPs
NC_cp = row.names(tab[(tab$cat == "NC") & (is.na(tab$cat)==F), ])
NC = results[results$cp %in% NC_cp, ]
NC = NC[NC$Type %in% c("NC", "NC>0.1"), ]

# add some statistics
NC$p_diff = NA
NC$slope = NA
NC$p_wave = NA
NC$Hs = NA
NC$Ht = NA
NC$Fst = NA

# get average over replicates - this is not really relevant for NC data, except for AICs
NC_means = aggregate(NC, by=list(NC$cp), FUN = function(x) mean(x, na.rm=T))
NC_means$cp=NULL
names(NC_means)[1] = "cp"
spl = data.frame(do.call('rbind', strsplit(as.character(NC_means$cp), "_", fixed=T)), stringsAsFactors=F)
NC_means$Contig = spl$X1
NC_means$Position = as.numeric(spl$X2)
NC_means$Type = "NC"
NC_means$Wave = NA


# merge clinal and NC data
means = rbind.data.frame(Clines_means, NC_means)
means$Index = means$left_break = means$se_lb = means$left_slope = means$se_ls = means$right_break = means$se_rb =
  means$right_slope = means$se_rs = NULL



################################################################################################################
##### ADD OTHER INFO ###########################################################################################

# add info about selection (based on threshold obtained from simulations)
means$sel = F
means$sel[means$Var.Ex>thresh] = T

# get map data
mapdata = read.table("./map/map_v11.txt", header=T, stringsAsFactors = F)
mapdata$cp = paste(mapdata$contig, mapdata$pos, sep="_")

# add map position for cline analysis SNPs that have one
means = merge(means, mapdata[, c("cp", "LG", "av")], by="cp", all.x=T, all.y=F)

# assign SNPs not on map to closest map position, if within 1000bp
for (number1 in seq(1,length(means$cp))){
  if(is.na(means[number1, "av"])){
    contig = means[number1, "Contig"]
    pos = means[number1, "Position"]
    focal = mapdata[mapdata$contig==contig, ]
    closest = focal[abs(focal$pos-pos) == min(abs(focal$pos-pos)), ][1,]

    if((abs(closest$pos-pos))<=1000 & (is.na(closest$pos)==F)){
      means[number1, "LG"] = closest$LG
      means[number1, "av"] = closest$av
    }
  print(number1)
  }
}

# check whether there are contigs associated with multiple LGs - should be none
a = aggregate(means$LG[is.na(means$LG)==F], by=list(means$Contig[is.na(means$LG)==F]), function(x) length(unique(x)))
hist(a$x, breaks=seq(0,10))


# add inversion info
means$inv = F
means[  (means$LG==6 & means$av<29.5) & (is.na(means$LG)==F)|
        (means$LG==14 & means$av<12.5) & (is.na(means$LG)==F)|
        (means$LG==17 & means$av>46) & (is.na(means$LG)==F), "inv"] = T



################################################################################################################
##### OUTPUT ###################################################################################################

means = means[is.na(means$Fst)==T | means$Fst<=1, ]
means = means[is.na(means$Var.Ex)==T | means$Var.Ex>=0, ]

write.table(means, paste("CLI04_", outp, "_means", thresh, ".txt", sep=""), append = F, quote=F, col.names=T,
            row.names=T)
