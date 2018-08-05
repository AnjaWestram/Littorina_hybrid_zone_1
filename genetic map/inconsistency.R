# version 1.0
# This script will read Roger's sex- sepcific maps
# Combine with the contig names and pos
# calculate the inconsistencies across a contig
# if inconsistencies are greater than 5 cM ( and ~ possible genotype errors) - report them and updates their LMids in the 'toremove' files

#####
# files needed: map, 'modified files' , existing 'toremove' marker file (if any) (!all in the working directory!)
# args: female/male , map file , 'modified' file , 'toRemove' file , output file - new 'toRemove' file
# 'modified file' : with linkage group, contig, pos, LMid information
# 'map' file : LM2 output files (for single chromosome)
# 'toRemove' file : LMids of the reported SNPs in a single line- to be used in LM2 position re-estimations

args = commandArgs(trailingOnly=TRUE)

library(dplyr)

mod <- read.table(args[3] , sep = "\t" ,  header = F , stringsAsFactors = F)
map <- read.table(args[2] , sep = "\t" , header = F , stringsAsFactors = F , comment.char = "#" , fill = T )
remove <- read.table(args[4] , sep = " " , h = F , stringsAsFactors = F)

map_merge <- merge(map , mod , by.x = "V1" , by.y = "V1" , sort = F) # merging the mod file - names of contig , pos known now

if(args[1] == "female"){
  map_merge2 <- map_merge[ , c(1,16,17,3,4)] # selecting only female position (Lmid , contig , pos, female map position, phase/error)
} else if(args[1] == "male"){
  map_merge2 <- map_merge[ , c(1,16,17,2,4)] # selecting only male position (Lmid , contig , pos, female map position, phase/error)
} else{
  stop("Pleae specify if this is a male or a female map")
}

names(map_merge2) <- c("Lmid" , "contig" , "pos" , "cM" , "phase")

##########################################################################################################################################
### look for inconsistent map positions

inconsist_f <- data.frame()
for(i in unique(map_merge2$contig)){
  x = map_merge2[map_merge2$contig == i, ]
  if(nrow(x)>1){
    x <- x[!duplicated(x[,4]),]
  }
  else{ # this will skip the markers which are standalone; ie, single contig. 
    next()
  }
  inconsist_f <- rbind(inconsist_f, x)
}

inconsist_f <- inconsist_f[!duplicated(inconsist_f[,c("contig","pos")]),]

fault <- inconsist_f[inconsist_f$contig %in% names(which(table(inconsist_f$contig)>1)),]

fault.1 <- map_merge2[map_merge2$contig %in% fault$contig , ] # extracting all the markers with such 'inconsistency errors'

fault.1$dif_cM <- ave(fault.1$cM , fault.1$contig , FUN = function(x) c(0 , diff(x))) # map position differences of markers on he same contig
fault.1$numCont <- ave(fault.1$contig , fault.1$contig , FUN= function(x) table(x)) # number of SNPs per contig
fault.1$numcM <- ave(fault.1$cM , fault.1[,c("contig" , "cM")],  FUN= function(x) table(x)) # number of SNPs with same map position per contig
############################################
# search for ~ genotype errors, only single numcM position for probable erroneous
genotypeErrors_probable_m <- data.frame()
for(i in unique(fault.1$contig)){
  x= fault.1[fault.1$contig == i ,]
  if(nrow(x) > 1)
  { y= x[x$numcM == 1 ,]}
  else
  {y= x}
  
  genotypeErrors_probable_m= rbind(genotypeErrors_probable_m, y)
} #   record their LMid and add to the 'toRemove' files

genotypeErrors_probable_m.2 <- subset(genotypeErrors_probable_m , dif_cM > 4.99) # selecting only those genotype problems which are more than 5 cM away

remove_new=data.frame()
if(nrow(genotypeErrors_probable_m.2) >= 1){
  remove_new <- cbind(remove , t(genotypeErrors_probable_m.2$Lmid))}


#############################################################################################################
crypt_errors <- fault.1[!(fault.1$Lmid %in% genotypeErrors_probable_m$Lmid),]

#   select which ones are greater than 5 cM - remove these ~ toRemove files
#   for the smaller than 5 cM if difference between numCont and numcM == 1, then okay, else problematic
#   the remaining should be averaged

#   if the difference is more than 5 cM for more than 1 marker in a contig - remove the entire contig- add all markers from the contig in the remove file
prob_5 <- subset(crypt_errors , dif_cM > 4.99)
prob_5_conts <- fault.1[fault.1$contig %in% prob_5$contig , ]

prob_5_conts.4 <- prob_5_conts[(prob_5_conts$Lmid %in% genotypeErrors_probable_m$Lmid),] # capturing those singular errors which escape the first problem , for eg "contig66426" lg10
if(nrow(prob_5_conts.4) > 0){
  if(length(remove_new) > 0){
    remove_new <- cbind(remove_new , t(prob_5_conts.4$Lmid))
  } else{
    remove_new <- cbind(remove , t(prob_5_conts.4$Lmid))
  }
}

prob_5_conts <- prob_5_conts[!(prob_5_conts$Lmid %in% genotypeErrors_probable_m$Lmid),] # removing the ones which are captured due to genotype errors issues; for eg "contig66426" lg10

if(nrow(prob_5_conts) > 0){
  prob_5_conts.2 = data.frame()
  for(i in unique(prob_5_conts$contig)){
    x=prob_5_conts[prob_5_conts$contig == i , ]
    x$numMarkers = as.numeric(x$numCont) - as.numeric(x$numcM)
    prob_5_conts.2 = rbind(prob_5_conts.2 , x)}
  prob_5_conts.3=subset(prob_5_conts.2 , numMarkers != 1) # retains only contigs with multiple markers with several map positions and diff bigger than 5 cM
  if(length(remove_new) >0){ 
    remove_new <- cbind(remove_new , t(prob_5_conts.3$Lmid)) # all markers from a contig with more than one 5 cM in the remove file
  } else{
    remove_new <- cbind(remove , t(prob_5_conts.3$Lmid))
  }
} 


#######################################################################################################################################
# write the new files : (i) the new remove file (ii) the fault file (iii) the prob > 5 file (if applicable)

if(length(remove_new) > 0){
  write.table(remove_new , args[5] , sep = " " , row.names = F , col.names = F , quote = F)
} else{
  print("No more 'errors' to remove")
  write.table(remove , args[5] , sep = " " , row.names = F , col.names = F , quote = F)
}
write.table(fault.1 , paste(args[2] , "fault" , sep = "_") , sep = "\t" , quote = F , row.names = F)
if((nrow(prob_5) >= 1) & (nrow(prob_5_conts) >=1)){
  write.table(rbind(prob_5 , prob_5_conts) , paste(args[2] , "prob_5" , sep = "_") , sep = "\t" , quote = F , row.names = F)}

