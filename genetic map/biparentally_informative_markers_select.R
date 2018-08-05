#rm(list=ls())

#setwd("C:/Users/Roger/Documents/Postgrad/Pragya_Chaube")

# corrected 10 May 17
#
# SNPs within contigs ordered so that all proper pairs are combined
#
# removed bug that output repeated SNPs
#
# and prevented from reusing a SNP in two synthetic biparental markers
#

# get original data
# NB contig and position columns separated in excel if necessary!!
#old <- read.table("lg1_geno_allmarkers.txt",sep="\t")

args=commandArgs(trailingOnly = T)

old <- read.table(args[1] , sep="\t")

# name some columns and create some useful variables
colnames(old)[1:14]<- c("ID","male_map","female_map","error","phase","ng1","ng2","ng3","nm","contig","pos",
                        "z","male_geno","female_geno")
off <- c(15:197)
p_off <- c(13:197)
old$type <- rep("both",length(old$ID))
old$type[old$male_geno != "1 2"] <- "female"
old$type[old$female_geno != "1 2"] <- "male"

old$ngen <- rep(183,length(old$ID))  # intended to get numbers of genotypes
for (i in 1:length(old$ID)){
  old$ngen[i] <- 183-sum(old[i,off]=="0 0")
}


# set up merging function
mark_merge <- function(gen1,gen2){
  ngen <- rep("0 0",length(gen1))
  ngen[gen1 == "1 2" & (gen2 == "1 1" | gen2 == "2 2")] <- "1 2"
  ngen[gen2 == "1 2" & (gen1 == "1 1" | gen1 == "2 2")] <- "1 2"
  ngen[gen1 == "1 2" & gen2 == "1 2"] <- "2 2"
  ngen[(gen1 == "1 1" | gen1 =="2 2") & (gen2 == "1 1" | gen2 == "2 2")] <- "1 1"
  return(ngen)
}

# loop through contigs
for (con in unique(old$contig)){
  temp <- old[((substr(old$error,1,1) != "-") & (old$contig == con)),] # df for this contig, excluding markers with errors
  temp <- temp[order(temp$pos),] # order SNPs by position
  if (sum(temp$type=="both") != 0){ # if biparental marker is present
    write.table(temp[temp$type=="both",1:197],file=paste(args[1] , "bi_only" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
    write.table(temp[temp$type=="both",1:197],file=paste(args[1] , "bi_male" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
    write.table(temp[temp$type=="both",1:197],file=paste(args[1] , "bi_female" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
    write.table(temp[temp$type=="both",1:197],file=paste(args[1] , "bi_both" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
  } else if (sum(temp$type=="male") != 0 & sum(temp$type=="female") != 0) {
    x <- 0 # count synthetic markers generated
    y <- 0
    for (c in 1:(length(temp$ID)-1)){
      if (temp$type[c] != temp$type[c+1] & abs(temp$pos[c+1]-temp$pos[c])<1000
          & sum(temp[c,off]=="0 0" | temp[c+1,off]=="0 0")<92
          & y != c) { 
        # find adjacent pairs that are male+female, <1kb apart, not too many missing genotypes
        x <- x+1
        y <- c+1
        synth <- temp[c,]
        synth[p_off] <- mark_merge(temp[c,p_off],temp[c+1,p_off])
        synth$ID <- paste(synth$ID,"S",x,sep="")
        #write.table(synth[,1:197],file=paste(args[1] , "bi_only" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
        write.table(synth[,1:197],file=paste(args[1] , "bi_male" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
        write.table(synth[,1:197],file=paste(args[1] , "bi_female" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
        write.table(synth[,1:197],file=paste(args[1] , "bi_both" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
      }
    } #end SNP loop
    if (x == 0){ # if no pair meets the criteria, output one male OR one female marker
      if (temp$type[which.max(temp$ngen)]=="male"){
        write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_male" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
        write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_both" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
      } else {
        write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_female" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
        write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_both" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
      }
    }
    
  } else if (length(temp$ID>0)) { # must be only one type of marker so output the one with the most offspring genotypes
    if (temp$type[which.max(temp$ngen)]=="male"){
      write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_male" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
      write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_both" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
    } else {
      write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_female" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
      write.table(temp[which.max(temp$ngen),1:197],file=paste(args[1] , "bi_both" , sep="_") , quote = F,sep="\t",append=T,col.names = F)
    }
  } # end one type of marker condition
  
} # end contig loop
