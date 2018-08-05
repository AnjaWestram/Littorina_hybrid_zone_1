rm(list=ls())

# aim is to combine maps from bi_male and bi_female LepMap output
#
# biparentally informative markers are in both files and should have similar positions
#
# singly informative markers are only in one file and need to be inserted into the map 
# according to their position relative to bi markers in the informative parent
#
# 23 June 2017
#
# modified to handle singly-informative markers at position 0.0 in either the bi-male or bi-female map
# and to handle singly-informative markers beyond the first and last biparental marker generally
# here the 'anchor' is the first (or last) biparental marker
# 

args=commandArgs(trailingOnly=T)

male <- c("chr" , args[1] , "_male.3_markerName")
male <- read.table(paste(male, collapse=""), header=F, sep="\t")
colnames(male) <- c("contig","pos","m_m","m_f","error","phase")
male$contig <- as.character(male$contig)

female <- c("chr" , args[1] , "_female.3_markerName")
female <- read.table(paste(female , collapse =""), header=F, sep="\t")
colnames(female) <- c("contig","pos","f_m","f_f","error","phase")
female$contig <- as.character(female$contig)

# get bi markers together

bi <- merge(male,female,by = c("contig","pos"),all=F)
#plot(bi$m_m,bi$f_f)

# check for reversed order and switch male map if necessary
if (cor(bi$m_m,bi$f_f)<0){
  male$m_m <- max(male$m_m)-male$m_m
  bi$m_m <- max(bi$m_m)-bi$m_m
}

# get average positions for bi markers

bi$av <- (bi$m_m+bi$f_f)/2
#plot(bi$av,bi$m_m)
#plot(bi$av,bi$f_f)

comb <- bi[,c(1,2,11)] # add bi markers to combined map

# work through female map, finding and positioning female informative markers

for (i in 1:length(female$contig)){
  # is this SNP female informative (ie not in the bi set)
  # NB a contig should not have bi and female markers
  # and it should only have one female marker
  if (!(female$contig[i] %in% bi$contig)){
    print(female$contig[i])
    x <- female$f_f[i] # map position of this contig
    dist <- bi$f_f-x # distances of bi markers from this SNP in the female map
    if (x < min(bi$f_f)){ # if female marker is before first bi marker
      rf <- min(bi$f_f)
      ra <- min(bi$av)
      xav <- x+ra-rf
    } else if (x>max(bi$f_f)){ # if female marker is after last bi marker
      lf <- max(bi$f_f)
      la <- max(bi$av)
      xav <- x+la-lf
    } else {
      lf <- bi$f_f[dist==min(dist[dist<=0])][1] # map position of nearest bi marker to the left in female map
      rf <- bi$f_f[dist==min(dist[dist>=0])][1] # map position of nearest bi marker to the right in female map
      la <- bi$av[dist==min(dist[dist<=0])][1] # map position of nearest bi marker to the left in average map
      ra <- bi$av[dist==min(dist[dist>=0])][1] # map position of nearest bi marker to the right in average map
      xav <- la+((ra-la)/(rf-lf))*(x-lf) # map position for this SNP on the average map
      if (rf==lf){xav <- la} # to deal with cases where the focal contig maps to the same position as a bi marker
    }
    comb<-rbind(comb, data.frame(contig=female$contig[i], pos=female$pos[i], av=xav)) # add marker to combined map
  }
}

# repeat for male map, finding and positioning male informative markers

for (i in 1:length(male$contig)){
  # is this SNP male informative (ie not in the bi set)
  # NB a contig should not have bi and male markers
  # and it should only have one male marker
  if (!(male$contig[i] %in% bi$contig)){
    print(male$contig[i])
    x <- male$m_m[i] # map position of this contig
    dist <- bi$m_m-x # distances of bi markers from this SNP in the male map
    if (x < min(bi$m_m)){ # if male marker is before first bi marker
      rf <- min(bi$m_m)
      ra <- min(bi$av)
      xav <- x+ra-rf
    } else if (x>max(bi$m_m)){ # if male marker is after last bi marker
      lf <- max(bi$m_m)
      la <- max(bi$av)
      xav <- x+la-lf
    } else {
      
      lf <- bi$m_m[dist==min(dist[dist<=0])][1] # map position of nearest bi marker to the left in male map
      rf <- bi$m_m[dist==min(dist[dist>=0])][1] # map position of nearest bi marker to the right in male map
      la <- bi$av[dist==min(dist[dist<=0])][1] # map position of nearest bi marker to the left in average map
      ra <- bi$av[dist==min(dist[dist>=0])][1] # map position of nearest bi marker to the right in average map
      xav <- la+((ra-la)/(rf-lf))*(x-lf) # map position for this SNP on the average map
      if (rf==lf){xav <- la} # to deal with cases where the focal contig maps to the same position as a bi marker
    }
    comb<-rbind(comb, data.frame(contig=male$contig[i], pos=male$pos[i], av=xav)) # add marker to combined map
  }
}

# adjust combined distances to start at zero
comb$av <- comb$av-min(comb$av)

# order the combined map
comb <- comb[order(comb$av),]
pl=c("FinalPlot_chr_" , args[1] , ".png")
png(filename=paste(pl , collapse=""))
plot(comb$av , main = paste("lg" , args[1] , sep = "_") , ylab = "Averaged distances (cM)")

write.table(comb,paste("SexAveraged_LG" , args[1] , sep = "_"), sep = "\t" , quote = F , col.names = T , row.names = F)

dev.off()
