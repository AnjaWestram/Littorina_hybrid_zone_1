rm(list=ls())

# estimate dispersal based on observed LD

#install.packages("genetics")
library(genetics)

#install.packages("dplyr")
library(dplyr)

# updates 1 May 2017
#
# focus on central demes (currently 75-110m, if using first deme option) 
# where there is most information (also allows more loci to be kept)
#
# now use D not r as measure of LD
#
# still need to update how map information is input and used (line 189)
# now just need to know which chromosome the SNP is on then exclude pairs within chromosome
#
# 2 May 2017
#
# added filters: width>2 (already present),145<centre>5,p_diff>0.1,centre-width>1,centre+width<150
# added Marina's counter correction
#
# replaced old map information with lg information
#
# adjusted to use Anja's genotype file and cline file with averages over jitters
#
# 10 June 2017
#
# add random shift to deme boundaries


# read genotype file
geno0 <- read.table("ANG14.filt.GT.format",header=T)

write.table("", "ANG15_LD.txt", append=F, quote=F, row.names=F, col.names=F)


# do everything 50 x
for (number01 in seq(1,50)){
  
  keep = sample(seq(1, length(unique(geno0$CHROM))), size = 2000, replace = F)
  keep = unique(as.character(geno0$CHROM))[keep]
  geno = geno0[geno0$CHROM %in% keep, ]
  
  # read snail data
  snail <- read.csv("ANG_SnailDataLitPath_20160702.csv")
  snail <- snail[is.na(snail$DistAlongPath)==F,]
  snail <- filter(snail, snail_ID %in% colnames(geno)) # keep only snails with genotype information
  
  snails.nm <- c("CHROM","POS",as.character(snail$snail_ID))
  geno <- select(geno, one_of(snails.nm)) # keep only genotypes with snail information
  
  # get cline fit and map information 
  
  fits <- read.table("CLI04_LD_meansNA.txt",header=T)
  
  # remove non-clinal SNPs
  fits <- fits[(substr(fits$Type,4,8)=="Cline" | fits$Type == "Cline"),]
  # optionally remove loci with very narrow clines (where slope estimate for average position in 'deme' does not make sense)
  fits <- fits[fits$Width>5 & fits$Width<40 & fits$Centre>60 & fits$Centre<120 & fits$p_diff>0.2
               & (fits$Centre+fits$Width)<150 & (fits$Centre-fits$Width)>1,]
  
  #make comparable SNP names then merge with map to get df with only relevant SNPs 
  fits$SNP <- paste(fits$Contig,fits$Position, sep=".")
  geno$SNP <- paste(geno$CHROM,geno$POS, sep=".")
  
  # add LG information
  lg <- read.table("map_v11.txt",head=T)
  lg$SNP <- paste(lg$contig,lg$pos, sep=".")
  fits <- merge(fits,lg,by="SNP",all=F) # keep only SNPs with slope and map information
  fits <- fits[fits$LG != 0,] # ADDED 20170803 Anja
  geno <- filter(geno, SNP %in% fits$SNP) # keep only genotypes with fit and LG information
  fits <- merge(fits,geno,by="SNP",all=F) # then keep only those with genotypes too
  
  #remove duplicates
  fits <- fits[!duplicated(fits$SNP),]
  
  # get genotypes in columns in <genetics> package format
  snail.map <- data.frame(colnames(geno)[3:(ncol(geno)-1)])
  colnames(snail.map) <- "snail_ID"
  
  for (snp in 1:length(fits$SNP)){
    if (fits$Wave[snp]==1){
      temp <- as.character(geno[snp,3:(ncol(geno)-1)])
      temp[temp=="2"]<-"CC"
      temp[temp=="3"]<-"CW"
      temp[temp=="4"]<-"WW"
      temp[temp=="1"]<-NA
    } else {
      temp <- as.character(geno[snp,3:(ncol(geno)-1)])
      temp[temp=="2"]<-"WW"
      temp[temp=="3"]<-"CW"
      temp[temp=="4"]<-"CC"
      temp[temp=="1"]<-NA
    }
    snail.map[,snp+1] <- temp
    colnames(snail.map)[snp+1]<-geno$SNP[snp]
  }
  
  
  
  #remove genotype columns where there are fewer than 100 genotypes
  gn <- 2
  for (g in 2:length(snail.map)){
    if (sum(is.na(snail.map[,g]))<500){
      snail.map[,gn] <- snail.map[,g]
      gn <- gn+1
    }
  }
  gn <- gn-1  # Marina's correction!
  snail.map <- snail.map[,1:gn]
  
  #create a factor than groups snails into demes
  breaks <- c(0,40,50,60,65,70,75,80,85,90,95,100,110,120,130,150,160)
  # slide boundaries left or right a bit
  breaks <- breaks+runif(1,-1,1)
  #breaks <- c(0,60,70,75,85,90,95,100,110,140,160) # option with bigger groups
  nd <- length(breaks)-1
  deme <- cut(snail$DistAlongPath,breaks,labels=F)
  #deme <- cut(snail$DistAlongPath,15,labels=F)
  
  # alternative way to set up demes - maybe better because they are roughly equal in terms of numbers
  # of snails, but loses control of deme width
  # rank_dist <- rank(snail$DistAlongPath,ties="average")
  # nd <- 16
  # deme <- cut(rank_dist,breaks=nd,labels=F)
  
  
  # create a genotype matrix from the relevant columns of snail.map, using one SNP per contig
  # choosing SNP with highest VarEx
  # assumes SNPs are grouped by contig in the genotype file
  snail_gen <- data.frame(snail.map$snail_ID)
  colnames(snail_gen)[1]<- "snail_ID"
  #snail_g#en <#- merge#(snail_gen,#snail[c#("snail_ID"#,"d#eme")],by="#snail#_ID",a#ll=F)
  
  
  i <- 0
  cn <- rep(NA,gn-1)
  
  # first get polymorphic loci
  for (g in 2:gn){
    
    print (g)
    snail_gen$l1 <- genotype(snail.map[,g],sep=1)
    n <- rep(0,nd)
    af <- rep(0,nd)
    for (c in 1:nd){
      n[c] <- summary(snail_gen$l1[deme==c])[[2]][1,1]+summary(snail_gen$l1[deme==c])[[2]][2,1]
      af[c] <- summary(snail_gen$l1[deme==c])[[2]][1,2]
    }
    #only keep the locus if it is present and polymorphic in central demes (necessary for the LD function)
    #save contig and position for loci retained
    if (min(n[7:12])>2 & min(af[7:12])>0 & max(af[7:12])<1){
      names(snail_gen)[names(snail_gen)=="l1"] <- colnames(snail.map)[g]
      i <- i+1 # count up if this SNP is in a new contig
      cn[i] <- colnames(snail.map)[g]
    }
    else{
      snail_gen <- select(snail_gen,-l1)
      
    }
    
  }
  
  
  # then find SNP with highest var.ex per contig
  max.vx <- aggregate(fits$Var.Ex[fits$SNP %in% colnames(snail_gen)],by=list(fits$Contig[fits$SNP %in% colnames(snail_gen)]),FUN = max)
  names(max.vx) = c("Contig", "Var.Ex")
  max.vx = merge(fits[,c("Contig", "SNP", "Var.Ex")], max.vx, by=c("Contig", "Var.Ex"))
  max.vx = max.vx[(duplicated(max.vx$Contig))==F, ]
  
  snail_gen <- snail_gen[, names(snail_gen) %in% c("snail_ID", max.vx$SNP)]
  
  
  #make vector of SNP names
  SNP <- cn[cn %in% max.vx$SNP]
  slopes <- data.frame(SNP)
  
  #merge with fits so that same subset of loci are used throughout
  slopes <- merge(slopes,fits,by="SNP",all=F)
  
  
  # get LD matrix for each central deme and mean LD 
  ld_means <- rep(0,length(breaks)-1)
  
  for (dc in 8:12){ #central demes
    #for (dc in 9:9){  #deme 9 only
    nam <- paste("rm", dc, sep = "")
    lddc <- LD(snail_gen[deme==dc,])
    assign(nam, lddc$D)
    ld_means[dc] <- mean(lddc$D,na.rm=T)
    rm(lddc) #save space
    print(dc)
  }
  
  
  # quick look at means
  dm <- aggregate(snail$DistAlongPath, by=list(deme), FUN= "mean")[,2]
  plot(dm,ld_means)
  
  
  # calculate slopes for each central deme (i.e. at mean position of individuals)
  subsetcline <- (slopes$Type=="Cline" & !is.na(slopes$Type))
  
  for (dc in 8:12){
    slope <- rep(0,length(slopes$Index)) # start by setting slope=0 - it will stay like this if Type=NC
    nam <- paste("slope",dc,sep="")
    # slope for type Cline
    slope[subsetcline] <- slopes$p_diff[subsetcline]*4*exp(4*(dm[dc]-slopes$Centre[subsetcline])/slopes$Width[subsetcline])/
      (slopes$Width[subsetcline]*(1+exp(4*(dm[dc]-slopes$Centre[subsetcline])/slopes$Width[subsetcline]))^2)
    # slope for type lt_Cline to the right of the break
    subsetltcr <- slopes$Type=="lt_Cline" & slopes$Centre-slopes$left_break<dm[dc] & !is.na(slopes$Type) & !is.na(slopes$left_break)
    slope[subsetltcr] <- slopes$p_diff[subsetltcr]*4*exp(4*(dm[dc]-slopes$Centre[subsetltcr])/slopes$Width[subsetltcr])/
      (slopes$Width[subsetltcr]*(1+exp(4*(dm[dc]-slopes$Centre[subsetltcr])/slopes$Width[subsetltcr]))^2)
    # slope for lt_Cline to the left of the break
    subsetltcl <- slopes$Type=="lt_Cline" & slopes$Centre-slopes$left_break>dm[dc] & !is.na(slopes$Type) & !is.na(slopes$left_break)
    z <- 4*slopes$left_break/slopes$Width
    a <- 1/(1+exp(z))
    k <- slopes$left_slope*(1-a)
    slope[subsetltcl] <- (4*a[subsetltcl]*k[subsetltcl]/slopes$Width[subsetltcl])*
      exp(k[subsetltcl]*(z[subsetltcl]+4*(dm[dc]-slopes$Centre[subsetltcl])/slopes$Width[subsetltcl]))
    # slope for type rt_Cline to the left of the break
    subsetrtcl <- slopes$Type=="rt_Cline" & slopes$Centre-slopes$right_break<dm[dc] & !is.na(slopes$Type) & !is.na(slopes$right_break)
    slope[ subsetrtcl] <- slopes$p_diff[ subsetrtcl]*4*exp(4*(dm[dc]-slopes$Centre[ subsetrtcl])/slopes$Width[ subsetrtcl])/
      (slopes$Width[ subsetrtcl]*(1+exp(4*(dm[dc]-slopes$Centre[ subsetrtcl])/slopes$Width[ subsetrtcl]))^2)
    # slope for rt_Cline to the right of the break
    subsetrtcr <- slopes$Type=="rt_Cline" & slopes$Centre-slopes$right_break>dm[dc] & !is.na(slopes$Type) & !is.na(slopes$right_break)
    z <- 4*slopes$right_break/slopes$Width
    a <- 1/(1+exp(z))
    k <- slopes$right_slope*(1-a)
    slope[subsetrtcr] <- (4*a[subsetrtcr]*k[subsetrtcr]/slopes$Width[subsetrtcr])*
      exp(k[subsetrtcr]*(z[subsetrtcr]-4*(dm[dc]-slopes$Centre[subsetrtcr])/slopes$Width[subsetrtcr]))
    #print(max(slope))
    #print(slopes$SNP[which.max(slope)]) # if necessary, check that maximum slopes are reasonable
    assign(nam, slope)
  }
  
  # make matrix of recombination fractions, c
  # should be 0.5 if different chromosomes, NA otherwise
  c <- matrix(nrow=length(slopes$SNP), ncol=length(slopes$SNP))
  for (i in 1:length(slopes$SNP)){
    for (j in 1:length(slopes$SNP)){
      if (slopes$LG[i]!=slopes$LG[j]){c[i,j]<-0.5}
      else{
        c[i,j]<- NA
        
      }
    }
  }
  
  # LD for each deme (and concatenate), plus regression
  r <- NULL
  sp <- NULL
  sp_means <- rep(NA,nd) # means of slope product corrected by recombination (1+c)/c
  sigma_2 <- rep(NA,nd)
  se_sigma_2 <- rep(NA,nd)
  
  for (dc in 8:12){
    #for (dc in 9:9){   # or just for central deme
    # get a matrix of slope products, laid out like the LD matrix
    sl <- get(paste("slope",dc,sep=""))
    rm <- as.vector(get(paste("rm",dc,sep="")))
    r <- rbind(r,rm) # concatenation, if required
    
    spm <- as.vector(matrix(sl,nrow=i,ncol=1) %*% matrix(sl,nrow=1,ncol=i))
    
    sp <- rbind(sp,spm)  # concatenation, if required
    
    # get slope.product*((1+c)/c) to adjust for recombination rate
    spmr <- spm*(1+as.vector(c))/as.vector(c)
    sp_means[dc] <- mean(spmr,na.rm=T)
    
    # regress LD on slope product
    
    ld_lm <- lm(rm~0+spmr,na.action=na.omit)
    print(summary(ld_lm))
    sigma_2[dc] <- coef(ld_lm)
    se_sigma_2[dc] <- coef(summary(ld_lm))[2]
    #print(plot(spm,rm))
    print(dc)
    hist(spm, main=dc)
  }
  
  
  #plot(dm,sigma_2)
  #plot(dm,sp_means)
  
  weight <- 1/se_sigma_2
  sigma <- sqrt(sum(sigma_2*weight, na.rm = T)/sum(weight, na.rm = T))
  
  sq = sqrt(sum(sigma_2[dm>75 & dm<110]*weight[dm>75 & dm<110], na.rm = T)/
              sum(weight[dm>75 & dm<110], na.rm = T)) # central 5 demes only
  
  sp <- sp*(1+as.vector(c))/as.vector(c)
  ld_lm <- lm(as.vector(r)~0+as.vector(sp),na.action=na.omit)
  #plot(sp,r)
  
  write.table(cbind(rep(number01, length(dm)), dm, sigma_2, sp_means, weight),
              "ANG15_LD.txt", append=T, quote=F, row.names=F, col.names=F)
  
  print("sigma")
  print(sigma)
  print("sq")
  print(sq)
  print("\n")
  print(ld_lm)
  print(summary(ld_lm))
  print("dim snail_gen")
  print(dim(snail_gen))
}
