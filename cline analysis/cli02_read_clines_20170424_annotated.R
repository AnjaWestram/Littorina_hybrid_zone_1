rm (list=ls())

###################################################
# Script for fitting cline models to read-count data
# as used in Westram et al. (2018) Evolution Letters
###################################################
#
# written by Roger Butlin and Anja Westram
#
# contact r.k.butlin@shef.ac.uk or a.westram@ist.ac.at
# 
# This script was written for analysis of our data set
# It was not designed to be portable and comes with no
# guarantees! However, it should be readily adaptable
# for similar read-count data sets.
#

args = commandArgs(trailingOnly=TRUE)

# install bbmle package for maximum likelihood fitting of non-linear model
install.packages("bbmle")
library(bbmle)

# read input files 
# system specifc - need to include one dimensional distance along the transect,
# sex of individual (to check for sex linkage)
# and read counts per individual and SNP for the reference and alternate allele
# data frame 'reads' has SNPs in rows and individuals in columns: ref_count_i, alt_count_i
snail <- read.csv("ANG_SnailDataLitPath_20160702.csv")
snail <- merge(snail,read.csv("ANG_dissections_20150723.csv"),by="snail_ID", all=T)
reads <- read.table(args[1], header=T)
taskid <- args[2]
print(taskid)


# System specific: identify sex of each snail, using brood pouch and penis data
sex <- function(b, p) { 
  if(b=="Y" & p=="N" & (is.na(b)==F)) y <- "female"
  if(b=="N" & p=="Y" & (is.na(b)==F)) y <- "male"
  if((b %in% c("Y", "N"))==F | (p %in% c("Y", "N"))==F | b==p) y<-"NA"
  return(y)
}
snail$sex = apply(snail[, c("brood", "penis")], MARGIN = 1, FUN = function(x) sex(x[1], x[2]))

########################
#
# Cline functions
#
########################

# cline functions assume ascending frequency of the 'reference' allele
# ref and alt alleles are later swapped, if necessary, so that ref is more common in the wave
# ecotype (right end) than the crab ecotype (left end). 

# cline functions use the probability of reads (n_r for ref allele R, n_a for alt allele A)
# given allele frequency (px) and error (e). Sum P over 3 possible genotypes.
# Buerkle & Gompert 2012, weighted by expected frequencies.
# assumes Hardy-Weinberg equilibrium

# function for log-likelihood of cline parameters (c=centre, w=width, ends= p_crab, p_wave)
# given positions (x), read counts and error (for the locus in question, le = log(error rate))
# this incorporates the probability of reads, as above

cline_lle <- function(x,n_w,n_c,c,w,p_crab,p_diff,le){
  f <- 1/(1+exp(0-4*(x-c)/w)) # ascending cline function as in HZAR (Derryberry et al. 2014)
  fx <- p_crab+p_diff*f # set end frequencies
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  if (p_crab+p_diff>1 | p_crab<0){minusll <- minusll+4} # penalty for end frequencies outside 0,1
  return(minusll)
}

#cline fit that includes Fis
cline_lle_fis <- function(x,n_w,n_c,c,w,p_crab,p_diff,fis,le){
  f <- 1/(1+exp(0-4*(x-c)/w)) # ascending cline function as in HZAR
  fx <- p_crab+p_diff*f # set ends
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  P <- (fx^2)*(1-fis)+fx*fis # hom freq adjusted for Fis
  H <- 2*fx*(1-fx)*(1-fis)
  Q <- ((1-fx)^2)*(1-fis)+(1-fx)*fis
  Prp <- Q*exp(llAA)+H*exp(llAR)+P*exp(llRR) #probability of genotype given local frequency, fx  
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  if (p_crab+p_diff>1 | p_crab<0){minusll <- minusll+4} # penalty for end frequencies outside 0,1
  return(minusll)
}

#version for asymmetric cline, otherwise the same

as_cline_lle <- function(x,n_w,n_c,c,wc,ww,p_crab,p_diff,le){
  f <- 1/(1+exp(0-4*(x-c)/wc)) # ascending cline function as in HZAR, crab side
  f[x>c] <- 1/(1+exp(0-4*(x[x>c]-c)/ww)) # change width for wave side
  fx <- p_crab+p_diff*f # set ends
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  if (p_crab+p_diff>1 | p_crab<0){minusll <- minusll+4} # penalty for end frequencies outside 0,1
  return(minusll)
}

#version for cline with left tail (crab side), otherwise the same, dl is the break point, 

lt_cline_lle <- function(x,n_w,n_c,c,w,dl,tl,p_crab,p_diff,le){
  f <- 1/(1+exp(0-4*(x-c)/w)) # ascending cline function as in HZAR, crab side
  z <- 4*dl/w
  A <- 1/(1+exp(z))
  K <- tl*(1-A) # equations follow HZAR, tl is the ratio of the tail slope to the centre slope at c-dl
  f[x<(c-dl)] <- A*exp(K*(z+(4*(x[x<(c-dl)]-c)/w))) # exponential decline when distance < c-dl
  fx <- p_crab+p_diff*f # set ends
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  if (p_crab+p_diff>1 | p_crab<0){minusll <- minusll+4} # penalty for end frequencies outside 0,1
  return(minusll)
}

#version for cline with right tail (wave side), otherwise the same, dr is the break point, 

rt_cline_lle <- function(x,n_w,n_c,c,w,dr,tr,p_crab,p_diff,le){
  f <- 1/(1+exp(0-4*(x-c)/w)) # ascending cline function as in HZAR, crab side
  z <- 4*dr/w
  A <- 1/(1+exp(z))
  K <- tr*(1-A) # equations follow HZAR, tl is the ratio of the tail slope to the centre slope at c+dr
  f[x>(c+dr)] <- 1-A*exp(K*(z+(0-4*(x[x>(c+dr)]-c)/w))) # exponential decline when distance > c+dr
  fx <- p_crab+p_diff*f # set ends
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  if (p_crab+p_diff>1 | p_crab<0){minusll <- minusll+4} # penalty for end frequencies outside 0,1
  return(minusll)
}

#version for cline with both tails (no symmetry assumed), otherwise the same, dl,dr are the break points, 

bt_cline_lle <- function(x,n_w,n_c,c,w,dl,dr,tl,tr,p_crab,p_diff,le){
  f <- 1/(1+exp(0-4*(x-c)/w)) # ascending cline function as in HZAR, crab side
  zr <- 4*dr/w
  Ar <- 1/(1+exp(zr))
  Kr <- tr*(1-Ar) # equations follow HZAR, 
  zl <- 4*dl/w
  Al <- 1/(1+exp(zl))
  Kl <- tl*(1-Al)
  f[x>(c+dr)] <- 1-Ar*exp(Kr*(zr+(0-4*(x[x>(c+dr)]-c)/w))) # exponential decline when distance < c-dl
  f[x<(c-dl)] <- Al*exp(Kl*(zl+(4*(x[x<(c-dl)]-c)/w))) # exponential decline when distance < c-dl
  fx <- p_crab+p_diff*f # set ends
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  if (p_crab+p_diff>1 | p_crab<0){minusll <- minusll+4} # penalty for end frequencies outside 0,1
  return(minusll)
}


# linear effect of distance on frequency, used to provide initial estimate of cline end frequencies
linear_lle <- function(x,n_w,n_c,p_crab,p_wave,le){
  fx <- p_crab+(p_wave-p_crab)*(x-min(x))/(max(x)-min(x)) # straight line between ends
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  Prp[fx>1 | fx<0] <- 0.0001
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  return(minusll)
}


# and the simplest model - just fixed fx
NC_lle <- function(x,n_w,n_c,p_all,le){
  fx <- p_all 
  fx[fx==0]<- 0.001
  fx[fx==1]<- 0.999
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR) #probability of genotype given local frequency, fx
  Prp[fx>1 | fx<0] <- 0.0001
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  return(minusll)
}



# get individual IDs
snail_ID <- gsub(".1", "", colnames(reads)[seq(3, length(reads[1,]), 2)], fixed=T)   

# start data frames for output
# 'Type' defines the outcome of the cline fitting process
# (Type = SL[sex linked], Dup[duplicated], NC[frequency difference between ends<0.1],
# Cline[deltaAIC>4 vs constant frequency]), Xt_Cline[deltaAIC>4 vs Cline], NV[maf<0.05])
# main output ("out") contains lines for all loci, regardless of Type

rm(out,SLout,Dupout) # clear old output if necessary
out <- data.frame(Index=numeric(),Contig=character(),Position=numeric(),Type=character(),Wave=character(),Centre=numeric(),SEcentre=numeric(),
                  Width=numeric(),SEwidth=numeric(),
                  p_crab=numeric(),SEp_crab=numeric(),p_diff=numeric(),SEp_diff=numeric(),error_rate=numeric(),
                  Var.Ex=numeric(),left_break=numeric(),se_lb=numeric(),left_slope=numeric(),se_ls=numeric(),
                  right_break=numeric(),se_rb=numeric(),right_slope=numeric(),se_rs=numeric(),Fis=numeric(),AIC_NC=numeric(),AIC_NC_Fis=numeric(),
                  AIC_IBD=numeric(),AIC_linear=numeric(),AIC_cline=numeric(),AIC_cline_fis=numeric(),AIC_lt=numeric(),AIC_rt=numeric(),AIC_bt=numeric(),Cov.ratio=numeric())
write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=T, append=F)

SLout <- data.frame(Index=numeric(),Contig=character(),Position=numeric(),Type=character(),Wave=character())
write.table(SLout,paste("CLI02_cline_snps_ANG_SL-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=T, append=F)

Dupout <- data.frame(Index=numeric(),Contig=character(),Position=numeric(),Type=character(),Wave=character(),Hets=numeric(),Homs=numeric(),
                    Chi_Tot=numeric(),Het_p=numeric(),Bad_ratio=numeric())
write.table(Dupout,paste("CLI02_cline_snps_ANG_Dup-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=T, append=F)

##########################
#
# Main loop
#
##########################

# start main SNP loop, works through SNPs trying different models and writing output
# NB this can be slow if there are many loci!

# loop starts by exttacting data for one SNP in appropriate format:
# take line i (contig,position,base1,base2,snailID_ref_count...,snailID_alt_count) and make two columns 
# for read counts and one with line distance, saving contig and position info


for (i in 1:length(reads[,1])) {

  
  type <- "NC" # default SNP type is no change
  wave <- "a1" # start by assuming allele 1 (ref) is more common in wave than crab

  contig <- as.character(reads[i,1])
  position <- as.numeric(reads[i,2])
  
  a1_count <- as.numeric(reads[i,seq(3,length(reads[1,]),2)]) #number of reads for allele 1 (ref)
  a2_count <- as.numeric(reads[i,seq(4,length(reads[1,]),2)]) #number of reads for allele 2 (alt)
  
  # working data frame
  wf <- merge(data.frame(snail_ID,a1_count,a2_count),snail,by="snail_ID",all.x=T,all.y=F)
  # wf$DIstAlongPath is the one-dimensional transect position 
  wf <- wf[is.na(wf$DistAlongPath)==F,] #remove snails with no position info
  
  #####################################
  # start by identifying 'problem' SNPs
  #####################################
  
  # find read_count_ratio and total counts
  wf$ratio <- wf$a1_count/(wf$a1_count+wf$a2_count)
  wf$reads <- wf$a1_count+wf$a2_count
  
  # remove individuals with total reads <3
  wf <- wf[wf$reads>2, ]
  
  #get crude counts of hets and homs, and so rough overall mean allele frequency (ep)
  h1 <- length(wf$a1_count[wf$a2_count<2])
  h2 <- length(wf$a1_count[wf$a1_count<2])
  het <- length(wf$a1_count[wf$a2_count>1 & wf$a1_count>1])
  ep <- (2*h1+het)/(2*length(wf$snail_ID))
  
  # get separate heterozygote and total counts for males and females
  f_hets <- length(wf$a1_count[wf$a2_count>1 & wf$a1_count>1 & wf$sex=='female'])
  m_hets <- length(wf$a1_count[wf$a2_count>1 & wf$a1_count>1 & wf$sex=='male'])
  f_tot <- length(wf$a1_count[wf$sex=='female'])
  m_tot <- length(wf$a1_count[wf$sex=='male'])
  
  # are there more hets than expected under HWE?
  exp_tot <- 2*ep*(1-ep)*length(wf$ratio)
  chi_sq_tot <- (het-exp_tot)^2/exp_tot
  chi_sq_tot[het<exp_tot] <- 1 # set to 1 if obs hets < 2pq as expected
  if((f_hets==0) & (m_hets==0)){chi_sq_tot<-0}
  # change Type if heterozygote excess is significant at p<0.01
  if (pchisq(chi_sq_tot,1,lower.tail = F)<0.01){type <- "Dup>HWE"}
  
  # do heterozygote frequencies differ between males and females (ie is SNP sex linked)?
  exp_f <- (f_hets+m_hets)*(f_tot/(f_tot+m_tot))
  exp_m <- (f_hets+m_hets)*(m_tot/(f_tot+m_tot))
  chi_sq_m_f <- (f_hets-exp_f)^2/exp_f + (m_hets-exp_m)^2/exp_m
  if((f_hets==0) & (m_hets==0)){chi_sq_m_f<-0}
  # change Type is male-female difference is significant at p<0.01
  if (pchisq(chi_sq_m_f,1,lower.tail = F)<0.01){type <- "SL"}
  
  # is overall minor allele freq too low? If so, change Type
  if (ep<0.05 | ep>0.95){type=="NV"}
  
  # is there a central frequency peak?
  mep_c <- mean(wf$ratio[wf$DistAlongPath<max(wf$DistAlongPath)/3],na.rm=T)
  mep_h <- mean(wf$ratio[wf$DistAlongPath>(max(wf$DistAlongPath)/3) & wf$DistAlongPath<(2*max(wf$DistAlongPath)/3)],na.rm=T)
  mep_w <- mean(wf$ratio[wf$DistAlongPath>(2*max(wf$DistAlongPath)/3)],na.rm=T)
  # change Type if central part of transect has allele frequency higher (or lower)
  # than both ends by at least 0.2 (cline fitting will not work for such loci)
  if((mep_h>max(c(mep_w,mep_c))+0.2) | (mep_h<min(c(mep_w,mep_c))-0.2)){type <- "C_Peak"} 
  
  ###########################################
  # at this point 'good' loci have type=NC
  # these loci are included in cline fitting
  ###########################################
  
  
  #get 1,0 for het yes/no
  wf$het_10 <- as.numeric(wf$a1_count>1 & wf$a2_count>1)
  wf <- wf[is.na(wf$het_10)==F,] #remove snails with NAs
  
  # output info for "type=SL" loci
  if (type=="SL"){
    
    SLout <- data.frame(Index=i,Contig=contig,Position=position,Type=type,Wave=wave) 
    write.table(SLout,paste("CLI02_cline_snps_ANG_SL-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
    
    out <- data.frame(Index=i,Contig=contig,Position=position,Type=type,Wave=wave,Centre=NA,SEcentre=NA,
                      Width=NA,SEwidth=NA,
                      p_crab=NA,SEp_crab=NA,p_diff=NA,SEp_diff=NA,error_rate=NA,
                      Var.Ex=NA,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,
                      right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=NA,AIC_NC_Fis=NA,
                      AIC_IBD=NA,AIC_linear=NA,AIC_cline=NA,AIC_cline_fis=NA,AIC_lt=NA,AIC_rt=NA,AIC_bt=NA,Cov.ratio=cov.ratio)
    write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
    
  }
  
  #output info for 'Dup' loci
  if (substr(type,1,3)=="Dup"){
    Dupout <- data.frame(Index=i,Contig=contig,Position=position,Type=type,Wave=wave,
                         Hets=het,Homs=h1+h2,Chi_Tot=chi_sq_tot,Het_p=NA,Bad_ratio=NA)
    write.table(Dupout,paste("CLI02_cline_snps_ANG_Dup-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
    
    out <- data.frame(Index=i,Contig=contig,Position=position,Type=type,Wave=wave,Centre=NA,SEcentre=NA,
                      Width=NA,SEwidth=NA,
                      p_crab=NA,SEp_crab=NA,p_diff=NA,SEp_diff=NA,error_rate=NA,
                      Var.Ex=NA,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,
                      right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=NA,AIC_NC_Fis=NA,
                      AIC_IBD=NA,AIC_linear=NA,AIC_cline=NA,AIC_cline_fis=NA,AIC_lt=NA,AIC_rt=NA,AIC_bt=NA,Cov.ratio=cov.ratio)
    write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
  }
  
  # output info for 'C_Peak' loci
  if (substr(type,1,6)=="C_Peak"){ 
    
    out <- data.frame(Index=i,Contig=contig,Position=position,Type=type,Wave=wave,Centre=NA,SEcentre=NA,
                      Width=NA,SEwidth=NA,
                      p_crab=NA,SEp_crab=NA,p_diff=NA,SEp_diff=NA,error_rate=NA,
                      Var.Ex=NA,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,
                      right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=NA,AIC_NC_Fis=NA,
                      AIC_IBD=NA,AIC_linear=NA,AIC_cline=NA,AIC_cline_fis=NA,AIC_lt=NA,AIC_rt=NA,AIC_bt=NA,Cov.ratio=cov.ratio)
    write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
  }
  
  #####################
  #
  # for good loci, start model fits
  #
  #####################
  
  if(type=="NC"){
    
    #####################
    #
    # cline fitting is conducted for the observed data and then repeated 10x with transect positions jittered 
    # to check the effect of uncertainty in position estimation
    #
    #####################
    
    ### keep original data frame, because wf will be changed
    wf_original = wf
    
    for (number1 in seq(1,11)){
      
      type <- "NC" # reset
      
      # fit the 'no change' model (see bbmle vignette for further information)
      theta.init <- list(p_all=ep,le=-5)
      
      mle.NC <- mle2(NC_lle, theta.init, method="L-BFGS-B",
                     upper=list(p_all=0.999,e=-1),
                     lower=list(p_all=0.001,e=-10),
                     data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
      AIC.nc <- AIC(mle.NC)
      
      
      # fit linear model
      theta.init <- list(p_crab=ep,p_wave=ep,le=-5)
      
      mle.linear <- mle2(linear_lle, theta.init, method="L-BFGS-B",
                         upper=list(p_crab=0.999,p_wave=0.999,e=-1),
                         lower=list(p_crab=0.001,p_wave=0.001,e=-10),
                         data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
      AIC.linear<-AIC(mle.linear)
      
      # get fitted end frequencies
      pars.lin <- coef(mle.linear)
      ep_c <- as.numeric(pars.lin[1])
      ep_w <- as.numeric(pars.lin[2])
      
      # and so change in frequency
      ep_diff <- ep_w-ep_c   
      
      
      # if change is negative, swap alleles to ensure that cline is ascending, 
      # i.e a1_count is for the allele that is more common in wave
      if (ep_diff<0) {
        temp <- wf$a1_count
        wf$a1_count <- wf$a2_count
        wf$a2_count <- temp
        ep_w <- 1-ep_w
        ep_c <- 1-ep_c
        wave <- "a2"    # this indicator is output to show that the swap has happened
      }
      
      
      # if difference between ends is large enough, fit cline models
      if (abs(ep_diff)>0.1){
        
        # for the cline fit, starting values matter. 
        # therefore get rough centre and width estimates
        # centre based on heterozygote frequencies - should peak at cline centre
        
        wf$ld2 <- wf$DistAlongPath^2
        lm.het <- lm(wf$het_10 ~ wf$DistAlongPath+wf$ld2) # quadratic fit for p_het
        summary(lm.het)
        pars.het <- coef(lm.het)
        max.x <- pars.het[2]/(-2*pars.het[3]) # peak of quadratic = rough centre
        # protect against cases that give max out of range
        if ((max.x<10) | (pars.het[2]<0)){max.x<-10}
        if ((max.x>140) | ((pars.het[2]>0) & (pars.het[3]>0))){max.x<-140}
        
        # width based on total heterozygote count (assumes ends are 0,1)
        # width = 3sum(pq) for selection across env step. Use 10m demes, therefore 3*10*mean(het)/2
        sum.4pq <- 3*5*sum(aggregate(wf$het_10~as.factor(round(wf$DistAlongPath/10,0)),FUN=mean)[,2])
        
        # then do first cline fit
        theta.init <- list(c=max.x,w=sum.4pq,p_crab=ep_c+0.001,p_diff=ep_w-ep_c,le=-5) #use estimates to initialize, +/-0.001 to avoid 0 or 1
        
        mle.cline <- mle2(cline_lle, theta.init, method="L-BFGS-B",
                          upper=list(c=150,w=100,p_crab=0.9,p_diff=0.999-0.5*ep_c,le=-1),
                          lower=list(c=1,w=1,p_crab=0.001,p_diff=0.01,le=-10),       # minimum p_diff=0.01 - can be lower than 0.1 despite requirement for ep_diff>0.1
                          control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                          data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
        
        
        pars <- coef(mle.cline)
        
        # for fair comparison with tailed clines, re-run using output as starting values (tends to improve fit)
        theta.init <- list(c=pars[1],w=pars[2],p_crab=pars[3],p_diff=pars[4],le=pars[5]) 
        
        mle.cline <- mle2(cline_lle, theta.init, method="L-BFGS-B",
                          upper=list(c=150,w=1.5*min(150-pars[1],pars[1]-1),p_crab=0.9,p_diff=0.999-0.5*pars[3],le=-1),
                          lower=list(c=0.001,w=1,p_crab=0.001,p_diff=0.01,le=-10),
                          control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                          data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
        
        AIC.cline <- AIC(mle.cline)
        
        # compare cline fit to no change fit using AIC
        deltaAIC <- AIC(mle.NC)-AIC(mle.cline) # positive value means cline fit is better
        # change Type if cline fit is an improvement
        if(deltaAIC>4){type <- "Cline"}
        
        
        # test for left, right or both tails, starting from fitted values from cline_lle and tail starting at one cline width with no slope change
        # test for left tail
        theta.init <- list(c=pars[1],w=pars[2],dl=pars[2],tl=0.999,p_crab=pars[3],p_diff=pars[4],le=pars[5]) 
        
        mle.lt.cline <- mle2(lt_cline_lle, theta.init, method="L-BFGS-B",
                             upper=list(c=150,w=1.5*min(150-pars[1],pars[1]-1),dl=2*pars[2],tl=1,p_crab=0.9,p_diff=0.999-0.5*pars[3],le=-1),
                             lower=list(c=0.001,w=1,dl=0,tl=0.001,p_crab=0.001,p_diff=0.01,le=-10),
                             control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                             data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
        
        AIC.lt.cline <- AIC(mle.lt.cline) 
        
        # test for right tail
        theta.init <- list(c=pars[1],w=pars[2],dr=pars[2],tr=0.999,p_crab=pars[3],p_diff=pars[4],le=pars[5]) 
        
        mle.rt.cline <- mle2(rt_cline_lle, theta.init, method="L-BFGS-B",
                             upper=list(c=150,w=1.5*min(150-pars[1],pars[1]-1),dr=2*pars[2],tr=1,p_crab=0.9,p_diff=0.999-0.5*pars[3],le=-1),
                             lower=list(c=0.001,w=1,dr=0,tr=0.001,p_crab=0.001,p_diff=0.01,le=-10),
                             control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                             data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
        
        AIC.rt.cline <- AIC(mle.rt.cline) 
        
        # test for both tails
        theta.init <- list(c=pars[1],w=pars[2],dl=pars[2],tl=0.999,dr=pars[2],tr=0.999,p_crab=pars[3],p_diff=pars[4],le=pars[5]) 
        
        mle.bt.cline <- mle2(bt_cline_lle, theta.init, method="L-BFGS-B",
                             upper=list(c=150,w=1.5*min(150-pars[1],pars[1]-1),dl=2*pars[2],tl=1,dr=2*pars[2],tr=1,p_crab=0.9,p_diff=0.999-0.5*pars[3],le=-1),
                             lower=list(c=0.001,w=1,dl=0,tl=0.001,dr=0,tr=0.001,p_crab=0.001,p_diff=0.01,le=-10),
                             control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                             data=list(x=wf$DistAlongPath,n_w=wf$a1_count,n_c=wf$a2_count))
        
        AIC.bt.cline <- AIC(mle.bt.cline) 
        
        # compare fits, changing Type as necessary
        
        if((AIC.cline-AIC.lt.cline)>4 & (AIC.lt.cline<=AIC.rt.cline)){type <- "lt_Cline"}
        if((AIC.cline-AIC.rt.cline)>4 & (AIC.rt.cline<AIC.lt.cline)){type <- "rt_Cline"}
        if((AIC.lt.cline-AIC.bt.cline)>4 & (AIC.rt.cline-AIC.bt.cline)>4){type <- "bt_Cline"}
        
        # get output values and calculate % var for simple cline (must be better if tailed fits better)
        
        pars <- coef(mle.cline) # replace pars with the estimates from the second fit
        se <- summary(mle.cline)@coef[1:5,2] # get standard errors
        
        # estimate % var explained by the cline
        fitted <- pars[3]+pars[4]*1/(1+exp(0-4*(wf$DistAlongPath-pars[1])/pars[2]))
        var.res <- var((wf$a1_count/(wf$a1_count+wf$a2_count))-fitted)
        var.tot <- var(wf$a1_count/(wf$a1_count+wf$a2_count))
        R2 <- (var.tot-var.res)*100/var.tot
        
        # output cline fitting results
        if(type=="Cline"){
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=position,Type=type,Wave=wave,Centre=as.numeric(pars[1]),SEcentre=as.numeric(se[1]),
                            Width=as.numeric(pars[2]),SEwidth=as.numeric(se[2]),
                            p_crab=as.numeric(pars[3]),SEp_crab=as.numeric(se[3]),p_diff=as.numeric(pars[4]),SEp_diff=as.numeric(se[4]),error_rate=as.numeric(coef(mle.cline)[5]),
                            Var.Ex=R2,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=AIC.nc,AIC_NC_Fis=NA,
                            AIC_IBD=AIC.IBD,AIC_linear=AIC.linear,AIC_cline=AIC.cline,AIC_cline_fis=NA,AIC_lt=AIC.lt.cline,AIC_rt=AIC.rt.cline,AIC_bt=AIC.bt.cline,Cov.ratio=cov.ratio)
          write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
        } #end of cline output
        
        if(type=="lt_Cline"){ #lt cline output
          pars <- coef(mle.lt.cline)
          se <- summary(mle.lt.cline)@coef[1:7,2]
          
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=position,Type=type,Wave=wave,Centre=as.numeric(pars[1]),SEcentre=as.numeric(se[1]),
                            Width=as.numeric(pars[2]),SEwidth=as.numeric(se[2]),
                            p_crab=as.numeric(pars[5]),SEp_crab=as.numeric(se[5]),p_diff=as.numeric(pars[6]),SEp_diff=as.numeric(se[6]),error_rate=as.numeric(pars[7]),
                            Var.Ex=R2,left_break=as.numeric(pars[3]),se_lb=as.numeric(se[3]),left_slope=as.numeric(pars[4]),se_ls=as.numeric(se[4]),right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=AIC.nc,AIC_NC_Fis=NA,
                            AIC_IBD=AIC.IBD,AIC_linear=AIC.linear,AIC_cline=AIC.cline,AIC_cline_fis=NA,AIC_lt=AIC.lt.cline,AIC_rt=AIC.rt.cline,AIC_bt=AIC.bt.cline,Cov.ratio=cov.ratio)
          write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
        }
        
        if(type=="rt_Cline"){ #rt cline output
          pars <- coef(mle.rt.cline)
          se <- summary(mle.rt.cline)@coef[1:7,2]
          
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=position,Type=type,Wave=wave,Centre=as.numeric(pars[1]),SEcentre=as.numeric(se[1]),
                            Width=as.numeric(pars[2]),SEwidth=as.numeric(se[2]),
                            p_crab=as.numeric(pars[5]),SEp_crab=as.numeric(se[5]),p_diff=as.numeric(pars[6]),SEp_diff=as.numeric(se[6]),error_rate=as.numeric(pars[7]),
                            Var.Ex=R2,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,right_break=as.numeric(pars[3]),se_rb=as.numeric(se[3]),right_slope=as.numeric(pars[4]),se_rs=as.numeric(se[4]),Fis=NA,AIC_NC=AIC.nc,AIC_NC_Fis=NA,
                            AIC_IBD=AIC.IBD,AIC_linear=AIC.linear,AIC_cline=AIC.cline,AIC_cline_fis=NA,AIC_lt=AIC.lt.cline,AIC_rt=AIC.rt.cline,AIC_bt=AIC.bt.cline,Cov.ratio=cov.ratio)
          write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
        }
        
        if(type=="bt_Cline"){ #bt cline output
          pars <- coef(mle.bt.cline)
          se <- summary(mle.bt.cline)@coef[1:9,2]
          
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=position,Type=type,Wave=wave,Centre=as.numeric(pars[1]),SEcentre=as.numeric(se[1]),
                            Width=as.numeric(pars[2]),SEwidth=as.numeric(se[2]),
                            p_crab=as.numeric(pars[7]),SEp_crab=as.numeric(se[7]),p_diff=as.numeric(pars[8]),SEp_diff=as.numeric(se[8]),error_rate=as.numeric(pars[9]),
                            Var.Ex=R2,left_break=as.numeric(pars[3]),se_lb=as.numeric(se[3]),left_slope=as.numeric(pars[4]),se_ls=as.numeric(se[4]),
                            right_break=as.numeric(pars[5]),se_rb=as.numeric(se[5]),right_slope=as.numeric(pars[6]),se_rs=as.numeric(se[6]),Fis=NA,AIC_NC=AIC.nc,AIC_NC_Fis=NA,
                            AIC_IBD=AIC.IBD,AIC_linear=AIC.linear,AIC_cline=AIC.cline,AIC_cline_fis=NA,AIC_lt=AIC.lt.cline,AIC_rt=AIC.rt.cline,AIC_bt=AIC.bt.cline,Cov.ratio=cov.ratio)
          write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
        }
        
        # output null result if no cline fit was an improvement over 'no change' model
        if(type=="NC"){
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=position,Type="NC>0.1",Wave=wave,Centre=NA,SEcentre=NA,
                            Width=NA,SEwidth=NA,
                            p_crab=NA,SEp_crab=NA,p_diff=NA,SEp_diff=NA,error_rate=NA,
                            Var.Ex=NA,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,
                            right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=AIC.nc,AIC_NC_Fis=NA,
                            AIC_IBD=AIC.IBD,AIC_linear=AIC.linear,AIC_cline=NA,AIC_cline_fis=NA,AIC_lt=NA,AIC_rt=NA,AIC_bt=NA,Cov.ratio=cov.ratio)
          write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
        }
        
      } #end of loop for 'good' loci with allele frequency difference >0.1
      
      
      else{
        # output for allele frequency difference <0.1
        out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=position,Type=type,Wave=wave,Centre=NA,SEcentre=NA,
                          Width=NA,SEwidth=NA,
                          p_crab=NA,SEp_crab=NA,p_diff=NA,SEp_diff=NA,error_rate=NA,
                          Var.Ex=NA,left_break=NA,se_lb=NA,left_slope=NA,se_ls=NA,
                          right_break=NA,se_rb=NA,right_slope=NA,se_rs=NA,Fis=NA,AIC_NC=AIC.nc,AIC_NC_Fis=NA,
                          AIC_IBD=AIC.IBD,AIC_linear=AIC.linear,AIC_cline=NA,AIC_cline_fis=NA,AIC_lt=NA,AIC_rt=NA,AIC_bt=NA,Cov.ratio=cov.ratio)
        write.table(out,paste("CLI02_cline_snps_ANG-", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append = T)
      }
      
      ######################
      #
      # jittering
      #
      ######################
      
      # set to orginal data frame again
      wf = wf_original
      
      # generate jittered distances (Gaussian, standard devaition = 0.5m)
      wf$DistAlongPath = rnorm(length(wf_original$DistAlongPath), mean = wf_original$DistAlongPath, sd = 0.5)
      
    } # end of "repeat 10x with jittered positions" loop
    
  } # end of if (type=NC), i.e. 'good' loci
  
  print(c(as.character(i),ep,type)) # print counter, allele frequency and Type, to monitor progress
  
} # end of snp loop

