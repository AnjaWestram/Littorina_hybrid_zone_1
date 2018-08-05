rm (list=ls())

# look at histograms of cline parameters in simulated and observed data


################################################################################################################
##### GET INPUT DATA ###########################################################################################

thresh = 35.69 # var.ex threshold used to distinguish between neutral and non-neutral clines
means = read.table(paste("CLI04_general_means", thresh, ".txt", sep=""), header=T, stringsAsFactors=F)
Clines_means = means[means$Type == "Cline", ]# keep only SNPs with significant clines

# get simulated data; cline fits analogous to observed data
# file "cline_snps_ANG-40.txt" is for 200 loci, sigma=1.46, 4000 generations.
sim_data = read.table("./simulations/cline_snps_ANG-40.txt", header=T, stringsAsFactors=F)
sim_data$Wave[sim_data$Wave=="a1"] = 1
sim_data$Wave[sim_data$Wave=="a2"] = 2
sim_data$Wave = as.numeric(sim_data$Wave)

# add and remove columns to simulated data as for observed data
sim_data$cp = paste(sim_data$Contig, sim_data$Position, sep="_")
sim_data$slope = sim_data$p_diff/sim_data$Width
sim_data$p_wave = sim_data$p_crab+sim_data$p_diff
sim_data$Hs = (2*sim_data$p_crab*(1-sim_data$p_crab) + 2*sim_data$p_wave*(1-sim_data$p_wave))/2
sim_data$Ht = 2 * ((sim_data$p_crab + sim_data$p_wave)/2) * (((1-sim_data$p_crab) + (1-sim_data$p_wave))/2)
sim_data$Fst = (sim_data$Ht - sim_data$Hs) / sim_data$Ht
sim_data[, c("left_break", "se_lb", "left_slope", "se_ls", "right_break", "se_rb",
             "right_slope", "se_rs", "Index")] = NULL
sim_data = sim_data[is.na(sim_data$Fst)==T | (sim_data$Fst<=1 & sim_data$Fst>=0), ]
sim_data = sim_data[is.na(sim_data$Var.Ex)==T | sim_data$Var.Ex>=0, ]

# split into simulated neutral and non-neutral ("sel") SNPs
sel_data = sim_data[sim_data$Position<=200, ]
neu_data = sim_data[sim_data$Position>200, ]



################################################################################################################
##### HISTOGRAMS OF VAR.EX PLOTS (WESTRAM ET AL. 2018 FIG 2) ###################################################

# compare var.ex between observed, and simulated neutral and selected loci
par(mar=c(4.1,3.9,4.1,0.2), mfrow=c(1,2), oma=c(2,1,1,2.1)) # plotting parameters

h = hist(Clines_means$Var.Ex, plot=F, breaks=seq(0,100,1))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, xlab="% Variance explained", col=rgb(0,1,1,0.8), main="", 
     xaxt='n', yaxt='n', ylab="")
axis(1, at=seq(0,100,20), seq(0,100,20), pos=0)
axis(2, at=seq(0,12,2), seq(0,12,2), pos=0, las=2)
title(ylab="% of SNPs", line=1.5)

h = hist(sel_data$Var.Ex, plot=F, breaks=seq(0,100,1))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=rgb(213, 94, 0, maxColorValue = 255, alpha=170), add=T,
     border=rgb(213, 94, 0, maxColorValue = 255))

h = hist(neu_data$Var.Ex, plot=F, breaks=seq(0,100,1))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=rgb(0.7,0.7,0.7,0.67), add=T, border=rgb(0.5,0.5,0.5))

lines(c(thresh,thresh), c(0,12), col=rgb(213, 94, 0, maxColorValue = 255), lwd=2)
legend(35, 12, c("simulated selected", "simulated neutral", "observed") , 
  pt.bg=c(rgb(213, 94, 0, maxColorValue = 255, alpha=170), rgb(0.7,0.7,0.7,0.67), rgb(0,1,1,0.8)),
  pch = 21, lty=NULL, y.intersp = 1.8, x.intersp=0.5, bty='n',
  col=c(rgb(213, 94, 0, maxColorValue = 255), rgb(0.5,0.5,0.5), "black"), pt.cex=1.5, cex=0.8)


# blow-up of tail of distribution to see difference between observed and simulated neutral data
h = hist(Clines_means$Var.Ex, plot=F, breaks=seq(0,100,1))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, xlab="% Variance explained", col=rgb(0,1,1,0.8), main="", 
     xaxt='n', yaxt='n', ylab="% of SNPs", ylim=c(0,0.3), xlim=c(thresh,100), xaxs='i')
axis(1, at=seq(0,100,20), seq(0,100,20), pos=0)
axis(2, at=seq(0,0.3,0.05), seq(0,0.3,0.05), pos=thresh, las=2)

h = hist(neu_data$Var.Ex, plot=F, breaks=seq(0,100,1))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=rgb(0.7,0.7,0.7,0.67), add=T, border=rgb(0.5,0.5,0.5))

lines(c(thresh,thresh), c(0,0.3), col=rgb(213, 94, 0, maxColorValue = 255), lwd=2)



################################################################################################################
##### PLOT CLINE PARAMETERS FOR NEUTRAL AND nnSNPS #############################################################

par(mar=c(4.1,4.3,4.1,0.2), mfrow=c(1,4)) # plotting parameters

h = hist(Clines_means$Width[Clines_means$sel==F], breaks=seq(0,120,4), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, xlab="Cline width (m)", col=rgb(0.7,0.7,0.7,1), main="", cex.lab=1.4, cex.axis=1.4,
     yaxt='n', ylab="% of SNPs", xaxt="n", ylim=c(0,50))
axis(1, at=seq(0,120,20), seq(0,120,20), pos=0, cex.axis=1.4, col="black")
axis(2, at=seq(0,50,10), seq(0,50,10), las=1, pos=0, cex.axis=1.4, col="black")

h = hist(Clines_means$Width[Clines_means$sel==T], breaks=seq(0,120,4), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, xlab="Cline width (m)", col=rgb(213, 94, 0, 100, maxColorValue = 255), main="", cex.lab=1.4, cex.axis=1.4,
     yaxt='n', ylab="% of non-neutral SNPs", ylim=c(0,30), add=T)

legend(7, 50, c("non-neutral SNPs", "neutral SNPs"), 
       pt.bg=c(rgb(213, 94, 0, maxColorValue = 255, alpha=170), col=rgb(0.7,0.7,0.7,1)),
       pch = 21, lty=NULL, y.intersp = 1.8, x.intersp=0.5, bty='n', pt.cex=1.5, cex=1.05)


h = hist(Clines_means$slope[Clines_means$sel==F], breaks=seq(0,1,0.05), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=rgb(0.7,0.7,0.7,1), xaxt='n', yaxt='n', main="", cex.lab=1.4, xlab="Cline slope",
     ylab="% of SNPs", ylim=c(0,50))
axis(1, at=seq(0,1,0.2), seq(0,1,0.2), pos=0, cex.axis=1.4)
axis(2, at=seq(0,50,10), seq(0,50,10), las=1, pos=0, cex.axis=1.4)

h = hist(Clines_means$slope[Clines_means$sel==T], breaks=seq(0,1,0.05), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=rgb(213, 94, 0, 100, maxColorValue = 255), xaxt='n', yaxt='n', main="", cex.lab=1.4, xlab="Cline slope",
     ylab="% of non-neutral SNPs", add=T)


h = hist(Clines_means$p_diff[Clines_means$sel==F], breaks=seq(0,1,0.05), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, xlab="Allele frequency difference\nbetween cline ends",
     col=rgb(0.7,0.7,0.7,1), freq=F, main="", xaxt='n', yaxt='n', ylim=c(0,25),
     cex.lab=1.4, ylab="% of SNPs")
axis(1, at=seq(0,1,0.2), seq(0,1,0.2), pos=0, cex.axis=1.4)
axis(2, at=seq(0,25,5), seq(0,25,5), las=1, pos=0, cex.axis=1.4)

h = hist(Clines_means$p_diff[Clines_means$sel==T], breaks=seq(0,1,0.05), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, xlab="Allele frequency difference\nbetween cline ends",
     col=rgb(213, 94, 0, 100, maxColorValue = 255), freq=F, main="", xaxt='n', yaxt='n', ylim=c(0,15),
     cex.lab=1.4, ylab="% of non-neutral SNPs", add=T)


h = hist(Clines_means$Centre[Clines_means$sel==F], breaks=seq(0,152,4), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, xlab="Cline centre (m)", col=rgb(0.7,0.7,0.7,1),
     freq=F, main="", xaxt='n', yaxt='n', ylim=c(0,70), cex.lab=1.4, ylab="% of SNPs")
axis(1, at=seq(0,150,20), seq(0,150,20), pos=0, cex.axis=1.4)
axis(2, at=seq(0,70,10), seq(0,70,10), las=1, pos=0, cex.axis=1.4)

h = hist(Clines_means$Centre[Clines_means$sel==T], breaks=seq(0,152,4), plot=F)
h$density = h$counts/sum(h$counts)*100
plot(h, xlab="Cline centre (m)", col=rgb(213, 94, 0, 100, maxColorValue = 255),
     freq=F, main="", xaxt='n', yaxt='n', ylim=c(0,60), cex.lab=1.4, ylab="% of non-neutral SNPs", add=T)
