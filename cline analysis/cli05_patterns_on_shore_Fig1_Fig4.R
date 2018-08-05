rm (list=ls())

# look at spatial patterns along the hybrid zone (genotypes and phenotypes)


################################################################################################################
##### GET INPUT DATA ###########################################################################################

# get snail data (location of each snail along shore)
snails = read.table("ANG_SnailDataLitPath_20160702.csv", sep=",", stringsAsFactors = F, header=T)

# get cline fitting results
thr = 35.69 # var.ex threshold used to distinguish between neutral and non-neutral clines
means = read.table(paste("CLI04_general_means", thr, ".txt", sep=""), header=T, stringsAsFactors=F)
Clines_means = means[means$Type=="Cline", ] # keep only SNPs with significant clines

# get habitat data
# and make habitat df: R/S: rock/stone; B/N: barnacles/none; F/N: Fucus/none
spatial2 = read.table("ANG_spatial.csv", sep=",", stringsAsFactors = F, header=T)
spatial2$x = spatial2$x - 2000 # make coordinates comparable with snail data
spatial2$y = spatial2$y - 2000
spatial2$z = spatial2$z - 200
habitat_categories = c("RBF", "RBN", "RNF", "RNN", "SBF", "SBN", "SNF", "SNN")
habitat = spatial2[spatial2$type %in% habitat_categories, ]
habitat_sep = data.frame(do.call('rbind', strsplit(as.character(habitat$type), "", fixed=T)))
habitat$rock = habitat_sep$X1
habitat$barnacle = habitat_sep$X2
habitat$fucus = habitat_sep$X3



################################################################################################################
##### FUNCTION TO IDENTIFY COORDINATES OF A POINT ON PATH ###################################################

# function to get point for any distance on path
getpoint = function(dist){
  lengths = c()
  # get lengths of all path segments
  for (number1 in seq(1, length(x_lim)-1)){
    len = ((x_lim[number1+1] - x_lim[number1])^2 + (y_lim[number1+1] - y_lim[number1])^2)^0.5
    lengths = append(lengths, len)  
  }
  
  # go through path segments until total length exceeds focal distance
  summe = 0
  for (number1 in seq(1, length(x_lim)-1)){
    summe = summe + lengths[number1]
    if (summe > dist){
      break()
    }
  }
  
  # get point by looking at only the segment identified above
  leftpoint_x = x_lim[number1]
  leftpoint_y = y_lim[number1]
  distpoint = dist - sum(lengths[0:(number1-1)]) # distance between leftpoint and focal point
  slope = (y_lim[number1+1] - y_lim[number1]) / (x_lim[number1+1] - x_lim[number1])
  a = (distpoint^2/(1+slope^2))^0.5
  point_x = leftpoint_x + a
  point_y = leftpoint_y + slope*a
  return(c(point_x, point_y))
}


################################################################################################################
##### PLOT SAMPLING AREA (WESTRAM ET AL. 2018 FIG. 1A)  ########################################################

trans1 = 67 # path position of transition boulder field -> platform
trans2 = 85 # path position of transition platform -> continuous cliff
trans3 = mean(means$Centre[means$sel==T], na.rm = T) # average centre of non-neutral clines


# plot habitat, path and snails
plot(habitat$x, habitat$y, col=as.numeric(habitat$rock)+7, pch=16,
     xlab="x coordinate (m)", ylab="y coordinate (m)", asp=1)
# one-dimensional path through the sampling habitat:
x_lim <- c(-31,5,11,13,16,17,21,30,39,49,63,75)
y_lim <- c(-62,-20,-5,3,3,5,2.5,-1,-1,15,10,10)
points(x_lim, y_lim, col=rgb(230,159,0,maxColorValue=255), type = "l", lwd=2) # plot path
points(snails$x, snails$y, col=rgb(240,228,66, maxColorValue=255), pch=16, cex=0.5) # plot snails
# add habitat transitions
arrows(getpoint(trans1)[1], getpoint(trans1)[2]+18, getpoint(trans1)[1], getpoint(trans1)[2],
       length=0.1, col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
arrows(getpoint(trans2)[1], getpoint(trans2)[2]+4, getpoint(trans2)[1], getpoint(trans2)[2],
       length=0.1, col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
arrows(getpoint(trans3)[1], getpoint(trans3)[2]+7, getpoint(trans3)[1], getpoint(trans3)[2],
       length=0.1, col=rgb(213, 94, 0, maxColorValue = 255), lwd=2)



################################################################################################################
##### PLOT SOME CLINES (WESTRAM ET AL. 2018, FIG. 1B) ##########################################################

par(mar=c(1,6.1,0,2.1), mfrow=c(2,1), oma=c(4,0,1,2.1)) # plotting parameters

x = seq(0, 150, 0.1) # points along shore
cline_plot = function(y_label) { 
  plot(x, x, ylim=c(0,1), type = "n", ylab=y_label, yaxt='n', xaxt='n', xlab="")
  axis(2, at=seq(0,1,0.2), seq(0,1,0.2), las=1)
  lines(c(trans1,trans1), c(-1,2), col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
  lines(c(trans2,trans2), c(-1,2), col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
  lines(c(trans3,trans3), c(-1,2), col=rgb(213, 94, 0, maxColorValue = 255), lwd=2)
}

# plot shape cline
plot(seq(0, 150, 0.1), seq(0, 150, 0.1), ylim=c(0,1), type = "n", ylab="Phenotype", yaxt='n', xaxt='n', xlab="")
axis(2, at=seq(0,1,0.2), seq(0,1,0.2), las=1)
lines(c(trans1,trans1), c(-1,2), col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
lines(c(trans2,trans2), c(-1,2), col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
lines(c(trans3,trans3), c(-1,2), col=rgb(213, 94, 0, maxColorValue = 255), lwd=2)
f <- 1/(1+exp(0-4*(x-70.96)/7.15)) # ascending cline function; parameters obtained from phenotypic cline fitting
fx <- 0.32 + 0.3 * f # set ends
points(x, fx, type = "l", lwd=1, col="black") 

# plot size cline
f <- 1/(1+exp(0-4*(x-88.3)/29.0)) # ascending cline function, paramteres obtained from phenotypic cline fitting
fx <- 0.24 + 0.55 * f # set ends
points(x, 1-fx, type = "l", lwd=1, col="black")

# plot beige cline
f <- 1/(1+exp(0-4*(x-91.6)/0.5)) # ascending cline function, paramteres obtained from phenotypic cline fitting
fx <- 0 + 1 * f # set ends
points(x, 1-fx, type="l", lwd=3, col="bisque3", lty=1)

# plot black cline
f <- 1/(1+exp(0-4*(x-100)/28.5)) # ascending cline function, paramteres obtained from phenotypic cline fitting
fx <- 0 + 0.29 * f # set ends
points(x, fx, type="l", lwd=3, col="black", lty=1)

# plot banded cline
f <- 1/(1+exp(0-4*(x-94.8)/12.3)) # ascending cline function, paramteres obtained from phenotypic cline fitting
fx <- 0 + 0.27 * f # set ends
points(x, fx, type="l", lwd=3, col="black", lty="dashed")


# plot some random neutral genetic clines
cline_plot("Allele frequency\n(\"Wave\" allele)")
for (number1 in (seq(1, 10))){
  means_filt = means[means$Type=="Cline" & means$sel==F, ]
  rand = round(runif(1, 1, length(means_filt$cp)))
  focal = means_filt[rand, ]
  f <- 1/(1+exp(0-4*(x-focal$Centre)/focal$Width)) # ascending cline function as in HZAR
  fx <- focal$p_crab + focal$p_diff * f # set ends
  points(x, fx, type = "l", col=rgb(0.2,0.2,0.2,0.3), lwd=1)
}

# plot some random non-neutral genetic clines
for (number1 in (seq(1, 10))){
  means_filt = means[means$Type=="Cline" & means$sel==T, ]
  rand = round(runif(1, 1, length(means_filt$cp)))
  focal = means_filt[rand, ]
  f <- 1/(1+exp(0-4*(x-focal$Centre)/focal$Width)) # ascending cline function as in HZAR
  fx <- focal$p_crab + focal$p_diff * f # set ends
  points(x, fx, type = "l", col=rgb(213, 94, 0, maxColorValue = 255, alpha=100), lwd=1)
}

axis(1, at=seq(0,152,50), seq(0,152,50), las=1)
mtext("Distance along shore (m)", 1, 3)



################################################################################################################
##### PLOT CLINE PATTERNS ALONG SHORE (WESTRAM ET AL. 2018 FIG. 4)  ############################################

dev.off()
# make empty plot
plot(0, xlim=c(67,125), xlab="Cline centre (m)", ylab="Cline slope", xaxt='n', yaxt='n',
     cex.axis=1.5, cex.lab=1.5, ylim=c(0.02,0.7))
axis(1, at=seq(60,125,10), seq(60,125,10), las=1, cex.axis=1.5)
axis(2, seq(0,0.7,0.1), seq(0,0.7,0.1), las=1, cex.axis=1.5, cex.lab=1.5)

# clinal SNPs in nnBlocks
points(Clines_means$Centre[Clines_means$LG==6&Clines_means$sel==T&Clines_means$inv==T], Clines_means$slope[Clines_means$LG==6&Clines_means$sel==T&Clines_means$inv==T],
       col=rgb(0, 158, 115, 255, maxColorValue = 255), cex=1, pch=16)
points(Clines_means$Centre[Clines_means$LG==14&Clines_means$sel==T&Clines_means$inv==T], Clines_means$slope[Clines_means$LG==14&Clines_means$sel==T&Clines_means$inv==T],
       col=rgb(240, 228, 66, 255, maxColorValue = 255), cex=1, pch=16)
points(Clines_means$Centre[Clines_means$LG==17&Clines_means$sel==T&Clines_means$inv==T], Clines_means$slope[Clines_means$LG==17&Clines_means$sel==T&Clines_means$inv==T],
       col=rgb(204, 121, 167, 255, maxColorValue = 255), cex=1, pch=16)

# nnSNPs outside nnBlocks
points(Clines_means$Centre[Clines_means$sel==T & Clines_means$inv==F],
     Clines_means$slope[Clines_means$sel==T & Clines_means$inv==F],
     col=rgb(213, 94, 0, maxColorValue = 255), cex=0.4, pch=16)

# add lines at habitat transitions
lines(c(trans1,trans1), c(-1,1), col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
lines(c(trans2,trans2), c(-1,1), col=rgb(86, 180, 233, maxColorValue = 255), lwd=2)
lines(c(trans3,trans3), c(-1,1), col=rgb(213, 94, 0, maxColorValue = 255), lwd=2)

# legend
legend(105, 0.7, c("nnBlock LG6", "nnBlock LG14", "nnBlock LG17", "all other non-neutral SNPs"), 
       col=c(rgb(0, 158, 115, 255, maxColorValue = 255), rgb(240, 228, 66, 255, maxColorValue = 255),
             rgb(204, 121, 167, 255, maxColorValue = 255), rgb(213, 94, 0, maxColorValue = 255)),
       pch = 16, lty=NULL, y.intersp = 1.8, x.intersp=0.5, bty='n', pt.cex=c(1.5,1.5,1.5,0.8), cex=1.05)
