rm (list=ls())

# analyse output of LD analysis

# get results of LD analysis
tab = read.table("ANG15_LD.txt", header=F)

names(tab) = c("replicate", "deme", "sigma_2", "sp_means", "weight")
tab = tab[is.na(tab$sigma_2) == F, ]

plot(tab$deme, tab$sigma_2)
plot(tab$deme, tab$weight)
plot(tab$deme, tab$sp_means)

tab$sigwei = tab$sigma_2 * tab$weight

sums = aggregate(tab[, c("sigwei", "weight")], by=list(tab$replicate), function(x) sum(x, na.rm=T))
sums$sigma = sqrt(sums$sigwei / sums$weight)

plot(sums$sigma)

sort(sums$sigma)[3]
median(sums$sigma)
sort(sums$sigma)[48]
