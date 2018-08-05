#!/bin/bash -ve

# get list of commonly covered contigs

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=27:59:00
##$ -P molecol
##$ -q molecol.q

# get contigs that appear at least 10 times from all individuals
cat *_paired.sam-contigs>allcontigs.txt

# sort by number of occurrences, get number of files (individuals) where contig occurs at least 10 times and sort by that
sort allcontigs.txt | uniq -c | sort -nr > allcontigs_sorted.txt

# get only contigs with sufficient reads in >= 150 files (CHANGE NUMBER HERE IF NECESSARY)
head -n 56293 allcontigs_sorted.txt | sed "s/ .\+ //g" > allcontigs_sorted_common.txt
