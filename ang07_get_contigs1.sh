#!/bin/bash -ve

# make one file per individual that lists all contigs with at least 10 reads mapping to them

# request memory for job (default 6G, max 72G)
#$ -l mem=12G
#$ -l rmem=12G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=27:59:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-567
##$ -tc 20

taskid=${SGE_TASK_ID}

files=`ls *_paired.sam`
file1=`echo $files | cut -d" " -f $taskid` # get focal file

# get only contigs that have at least 10 reads mapping to them
grep -v "@" $file1 | cut -f 3 | sort | uniq -c | gawk '$1>=10{print $2}' > $file1-contigs
