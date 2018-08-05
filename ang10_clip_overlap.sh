#!/bin/bash -ve

# remove overlap between PE reads from one of the two reads (only for paired files...)

# request memory for job (default 6G, max 72G)
#$ -l mem=3G
#$ -l rmem=3G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:30:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-567
##$ -tc 25

taskid=${SGE_TASK_ID}

files=`ls *_paired_dedup.bam` # all PE dedup files, e.g. ANG370_000_paired_dedup.bam
samples=`for n in $files;do echo $n | cut -d"_" -f 1,2;done | sort | uniq` # all samples, e.g. ANG370_000
sample=`echo $samples | cut -d" " -f $taskid` # get focal sample, e.g. ANG370_000

/home/bo1awx/programs/bamUtil_1.0.13/bamUtil/bin/bam clipOverlap --in $sample\_paired_dedup.bam --out $sample\_pairedNoOverl_dedup.bam --stats
