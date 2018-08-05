#!/bin/bash -ve

# trim SE read files

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=03:59:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-80
##$ -tc 20

taskid=${SGE_TASK_ID}

files=`ls *_S.fastq` # all SE read files, e.g. ANG370_207_S.fastq
samples=`for n in $files;do echo $n | cut -d"_" -f 1,2;done | sort | uniq` # all SE samples, e.g. ANG370_207
sample=`echo $samples | cut -d" " -f $taskid` # get focal sample, e.g. ANG370_207

infile=$sample\_S.fastq

# trim with Trimmomatic
module add apps/java/1.7
java -jar /home/bo1awx/programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $infile $sample\_S_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:20 MINLEN:70
