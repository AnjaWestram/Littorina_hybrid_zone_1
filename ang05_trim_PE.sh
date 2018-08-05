#!/bin/bash -ve

# Trim PE read files. Read file names are in the format IndividualID_RunID_L1.fastq (left reads) and IndividualID_RunID_R1.fastq (right reads).
# The run ID was included because for some individuals there was a separate single-end run.

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=07:59:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-567
##$ -tc 40

taskid=${SGE_TASK_ID}

files=`ls *_L1.fastq` # all left PE read files, e.g. ANG370_000_L1.fastq
samples=`for n in $files;do echo $n | cut -d"_" -f 1,2;done | sort | uniq` # all PE samples, e.g. ANG370_000
sample=`echo $samples | cut -d" " -f $taskid` # get focal sample, e.g. ANG370_000

infileL=$sample\_L1.fastq # left reads
infileR=$sample\_R1.fastq # right reads

# trim with Trimmomatic
module add apps/java/1.7
java -jar /home/bo1awx/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $infileL $infileR $sample\_L_paired_trimmed.fastq $sample\_L_unpaired_trimmed.fastq $sample\_R_paired_trimmed.fastq $sample\_R_unpaired_trimmed.fastq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:20 MINLEN:70


