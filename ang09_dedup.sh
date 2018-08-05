#!/bin/bash -ve

# remove PCR duplicates

# request memory for job (default 6G, max 72G)
#$ -l mem=16G
#$ -l rmem=16G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=12:30:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-1781
##$ -tc 25

taskid=${SGE_TASK_ID}

files=`ls *_filtered.bam` # all filtered bam files, e.g. ANG370_000_paired_filtered.bam, ANG370_207_S_filtered.bam, ANG370_000_L_unpaired_filtered.bam, ANG370_000_R_unpaired_filtered.bam
samples=`for n in $files;do echo $n | sed "s/\_filtered\.bam//g";done` # all file names without ending, e.g. ANG370_000_paired
sample=`echo $samples | cut -d" " -f $taskid` # get focal name

# get java (needed for running Picard)
module add apps/java/1.7

# sort with Picard
java -Xmx10g -jar /home/bo1awx/programs/picard-tools-1.129/picard.jar SortSam INPUT=$sample\_filtered.bam OUTPUT=$sample\_sorted.bam SORT_ORDER=coordinate 

# remove potential PCR duplicates
java -Xmx8g -jar /home/bo1awx/programs/picard-tools-1.129/picard.jar MarkDuplicates INPUT=$sample\_sorted.bam OUTPUT=$sample\_dedup.bam METRICS_FILE=$sample.metrics READ_NAME_REGEX=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=True
	# REMOVE_DUPLICATES: remove duplicates, rather than just marking them
	# READ_NAME_REGEX: used to estimate amount of optical duplicates, and thereby true library size. If set to null, no optical duplicate detection.
	# ASSUME_SORTED=True: file has been sorted before
	# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: reduce the number of concurrently open files (trading off execution speed)
