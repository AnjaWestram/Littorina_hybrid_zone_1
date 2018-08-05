#!/bin/bash -ve

# map properly paired PE reads to genome

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=07:59:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-567
##$ -tc 30

taskid=${SGE_TASK_ID}

files=`ls *_L_paired_trimmed.fastq` # all paired left PE files, e.g. ANG370_000_L_paired_trimmed.fastq
samples=`for n in $files;do echo $n | cut -d"_" -f 1,2;done | sort | uniq` # all PE samples, e.g. ANG370_000
sample=`echo $samples | cut -d" " -f $taskid` # get focal sample, e.g. ANG370_000
sample_short=`echo $sample | cut -d "_" -f1` # get individual ID, e.g. ANG370 - will be used as read group ID

# now map the reads
/home/bo1awx/programs/bwa-0.7.12/bwa mem -M -R "@RG\tID:$sample_short\tSM:$sample_short" Littorina_scaffolded_PacBio_run2_7_Oct_2016_hard_masked.fasta $sample\_L_paired_trimmed.fastq $sample\_R_paired_trimmed.fastq > $sample\_paired.sam
	# -M: Mark shorter split hits as secondary (for Picard compatibility)
	# -R: read group
