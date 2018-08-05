#!/bin/bash -ve

# sort, filter and index all mapped read files
# "keeping only common contigs" is a way of speeding analyses up by discarding contigs with few reads mapping from the beginning

# request memory for job (default 6G, max 72G)
#$ -l mem=3G
#$ -l rmem=3G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=04:30:00
#$ -P molecol
#$ -q molecol.q
#$ -t 1-1781
#$ -tc 30

taskid=${SGE_TASK_ID}

files=`ls *.sam` # all mapping files, e.g. ANG370_000_paired.sam, ANG370_207_S.sam, ANG370_000_L_unpaired.sam, ANG370_000_R_unpaired.sam
samples=`for n in $files;do echo $n | sed "s/\.sam//g";done` # all file names without sam ending, e.g. ANG370_000_paired etc.
sample=`echo $samples | cut -d" " -f $taskid` # get focal name

# convert individual sam to bam, keeping only common contigs; sort and index
/home/bo1awx/programs/samtools-1.3.1/samtools view -L allcontigs_sorted_common.bed -u -h -S $sample.sam | /home/bo1awx/programs/samtools-1.3.1/samtools sort -T $sample-unfiltered -o $sample.bam
	# view -u: output uncompressed bam
	# view -h: keep header
	# sort -T: write temporary files to PREFIX.nnnn.bam

# index bam file
/home/bo1awx/programs/samtools-1.3.1/samtools index $sample.bam

# filter bam file, then sort
/home/bo1awx/programs/samtools-1.3.1/samtools view -b -h -F 0x100 -q 20 $sample.bam | /home/bo1awx/programs/samtools-1.3.1/samtools sort -T $sample\_filtered -o $sample\_filtered.bam
	# view -b: output bam format
	# view -F: filter out secondary hits (-F 0x100 )
	# view -q: MAPQ filter
	# (there is no properly-paired filter because we use SE reads anyway)

# remove unfiltered versions
rm $sample.bam $sample.bam.bai
