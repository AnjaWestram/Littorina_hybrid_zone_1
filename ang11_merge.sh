#!/bin/bash -ve

# merge all files from the same individual

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:59:00
##$ -P molecol
##$ -q molecol.q
#$ -t 1-567
##$ -tc 25

taskid=${SGE_TASK_ID}

files=`ls *dedup.bam` # all dedup files, e.g. ANG370_000_pairedNoOverl_dedup.bam, ANG370_207_S_dedup.bam, ANG370_000_L_unpaired_dedup.bam   ANG370_000_R_unpaired_dedup.bam - careful: PE files without overlap clipping must be removed before running this
samples=`for n in $files;do echo $n | cut -d"_" -f 1;done | sort | uniq` # all individuals, e.g. ANG370
sample=`echo $samples | cut -d" " -f $taskid` # get focal individual, e.g. ANG370

# merge all bam files for this individual
all=`ls $sample\_*bam`
/home/bo1awx/programs/samtools-1.3.1/samtools merge -c $sample.bam $all
	# -c: output only the first @RG header
