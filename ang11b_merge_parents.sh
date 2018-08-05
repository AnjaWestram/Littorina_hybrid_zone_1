#!/bin/bash -ve

# merge replicates for the two Crab family parents (2 replicates each)

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:59:00

module add apps/java/1.7

# merge parent replicates
/home/bo1awx/programs/samtools-1.3.1/samtools merge -c 56-female.bam.tmp 56-F.bam 56-F2.bam
/home/bo1awx/programs/samtools-1.3.1/samtools merge -c 56-male.bam.tmp 56-M.bam 56-M2.bam

# give the same RG to every read
java -jar /home/bo1awx/programs/picard-tools-1.129/picard.jar AddOrReplaceReadGroups I=56-female.bam.tmp O=56-female.bam RGID=56-female RGLB=56-female RGPL=illumina RGPU=56-female RGSM=56-female VALIDATION_STRINGENCY=SILENT
java -jar /home/bo1awx/programs/picard-tools-1.129/picard.jar AddOrReplaceReadGroups I=56-male.bam.tmp O=56-male.bam RGID=56-male RGLB=56-male RGPL=illumina RGPU=56-male RGSM=56-male VALIDATION_STRINGENCY=SILENT

mv 56-F.bam 56-F.bamOLD
mv 56-F2.bam 56-F2.bamOLD
mv 56-M.bam 56-M.bamOLD
mv 56-M2.bam 56-M2.bamOLD

