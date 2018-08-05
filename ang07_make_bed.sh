#!/bin/bash -ve

# make bed file containing contigs to be included in following analyses

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=15:00:00
# current environment settings are used for the job???
#$ -V
##$ -P molecol
##$ -q molecol.q

# make bed file of all contigs:
cat Littorina_scaffolded_PacBio_run2_7_Oct_2016_hard_masked.fasta | awk '$0 ~ ">" {print 0,c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'> Littorina_scaffolded_lengths.bed
sed -i "s/ /\t/g" Littorina_scaffolded_lengths.bed

rm -f allcontigs_sorted_common.bed

while read contig;do
pcregrep "$contig\t" Littorina_scaffolded_lengths.bed >> allcontigs_sorted_common.bed
done < allcontigs_sorted_common.txt


