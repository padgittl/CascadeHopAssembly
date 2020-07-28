################################
## calculateIndelErrorRate.py ##
################################

#### the purpose of this script is to calculate the error rate for indels

# each gff file is contig specific, e.g. 000000F.gff
# alternatively, can be run with a genome-scale gff file containing genes from all contigs

for f in *gff; do echo python calculateIndelErrorRate.py $f variants.vcf contigLengths.txt; done > launch.sh
--> example output file is contigID.error.txt

ls -1 *error.txt > errorFileList.txt


##############################
## calculateSNPErrorRate.py ##
##############################

#### the purpose of this script is to calculate the error rate for substitutions only

# each gff file is contig specific, e.g. 000000F.gff
# alternatively, can be run with a genome-scale gff file containing genes from all contigs

for f in *gff; do echo python calculateSNPErrorRate.py $f variants.vcf contigLengths.txt; done > launch.sh
--> example output file is contigID.snpError.txt

ls -1 *error.txt > errorFileList.txt


################################
## calculateTotalErrorRate.py ##
################################

#### the purpose of this script is to calculate the genome-wide error rate for substitutions or indels

# the errorFileList.txt can be generated from either calculateSNPErrorRate.py or calculateIndelErrorRate.py

python calculateTotalErrorRate.py errorFileList.txt > totalErrorRate.txt

