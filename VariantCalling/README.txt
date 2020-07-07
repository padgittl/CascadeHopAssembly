##########################
## filterNNNsFromVCF.py ##
##########################

python filterNNNsFromVCF.py geneModels.gff variants.vcf contigLengths.txt


###########################
## calculateErrorRate.py ##
###########################

for f in *gff; do echo python calculateErrorRate.py $f variants.vcf contigLengths.txt; done > launch.sh
--> output file is contigID.error.txt

ls -1 *error.txt > errorFileList.txt


################################
## calculateTotalErrorRate.py ##
################################

python calculateTotalErrorRate.py errorFileList.txt > totalErrorRate.txt

