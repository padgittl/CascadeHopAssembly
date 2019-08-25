############################
##### buildMAFIndex.py #####
############################
python buildMAFIndex.py lastzOutputFile.txt

# lastzOutputFile.txt --> primaryContig_vs_associateContig.txt
# outputFile --> primaryContig_vs_associateContig.index.txt

##########################################
##### getAlignmentSequenceFromMAF.py #####
##########################################

python getAlignmentSequenceFromMAF.py lastzOutputFile.txt outputFileName primaryContigAugustusOutput.gff primaryContig_vs_associateContig.index.txt

# outputFile --> primaryContig_vs_associateContig.aligned.cds.fasta

#########################
##### splitFasta.py #####
#########################
python splitFasta.py seqFileList.txt

# seqFileList.txt --> file list of all primaryContig_vs_associateContig.aligned.cds.fasta files

###############################
##### run exportAlignment #####
###############################

java -jar macse_v2.03.jar -prog exportAlignment -align primary_vs_associateContig.aligned.cds.fasta -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS - -out_NT primary_vs_associateContig.aligned.cds.fasta -out_AA primary_vs_associateContig.aligned.cds.fasta

# outputFile --> exportAlignmentOutFile.fasta

############################################
##### processExportAlignmentOutFile.py #####
############################################

python processExportAlignmentOutFile.py exportAlignmentOutFile.fasta 

