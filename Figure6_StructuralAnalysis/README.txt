# Get gap lengths first

###################################
##### getGapLengthsFromMAF.py #####
###################################

python getGapLengthsFromMAF.py lastzOutputFile.txt associateContigGapRate.txt	primaryContigGapRate.txt

# lastzOutputFile.txt --> associateContig_vs_primaryContig.lastz.txt

# concatenate all associateContigGapRate.txt files --> allAssociateContigGapRates.txt

###################################
##### createGapLenDistHist.py #####
###################################

python createGapLenDistHist.py allAssociateContigGapRates.txt contigMap.txt

##############################
##### kimuraHistogram.py #####
##############################

python kimuraHistogram.py contigMap.txt	mafOutputFileList.txt

##########################################
##### getScatterDataFromHaplotigs.py #####
##### getScatterDataFromHPCs.py ##########
##########################################

python getScatterDataFromHaplotigs.py extractedCDSFileList.txt  overlapFilteredContigMap.txt > haplotigScatterData.txt
# extractedCDSFileList.txt --> from ExtractCDS pipeline composed of all primaryContig_vs_associateContig.aligned.cds.fasta files

python getScatterDataFromHPCs.py extractedCDSFileList.txt  overlapFilteredContigMap.txt > hpcScatterData.txt

##############################
##### mutationScatter.py #####
##############################

python mutationScatter.py haplotigScatterData.txt hpcScatterData.txt 


