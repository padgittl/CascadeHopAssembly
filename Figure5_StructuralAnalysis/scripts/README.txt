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

