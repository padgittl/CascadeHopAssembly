####################################################
##### printLowDuplicationFrequencyBlastHits.py #####
####################################################

python printLowDuplicationFrequencyBlastHits.py outfmt8.megablast.txt > lowDuplicationFrequencyBlastHits.txt

################################################
##### computeNonOverlappingScoreDensity.py #####
################################################

python computeNonOverlappingScoreDensity.py lastzOutputFileList.txt > results.txt

###########################
##### processLASTZ.pl #####
###########################

perl processLASTZ.pl results.txt primaryContigLengths.txt

primaryContigLengths.txt format --> primaryContigID, primaryContigLength

outputFile --> clusters.txt

#################################
##### deduplicateContigs.py #####
#################################

python deduplicateContigs.py clusters.txt primaryContigLengths.txt lastzOutputFileList.txt purge_haplotigs_lastzOutputFileList.txt listOfPrimaryContigs.txt listOfHaplotigs.txt curated.contig_associations.log outBaseName

# curated.contig_associations.log --> output from purge_haplotigs

outputFiles --> 

outBaseName.log.txt
outBaseName.fasta
outBaseName_contigMap.txt

outBaseName.log.txt format --> 
filtered by haplotig/HPC: primaryID, associateID, primaryContigCoverageDensity, primaryContigAlignmentStartPos, primaryContigStopPos

outBaseName_contigMap.txt format --> 
primaryID, associateID, primaryContigAlignmentStartPos, primaryContigStopPos, primaryContigCoverageDensity, source 

# source will be Falcon-unzip, purgehaplotigs, or LASTZ
