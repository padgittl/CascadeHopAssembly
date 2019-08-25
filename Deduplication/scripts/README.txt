##### printLowDuplicationFrequencyBlastHits.py #####
python printLowDuplicationFrequencyBlastHits.py outfmt8.megablast.txt > lowDuplicationFrequencyBlastHits.txt

##### computeNonOverlappingScoreDensity.py #####
python computeNonOverlappingScoreDensity.py lastzOutputFileList.txt > results.txt

##### processLASTZ.pl #####
perl processLASTZ.pl results.txt primaryContigLengths.txt

primaryContigLengths.txt format --> primaryContigID\tprimaryContigLength

outputFile --> clusters.txt

##### deduplicateContigs.py #####
python deduplicateContigs.py clusters.txt primaryContigLengths.txt lastzOutputFileList.txt purge_haplotigs_lastzOutputFileList.txt listOfPrimaryContigs.txt listOfHaplotigs.txt curated.contig_associations.log outBaseName

# curated.contig_associations.log --> output from purge_haplotigs
