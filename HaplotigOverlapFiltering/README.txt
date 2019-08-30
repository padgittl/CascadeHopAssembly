####################################
##### assessHaplotigOverlap.py #####
####################################

python assessHaplotigOverlap.py contigMap.txt primaryContigLengths.txt haplotigLengths.txt <int> > overlapFilteredContigMap.txt

# contigMap.txt is an output file of deduplicateContigs.py, located in in directory 'Deduplication/'
# <int> corresponds to the percentage of allowed associate contig overlap. A value of zero was set for this study, signifying that no overlap was allowed

