###############################
##### createGenomePlot.py #####
###############################

python createGenomePlot.py heterozygosity.csv repeatFilteredHopCascadeGeneModels.gff linkageMap.csv primaryContigLengths.txt primaryAssembly_contigMap.txt associateContigLengths.txt linkageGroupNumber ltrs.gff overlapFilteredPrimaryAssembly_contigMap.txt <'yes' or 'no' for overlap filtering>

# createGenomePlot.py requires svgwrite

# heterozygosity.csv format -->
Site Number,Site Name,Chromosome,Physical Position,Number of Taxa,Ref,Alt,Major Allele,Major Allele Gametes,Major Allele Proportion,Major Allele Frequency,Minor Allele,Minor Allele Gametes,Minor Allele Proportion,Minor Allele Frequency,Allele 3,Allele 3 Gametes,Allele 3 Proportion,Allele 3 Frequency,Proportion Missing,Number Heterozygous,Proportion Heterozygous,Inbreeding Coefficient,Inbreeding Coefficient Scaled by Missing

# linkageMap.csv information used --> linkage group, contig ID, and linkage position

# primaryContigLengths.txt and associateContigLengths.txt format --> contigID, contigLength

# primaryAssembly_contigMap.txt (output from deduplicateContigs.py script) format --> primaryContigID, associateContigID, primaryAlignmentStartPos, primaryAlignmentStopPos, source

# overlapFilteredPrimaryAssembly_contigMap.txt --> same format as primaryAssembly_contigMap.txt
