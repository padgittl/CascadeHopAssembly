###############################
##### createGenomePlot.py #####
###############################

# the purpose of this script is to visualize ordered primary and associate contigs, gene, LTR, SNP, and average heterozygosity values for an entire linkage group. Primary contigs are ordered by linkage map genetic positions

python createGenomePlot.py primaryAssembly_contigMap.txt primaryContigLengths.txt associateContigLengths.txt heterozygosity.csv linkageMap.csv repeatFilteredHopCascadeGeneModels.gff ltrs.gff linkageGroupNumber windowSize


# createGenomePlot.py requires svgwrite

# heterozygosity.csv format -->
Site Number,Site Name,Chromosome,Physical Position,Number of Taxa,Ref,Alt,Major Allele,Major Allele Gametes,Major Allele Proportion,Major Allele Frequency,Minor Allele,Minor Allele Gametes,Minor Allele Proportion,Minor Allele Frequency,Allele 3,Allele 3 Gametes,Allele 3 Proportion,Allele 3 Frequency,Proportion Missing,Number Heterozygous,Proportion Heterozygous,Inbreeding Coefficient,Inbreeding Coefficient Scaled by Missing

# linkageMap.csv information used --> linkage group, contig ID, and linkage position

# primaryContigLengths.txt and associateContigLengths.txt format --> contigID\scontigLength

# primaryAssembly_contigMap.txt (output from deduplicateContigs.py script) format --> primaryContigID, associateContigID, primaryAlignmentStartPos, primaryAlignmentStopPos, source

# windowSize is 10000000 (10 Mb)
