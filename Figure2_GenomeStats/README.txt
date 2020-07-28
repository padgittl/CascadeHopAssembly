##########################
##### buildFigure.py #####
##########################

python buildFigure.py fullTableFileList.txt

# fullTableFileList.txt format -->  assemblyName, full_table_assemblyName.tsv

# full_table_assemblyName.tsv is output file from BUSCO

##############################################
##### N50_vs_contigNumber_scatterPlot.py #####
##############################################

python N50_vs_contigNumber_scatterPlot.py <draft primary assembly N50 info> <draft haplotig assembly N50 info> < purge_haplotigs N50 info> <final, deduplicated primary assembly N50 info> <associate assembly after deduplication N50 info> <Shinshuwase assembly N50 info> <filtered associate assembly N50 info> <Teamaker assembly N50 info>

# each N50 assembly info file is formatted in the following way --> 
Total number of contigs: XX
Assembly size is: XX
50 Contig Sum is: XX
Partial Contig Sum is: XX, with counter: XX, and contig length XX

