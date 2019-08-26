#####################################
##### createCombinedPieChart.py #####
#####################################

python createCombinedPieChart.py primaryAssembly.rm.gff hopAssemblySize arabidopsis.rm.out arabidopsisAssemblySize maize.rm.out maizeAssemblySize

# hopAssemblySize --> 3711963939
# arabidopsisAssemblySize --> 119668634
# maizeAssemblySize --> 2075000000

# Arabidopsis TAIR10 RepeatMasker output file obtained from NCBI Assembly
# Maize B73 RepeatMasker output file obtained from --> https://www.maizegdb.org/genome/genome_assembly/Zm-B73-REFERENCE-GRAMENE-4.0

####################################
##### createStackedBarChart.py #####
####################################

python createStackedBarChart.py primaryAssembly.rm.gff	hopAssemblySize arabidopsis.rm.out arabidopsisAssemblySize maize.rm.out	maizeAssemblySize

###################################
##### plotLTRLengthBoxplot.py #####
###################################

python plotLTRLengthBoxplot.py denovoLTRs.fasta

# denovoLTRs.fasta generated from LTRharvest
