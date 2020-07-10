#### CALCULATE REPEAT CONTENT ####

#############################
## getRepeatPercentages.py ##
#############################

# the purpose of this script is to get the total base count of each contig in the Cascade assembly that is repeat-associated

for f in *gff; do echo python getRepeatPercentages.py $f contigLengths.txt '>' `basename $f gff`txt; done > getRepeatContent.sh

# GFF file is contig-specific 
# contigLengths.txt contains contig lengths in the format --> contigID\tcontigLength
# output file is named contigID.txt



#### STACKED BAR CHART ####

##########################################################
## getRepeatPercentageFromSingleGFFs4StackedBarChart.py ##
##########################################################

# the purpose of this script is to calculate the total repeat content for different repeat types across the Cascade assembly, relative to the total length of the assembly

for f in *.txt; do echo $f; done > repeatCountFileList.txt
python getRepeatPercentageFromSingleGFFs4StackedBarChart.py repeatCountFileList.txt 3711963939 > hop.stackedBarChartValues.txt

# 3711963939 is the total number of bases in the assembly

##################################################
## getOtherRepeatPercentages4StackedBarChart.py ## 
##################################################

# Arabidopsis TAIR10 RepeatMasker output file obtained from NCBI Assembly
# Maize B73 RepeatMasker output file obtained from --> https://www.maizegdb.org/genome/genome_assembly/Zm-B73-REFERENCE-GRAMENE-4.0

python getOtherRepeatPercentages4StackedBarChart.py arabidopsis.out 119668634 maize.out 2075000000

# 119668634 is the size of the TAIR10 assembly
# 2075000000 is the size of the maize assembly

## this script produces these files --> arabidopsis.stackedBarChartValues.txt // maize.stackedBarChartValues.txt

####################################
##### createStackedBarChart.py #####
####################################

python createStackedBarChart.py hop.stackedBarChartValues.txt 3711963939 arabidopsis.stackedBarChartValues.txt 119668634 maize.stackedBarChartValues.txt 2075000000

# hop assembly size --> 3711963939
# arabidopsis assembly size --> 119668634
# maize assembly size --> 2075000000



#### PIE CHART ####

###################################################
## getRepeatPercentageFromSingleGFFs4PieChart.py ##
###################################################

# the purpose this script is to calculate the total repeat content for different repeat types across the Cascade assembly, relative to the total repeat content for all repeat types

for f in *.txt; do echo $f; done > repeatCountFileList.txt
python getRepeatPercentageFromSingleGFFs4PieChart.py repeatCountFileList.txt > hop.pieChartValues.txt

###########################################
## getOtherRepeatPercentages4PieChart.py ##
###########################################

# the purpose of this script is to calculate the total repeat content for different repeat types in the Arabidopsis and maize assemblies

python getOtherRepeatPercentages4PieChart.py arabidopsis.out maize.out

## this script produces these files --> arabidopsis.pieChartValues.txt // maize.pieChartValues.txt

##########################
##### createPieChart.py ##
##########################

python createPieChart.py hop.pieChartValues.txt	arabidopsis.pieChartValues.txt maize.pieChartValues.txt



#### LTR LENGTH BOX PLOT ####

###################################
##### plotLTRLengthBoxplot.py #####
###################################

python plotLTRLengthBoxplot.py denovoLTRs.fasta

# denovoLTRs.fasta generated from LTRharvest
