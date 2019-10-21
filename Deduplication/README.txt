#########################################
##### Run megablast on every contig #####
#########################################

for f in *.fasta; do echo "perl runMegablast.pl $f draftPrimaryAssembly.fasta"; done > megablastJobs.sh

	##### example megablast command -->
	##### megablast -i $query -d $db -o megablastOutput/outfmt8.megablast.txt -e 1e-5 -m 8 -F F -W 28

####################################################
##### printLowDuplicationFrequencyBlastHits.py #####
####################################################

python printLowDuplicationFrequencyBlastHits.py outfmt8.megablast.txt > lowDuplicationFrequencyBlastHits.txt

	##### concatenate lowDuplicationFrequencyBlastHits.txt files for every contig -->
	##### allLowDuplicationFrequencyBlastHits.txt

###############################################################
##### Run LASTZ on candidate primary contig and HPC pairs #####
###############################################################

python prepareLASTZCommands.py primaryContigLengths.txt allLowDuplicationFrequencyBlastHits.txt > launchLASTZ.sh

	##### example LASTZ command -->
	##### lastz path2FastaFiles/candidateHPC.fasta path2FastaFiles/primaryContig.fasta --gfextend --hspthresh=20000 --chain --gapped --output=candidateHPC_vs_primaryContig_lastzPlusStrand.txt --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,continuity,coverage --inner=10000 --identity=80 --strand=plus
	##### create file list --> for f in *lastz*Strand.txt; do echo $f; done > lastzOutputFileList.txt

################################################
##### computeNonOverlappingScoreDensity.py #####
################################################

python computeNonOverlappingScoreDensity.py lastzOutputFileList.txt > results.txt

################################################################
##### Run mummer on candidate primary contig and HPC pairs #####
################################################################

python runMummer.py contigPairs.txt primaryContigLengths.txt > launchMummer.sh

	##### contigPairs.txt is a tab-delimited file containing pairs of IDs -->
	##### shorterContigID	longerContigID

####################################################################
##### Run mummerplot on candidate primary contig and HPC pairs #####
####################################################################

##### example command -->
##### mummerplot -postscript -p longerContigID_vs_shorterContigID longerContigID_vs_shorterContigID.mums

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
