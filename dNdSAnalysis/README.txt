# code for dNdS obtained from --> https://github.com/adelq/dnds.git

##################################
##### printSelectionRates.py #####
##################################

python printSelectionRates.py processedExportAlignmentFile.fasta

# processedExportAlignmentFile.fasta generated from ExtractCDS pipeline and each file is gene-specific

##########################
##### getGeneHits.py #####
##########################

python getGeneHits.py overlapFilteredContigMap.txt all_dNdS_rateRatios.txt hopTopHits.txt ungappedCDSLenFile.txt > allGenesWithUniProtIDsAndDNDS.txt

	# all_dNdS_rateRatios.txt is created by concatenating all gene-specific dN/dS files
 	# all_dNdS_rateRatios.txt format --> 
	primaryGeneID,haplotigGeneID,strand,nonSynonymousSubstitutions,nonSynonymousSites,pN,synSubstitutions,synSites,pS,dN,dS,dNdS 

	######################
	### hopBestHits.py ###
	######################
	
	python hopBestHits.py hop_vs_uniprot.txt uniprot_vs_hop.txt > hopTopHits.txt

	# hop_vs_uniprot.txt is obtained by aligning hop gene models to uniprot embryophyta genes -->
	blastp -query hopGenes.fasta -db uniprotEmbryophyta.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out hop_vs_uniprot.txt
	blastp -query uniprotEmbryophyta.fasta -db hopGenes.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out uniprot_vs_hop.txt
			
	# hopTopHits.txt format --> primaryGeneID,uniprotID,percentIdentity,eValue,bitScore,queryCoverage

	############################################
	##### printTransitionsTransversions.py #####
	############################################
	
	python printTransitionsTransversions.py overlapFilteredContigMap.txt processedAlignmentSeqFileList.txt repeatFilteredGeneModels.gff
	output files --> hapUngappedCDSLens.txt // hpcUngappedCDSLens.txt
	cat hapUngappedCDSLens.txt hpcUngappedCDSLens.txt > ungappedCDSLenFile.txt

	# ungappedCDSLenFile.txt format --> primaryProteinID,haplotigGeneID,transitions,transversions,totalUngappedCDSLen
	
	# processedAlignmentSeqFileList.txt comes from ExtractCDS pipeline

###########################
##### linkGenes2GO.py #####
###########################

python linkGenes2GO.py allGenesWithUniProtIDsAndDNDS.txt ../uniprotEmbryophytaGOTerms.tab > allGenesWithGOTerms.txt

# uniprotEmbryophytaGOTerms.tab obtained from: https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:3193
# Columns --> Entry   Entry name      Protein names   Gene names      Organism        Gene ontology (GO)      Gene ontology IDs

######################################
##### identifyEnrichedGOTerms.py #####
######################################

python identifyEnrichedGOTerms.py allGenesWithGOTerms.txt > goTermCentricGeneList.txt

#########################
##### getGOCount.py #####
#########################

python getGOCount.py biologicalProcessesGOTermIDs.tab goTermCentricGeneList.txt > geneGOInfo_biologicalProcesses.txt
python getGOCount.py cellularComponentsGOTermIDs.tab goTermCentricGeneList.txt > geneGOInfo_cellularComponents.txt
python getGOCount.py molecularFunctionGOTermIDs.tab goTermCentricGeneList.txt > geneGOInfo_molecularFunction.txt

# goTermCategorySpecificIDs.tab obtained from: https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:3193

#################################
##### hypergeometricTest.py #####
#################################

python hypergeometricTest.py geneGOInfo_biologicalProcesses.txt biologicalProcesses
python hypergeometricTest.py geneGOInfo_cellularComponents.txt cellularComponents
python hypergeometricTest.py geneGOInfo_molecularFunction.txt molecularFunction
