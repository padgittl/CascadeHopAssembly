# MBH
hop_vs_uniprot_MBH.txt

1. Create file with dN/dS rate ratios and uniprot IDs (filtering out overlapping haplotigs)
get_dNdS_withUniProtIDs_sge/get_dNdS_withUniProtIDs_sge_1_sge.sh
python /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/cdsFromAlignment/macseAnalysis/scripts/getGeneHits_filteredHaps.py ../dNdS_fileList.txt ../hop_vs_uniprot_MBH_withMatthews.txt ../filteredContigMap_v5.txt > dNdS_rateRatioWithUniProtIDs.txt

2. link hits to GO terms
genesWithGOTerms_len150_filteredHaps.txt
genesWithGOTerms_len150_filteredHaps.sh

3. GO term-focused file containing GO term, primary geneID, hap geneID, dnds, uniprotID
goTermCentricGeneList.txt
goTermCentricGeneList.sh

4. create GO category-specific files containing uniprotID, go term, go desc, go term count, primary geneId, hap geneID, dnds rate
genesWithBiologicalProcessesGOTerms_withBestHitsFromMBH_len150_filteredHaps.txt
genesWithBiologicalProcessesGOTerms_withBestHitsFromMBH_len150_filteredHaps.sh

5. Hypergeometric test
# produces below files
genesWithBiologicalProcessesGOTerms_withBestHitsFromMBH_len150_filteredHaps_hypergeometricAnalysis.txt
genesWithCellularComponentsGOTerms_withBestHitsFromMBH_len150_filteredHaps_hypergeometricAnalysis.txt
genesWithMolecularFunctionGOTerms_withBestHitsFromMBH_len150_filteredHaps_hypergeometricAnalysis.txt
hypergeometric.sh
