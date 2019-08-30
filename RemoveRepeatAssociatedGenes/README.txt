#########################################
##### getRepeatAssociatedDomains.py #####
#########################################

python getRepeatAssociatedDomains.py Pfam-A.hmm > repeatAssociatedDomains.txt

# run hmmscan on contig-specific files containing genes
# hmmscan --cpu 8 --domtblout contig.domtblout Pfam-A.hmm contigSpecificGenes.fasta > contig.err

# cat all *domtblout files
# for f in *.domtblout; do echo cat $f ">>" hopCascade.domtblout; done > catFiles.sh

###########################################
##### identifyRepeatDomains2Remove.py #####
###########################################

python identifyRepeatDomains2Remove.py hopCascade.domtblout repeatAssociatedDomains.txt > hopCascadeRepeatGenes.txt

################################
##### findGenesWithBoth.py #####
################################

python findGenesWithBoth.py hopCascade.domtblout repeatAssociatedDomains.txt hop_vs_uniprot.txt goTermIDs.tab hopCascadeGeneModels.gff > hopCascadeGeneModelsWithBoth.gff

# hop_vs_uniprot.txt --> align hop gene models and UniProt Embryophyta genes with blastp
# goTermIDs.tab --> obtained from: https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:3193

################################
##### removeRepeatGenes.py #####
################################

python removeRepeatGenes.py hopCascadeGeneModels.fasta hopCascadeRepeatGenes.txt

# output file --> hopCascadeGeneModels_repeatDomainsRemoved.fasta

########################
##### createGFF.py #####
########################

python createGFF.py hopCascadeGeneModels.gff hopCascadeGeneModels_repeatDomainsRemoved.fasta > hopCascadeGeneModels_repeatDomainsRemoved.gff


