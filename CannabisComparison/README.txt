# Align hop and cannabis gene sequences to UniProt Embryophyta genes

blastp –query hopGeneModels.pep.fasta –db embryophytaUniprot.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out hop_vs_uniprot.txt
blastp –query embryophytaUniprot.fasta –db hopGeneModels.pep.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out uniprot_vs_hop.txt

blastp –query cs10.protein_coding.2905.PEP.fasta –db embryophytaUniprot.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out cannabis_vs_uniprot.txt
blastp –query embryophytaUniprot.fasta –db cs10.protein_coding.2905.PEP.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out uniprot_vs_cannabis.txt

##################################
##### getHopGeneModelSeqs.py #####
##################################

python getHopGeneModelSeqs.py hopGeneModels.pep.fasta hop_vs_uniprot.txt uniprot_vs_hop.txt uniprotEmbryophyta.fasta

# output file --> hopHits.fasta

#######################################
##### getCannabisGeneModelSeqs.py #####
#######################################

python getCannabisGeneModelSeqs.py cannabisGeneModels.pep.fasta cannabis_vs_uniprot.txt uniprot_vs_cannabis.txt uniprotEmbryophyta.fasta

# output file --> cannabisHits.fasta

##### concatenate hit fasta files --> 
cat hopHits.fasta cannabisHits.fasta > allHits.fasta

##### run clustalw2

clustalw2 -infile=allHits.fasta -type=PROTEIN -outfile=allHits.aln > allHits.clustalw2.out

##### run trimal

trimal -in allHits.aln -gappyout > allHits_noGaps.aln

##### convert to phylip format

>>> from Bio import AlignIO
>>> alignment = AlignIO.parse("allHits_noGaps.aln","clustal")
>>> AlignIO.write(alignment,open("allHits_noGaps.phy","w"),"phylip-relaxed")

##### run phyml

phyml -i allHits_noGaps.phy -d aa -m Blosum62 -c 4 -a e

# output files --> allHits_noGaps.phy_phyml_stats.txt // allHits_noGaps.phy_phyml_tree.txt

##### create tree with createTree.py --> details in directory 'Figure4_Synteny/'

################################################################################

##### identify mutual best hits (MBH)

# align hop and cannabis top hits to each other
blastp -query hopHits.fasta -db cannabisHits.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out hop_vs_cannabis.txt
blastp -query cannabisHits.fasta -db hopHits.fasta -evalue 1e-5 -outfmt '6 std qcovs' -out cannabis_vs_hop.txt

#####################
##### getMBH.py #####
#####################

python getMBH.py hop_vs_cannabis.txt cannabis_vs_hop.txt > MBH.txt
