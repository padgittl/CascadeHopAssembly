##### Prepare data for tree

# Run clustalw2
clustalw2 -infile=hits.fasta -type=PROTEIN -outfile=hits.aln > hits.clustalw2.out

# Convert to phylip-relaxed format
>>> from Bio import AlignIO
>>> alignment = AlignIO.parse("hits.aln","clustal")
>>> AlignIO.write(alignment,open("hits.phy","w"),"phylip-relaxed")

# Run phyml
phyml -i hits.phy -d aa -m Blosum62 -c 4 -a e

#########################
##### createTree.py #####
#########################

python createTree.py phymlTree.txt hopTopHits.txt cannabisTopHits.txt 

# format of hopTopHits.txt and cannabisTopHits.txt --> geneID,uniprotID,uniprotDescription,percIdent,eValue,bitScore,queryCov

########################
##### domainViz.py #####
########################

python domainViz.py gene.domtblout

# gene.domtblout is specific to one gene of interest

