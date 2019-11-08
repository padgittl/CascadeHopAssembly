##########################################
##### create contig length histogram #####
##########################################

python createContigLengthHist.py primaryContigLengths.txt
	
	##### contigID	contigLength

#################################
##### create domain drawing #####
#################################

python domainDrawing.py contig.domtblout

	##### contig.domtblout is created by running hmmscan
	##### hmmscan --cpu 8 --domtblout contig.domtblout Pfam-A.hmm contig.fasta > contig.err
	
