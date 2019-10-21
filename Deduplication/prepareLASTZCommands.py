import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

def getContigLengths(contigLengthsFile):
    contigLengths = {}
    for line in open(contigLengthsFile):
        contigID,Len = line.strip().split(" ")
        contigLengths[contigID] = int(Len)
    return contigLengths
        
def readAllSigLowBlastHitsFile(allSigLowBlastHits,contigLengths):
    contigDict = {}
    minCount = 2
    for line in open(allSigLowBlastHits,"r"):
        ID1,ID2,count = line.strip().split(" ")
        ID1 = ID1.replace('|arrow','')
        ID2 = ID2.replace('|arrow','')
        if int(count) >= minCount:
            IDs = sorted([ID1,ID2], key=lambda x:contigLengths[x])
	    if (IDs[0],IDs[1]) not in contigDict:
                # add to contig dict to prevent double-running
                contigDict[(IDs[0],IDs[1])] = True
                # run LASTZ commands:
                print "lastz path2FastaFiles/" + IDs[0] + ".fa path2FastaFiles/" + IDs[1] + ".fa --gfextend --hspthresh=20000 --chain --gapped --output=" + IDs[0] + "_vs_" + IDs[1] + "_lastzPlusStrand.txt --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,continuity,coverage --inner=10000 --identity=80 --strand=plus\" -r run_" + IDs[0] + "_vs_" + IDs[1]
                print "lastz path2FastaFiles/" + IDs[0] + ".fa path2FastaFiles/" + IDs[1] + ".fa --gfextend --hspthresh=20000 --chain --gapped --output=" + IDs[0] + "_vs_" + IDs[1] + "_lastzMinusStrand.txt --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,continuity,coverage --inner=10000 --identity=80 --strand=minus\" -r run_" + IDs[0] + "_vs_" + IDs[1]
    
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <contig lengths file> <allSigLowBlastHits File>"
if len(sys.argv) != 3:
    print usage
    sys.exit()

contigLengthsFile = sys.argv[1]
blastHitsCountFile = sys.argv[2]

contigLengths = getContigLengths(contigLengthsFile)
readAllSigLowBlastHitsFile(blastHitsCountFile,contigLengths)


