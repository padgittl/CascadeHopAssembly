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
        
def readContigPairsFile(contigPairsFile,contigLengths):
    contigDict = {}
    for line in open(contigPairsFile,"r"):
        ID1,ID2 = line.strip().split(" ")
        IDs = sorted([ID1,ID2], key=lambda x:contigLengths[x], reverse=True)
        if (IDs[0],IDs[1]) not in contigDict:
            contigDict[(IDs[0],IDs[1])] = True

            print "mummer -maxmatch -b -c -l 100 path2FastaFiles/" + IDs[0] + ".fa path2FastaFiles/" + IDs[1] + ".fa > " + IDs[0] + "_vs_" + IDs[1] + ".mums\" -r run_" + IDs[0] + "_vs_" + IDs[1]


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <contig pairs file> <contig lengths file>"
if len(sys.argv) != 3:
    print usage
    sys.exit()

contigPairsFile = sys.argv[1]
contigLengthsFile = sys.argv[2]

contigLengths = getContigLengths(contigLengthsFile)
readContigPairsFile(contigPairsFile,contigLengths)
