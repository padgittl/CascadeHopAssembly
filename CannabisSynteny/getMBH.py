import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
import math

###############
# SUBROUTINES #
###############

def getBestHits(blastOutputFile,eValueThreshold):
    bestHits = {}
    bestScore = {}
    bestQueryCov = {}
    for line in open(blastOutputFile,"r"):
        queryID,uniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
        pIdent = float(pIdent)
        bitScore = float(bitScore)
        eValue = float(eValue)
        if eValue < eValueThreshold:
            if queryID in bestHits:
                if bitScore > bestScore[queryID]:
                    bestScore[queryID] = bitScore
                    bestHits[queryID] = uniprotID
            else:
                # store it. 
                bestScore[queryID] = bitScore
                bestHits[queryID] = uniprotID
    return bestHits

def printMutualBestHits(bestHits1,bestHits2):
    for queryID1 in bestHits1:
        subjectID1 = bestHits1[queryID1]
        # if subjectID1 is a key to bestHits2
        if subjectID1 in bestHits2:
            # subjectID2 is the subjectID in bestHits2 where subjectID1 is the query
            subjectID2 = bestHits2[subjectID1]
            #print subjectID1,subjectID2
            if subjectID2 == queryID1:
                # evm.model.09.1318       001039F.g9.t1  
                print("%s\t%s\t" % (subjectID1,queryID1))
                getContigID = re.search('(\d+F)',queryID1)
                contigID = getContigID.group(1)
                #print contigID
                
            
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <blast output 1> <blast output 2>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

blastOutputFile1 = sys.argv[1]
blastOutputFile2 = sys.argv[2]

eValueThreshold = 1e-5

bestHits1 = getBestHits(blastOutputFile1,eValueThreshold)
bestHits2 = getBestHits(blastOutputFile2,eValueThreshold)
printMutualBestHits(bestHits1,bestHits2)






