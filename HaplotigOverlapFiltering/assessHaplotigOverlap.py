import sys, os, re
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

def readContigLengthsFile(contigLengthsFile):
    lenDict = {}
    with open(contigLengthsFile,'r') as LEN:
        for line in LEN:
            contigID,contigLen = line.strip().split()
            if contigID not in lenDict:
                lenDict[contigID] = int(contigLen)
    return(lenDict)

def readContigMapFile(contigMapFile,lenDict,hapLenDict):
    contigMapDict = {}
    with open(contigMapFile,'r') as CMF:
        for line in CMF:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if primaryID in lenDict:
                primaryLen = lenDict[primaryID]
                if hapID in hapLenDict:
                    hapLen = hapLenDict[hapID]
                    hapAlignLen = int(primaryStop) - int(primaryStart)
                    #print hapAlignLen
                    if primaryID not in contigMapDict:
                        contigMapDict[primaryID] = []
                    contigMapDict[primaryID].append((hapID,int(primaryStart),int(primaryStop),source,int(primaryLen),int(hapLen),float(scoreDens),hapAlignLen))
        # edit in line below
        #contigMapDict[primaryID].sort(key=lambda x:x[5], reverse=True)
    return(contigMapDict)

def calcOverlap(hapID,start2,end2,coordList):
    maxOverlapPercent = 0
    for haplotigID,start1,end1,contigSource,primaryLength,hapLength,scoreDens,hapAlignLen in coordList:
        shorterLength = min(end1-start1+1,end2-start2+1)
        overlapLength = 0
        if hapID != haplotigID:
            if start1 <= start2 and start2 <= end1:
                if end1 < end2:
                    # partially overlapping
                    overlapLength = end1 - start2 + 1
                else:
                    overlapLength = end2 - start2 + 1
                    # start2,stop2 totally overlapping
            elif start1 <= end2 and end2 <= end1:
                if start2 < start1:
                    overlapLength = end2 - start1 + 1
                    # partially overlapping
                else:
                    overlapLength = end2 - start2 + 1
                    # start2,stop2 totally overlapping
            elif start2 <= start1 and start1 <= end2:
                if end2 < end1:
                    # partially overlapping 
                    overlapLength = end2 - start1 + 1
                else:
                    overlapLength = end1 - start1 + 1
                    # totally overlapping
            elif start2 <= end1 and end1 <= end2:
                if start1 < start2:
                    # partially overlapping
                    overlapLength = end1 - start2 + 1
                else:
                    overlapLength = end1 - start1 + 1
                    # totally overlapping
                
            else:
                overlapLength = 0.0

            overlapFraction = float(overlapLength)/float(shorterLength)
            overlapPercent = 100*overlapFraction
            if overlapPercent > maxOverlapPercent:
                maxOverlapPercent = overlapPercent
    return(maxOverlapPercent)

def assessOverlap(contigMapDict,hapLenDict,maxOverlapThresh):
    filteredContigMapDict = {}
    maxOverlapThresh = int(maxOverlapThresh)
    outName = open("overlappingHaplotigs_minThresh" + str(maxOverlapThresh) + ".txt",'w')
    # primaryID,hapID,primaryStart,primaryStop,scoreDens,source
    for primaryID in contigMapDict:
        contigMapDict[primaryID].sort(key=lambda x:x[5], reverse=True)
        #contigMapDict[primaryID].sort(key=lambda x:x[7], reverse=True)
        if primaryID not in filteredContigMapDict:
            filteredContigMapDict[primaryID] = []
        if len(contigMapDict[primaryID]) > 1:
            for hapID,primaryStart,primaryStop,source,primaryLen,hapLen,scoreDens,hapAlignLen in contigMapDict[primaryID]:
                overlapPercent = calcOverlap(hapID,primaryStart,primaryStop,filteredContigMapDict[primaryID])
                if overlapPercent <= maxOverlapThresh:
                    filteredContigMapDict[primaryID].append((hapID,primaryStart,primaryStop,source,primaryLen,hapLen,scoreDens,hapAlignLen))
                    #print primaryID,hapID,primaryStart,primaryStop,overlapPercent,source
                    #outName.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (primaryID,hapID,primaryStart,primaryStop,source,primaryLen,hapLen))
                    #print("%s\t%s\t%s\t%s\t%s\t%s\t" % (primaryID,hapID,primaryStart,primaryStop,scoreDens,source))
        #else:
            #print("%s\t%s\t%s\t%s\t%s\t%s\t" % (primaryID,hapID,primaryStart,primaryStop,scoreDens,source))
    return(filteredContigMapDict)

def printUpdatedContigMapFile(contigMapDict,filteredContigMapDict):
    for primaryID in filteredContigMapDict:
        for hapID,primaryStart,primaryStop,source,primaryLen,hapLen,scoreDens,hapAlignLen in filteredContigMapDict[primaryID]:
            print("%s\t%s\t%s\t%s\t%s\t%s\t" % (primaryID,hapID,primaryStart,primaryStop,scoreDens,source))
            #continue
    for primaryContigID in contigMapDict:
        if len(contigMapDict[primaryContigID]) == 1:
            for haplotigID,priStart,priStop,hapSource,priLen,haplotigLen,scoreDensity,hapAlignLen in contigMapDict[primaryContigID]:
                print("%s\t%s\t%s\t%s\t%s\t%s\t" % (primaryContigID,haplotigID,priStart,priStop,scoreDensity,hapSource))
                #continue

###########
# MAIN ####
###########

usage = "Usage: " + sys.argv[0] + " <contig map file> <primary contig lengths file> <haplotig lengths file> <max overlap threshold value>"
if len(sys.argv) != 5:
    print usage
    sys.exit()

contigMapFile = sys.argv[1]
contigLengthsFile = sys.argv[2]
haplotigLengthsFile = sys.argv[3]
maxOverlapThresh = sys.argv[4]

primaryLenDict = readContigLengthsFile(contigLengthsFile)
hapLenDict = readContigLengthsFile(haplotigLengthsFile)

contigMapDict = readContigMapFile(contigMapFile,primaryLenDict,hapLenDict)
filteredContigMapDict = assessOverlap(contigMapDict,hapLenDict,maxOverlapThresh)

printUpdatedContigMapFile(contigMapDict,filteredContigMapDict)
