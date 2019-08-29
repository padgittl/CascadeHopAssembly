#!/local/cluser/bin/python
import sys

def readBlastFile(blastFile):
    blastHits = []
    queryID = None
    with open(blastFile,'r') as F:
        for line in F:
            qSeqID,sSeqID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore = line.strip().split("\t")
            blastHits.append((sSeqID,int(qStart),int(qEnd),float(pIdent),float(eValue),int(float(bitScore))))
            if not queryID:
                queryID = qSeqID
    blastHits = sorted(blastHits, key=lambda x: x[5], reverse=True)
    return(queryID,blastHits)
    
def getQueryLength(queryID,blastHits):
    topSeqID,topQStart,topQEnd,topPIdent,topEValue,topBitScore = blastHits[0]
    if queryID != topSeqID:
        print "queryID is not top hit! qSeqID=" + queryID + " top sSeqID=" + topSeqID
        sys.exit()
    if int(topPIdent) != 100:
        print "percent identity of top hit not 100: " + str(topPIdent)
        sys.exit()
    if topQStart != 1:
        print "topQStart is not 1, topQStart=" + str(topQStart)
    if float(topEValue) != 0.0:
        print "top topEValue is not 0.0, topEValue=" + str(topEValue)
    return topQEnd

def printBlastHitDupFreq(queryID,queryLength,blastHits):
    countArray = [0]*queryLength
    hitCount = {}
    # count distinct hits for each query position
    for hit in blastHits:
        seqID,qStart,qEnd,pIdent,eValue,bitScore = hit
        if seqID != queryID:
            for i in range(qStart,qEnd):
                countArray[i] += 1
    for hit in blastHits:
        seqID,qStart,qEnd,pIdent,eValue,bitScore = hit
        if seqID != queryID:
            avgCount = 0.0
            length = qEnd - qStart + 1
            for i in range(qStart,qEnd):
                avgCount += float(countArray[i])
            avgCount /= float(length)
            #print queryID, seqID, avgCount, bitScore
            if avgCount < float(10.0) and float(eValue) <= eValueThresh:
                #print queryID, seqID, avgCount, eValue
                if seqID not in hitCount:
                    hitCount[seqID] = 0
                hitCount[seqID] += 1
    for seqID in hitCount:
        print queryID,seqID,hitCount[seqID]

############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <contig vs all blast output file>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

blastFile = sys.argv[1]
eValueThresh = 1e-100

queryID,blastHits = readBlastFile(blastFile)
queryLength = getQueryLength(queryID,blastHits)
printBlastHitDupFreq(queryID,queryLength,blastHits)
