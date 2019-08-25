import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

def readLastzOutputFile(lastzOutputFile):
    alignmentList = []
    for line in open(lastzOutputFile,"r"):
        if not line.startswith('#'):
            score,refName,refStrand,refLength,refStartPos1,refEnd,queryName,queryStrand,queryLength,queryStartPos1,queryEnd,identity,percentIdentity,continuity,percentContinuity,coverage,percentCoverage = line.strip().split("\t")
            parsePercentIdentity = re.search('(\d+\.\d+)\%',percentIdentity)
            parsePercentContinuity = re.search('(\d+\.\d+)\%',percentContinuity)
            parsePercentCoverage = re.search('(\d+\.\d+)\%',percentCoverage)
            if  parsePercentIdentity and parsePercentContinuity and parsePercentCoverage:
                percentIdentityFloat = parsePercentIdentity.group(1)
                percentContinuityFloat = parsePercentContinuity.group(1)
                percentCoverageFloat = parsePercentCoverage.group(1)
                alignmentList.append((int(score),refName,refStrand,int(refLength),int(refStartPos1),int(refEnd),queryName,queryStrand,int(queryLength),int(queryStartPos1),int(queryEnd),identity,float(percentIdentityFloat),continuity,float(percentContinuityFloat),coverage,float(percentCoverageFloat)))
            else:
                print "ERROR: Could not parse percent coverage for "+ refName + " vs " + queryName
                sys.exit()
    # sort by score in desending order
    alignmentList.sort(key=lambda x: x[0], reverse=True)
    return alignmentList

def filterAlignmentList(alignmentList):    
    filteredAlignmentList = []
    for alignment in alignmentList:
        if not checkListOverlap(alignment,filteredAlignmentList):
            if passesFilter(alignment):
                filteredAlignmentList.append(alignment)
                #print alignment
    return filteredAlignmentList
    
def passesFilter(alignment):
    # alignment is from alignmentList
    score,refName,refStrand,refLength,refStart,refEnd,queryName,queryStrand,queryLength,queryStart,queryEnd,identity,percentIdentity,continuity,percentContinuity,coverage,percentCoverage = alignment
    if percentIdentity >= minPercentIdentity:
        return True
    else:
        return False

def checkListOverlap(alignment,filteredAlignmentList):
    # alignment is from alignmentList
    score,refName,refStrand,refLength,refStart,refEnd,queryName,queryStrand,queryLength,queryStart,queryEnd,identity,percentIdentity,continuity,percentContinuity,coverage,percentCoverage = alignment
    for otherAlignment in filteredAlignmentList:
        oScore,oRefName,oRefStrand,oRefLength,oRefStart,oRefEnd,oQueryName,oQueryStrand,oQueryLength,oQueryStart,oQueryEnd,oIdentity,oPercentIden,oContinuity,oPercentContinuity,oCoverage,oPercentCoverage = otherAlignment
        if percentOverlap(refStart,refEnd,oRefStart,oRefEnd) >= 50:
            return True
    return False

def checkTopHitsOverlap(hit,filteredHitList):
    # hit is from hitList
    refName,refStrand,minRefStart,maxRefEnd,refLength,queryName,queryStrand,minQueryStart,maxQueryEnd,queryLength,scoreSum,scoreDensity,coverageSum,identitySum,continuitySum=hit
    for otherHit in filteredHitList:
        oRefName,oRefStrand,oMinRefStart,oMaxRefEnd,oRefLength,oQueryName,oQueryStrand,oMinQueryStart,oMaxQueryEnd,oQueryLength,oScoreSum,oScoreDensity,oCoverageSum,oIdentitySum,oContinuitySum=otherHit
        if refName == oRefName:
            #print "***", refName, queryName, oQueryName, minRefStart,maxRefEnd,oMinRefStart,oMaxRefEnd, percentOverlap(minRefStart,maxRefEnd,oMinRefStart,oMaxRefEnd)
            if percentOverlap(minRefStart,maxRefEnd,oMinRefStart,oMaxRefEnd) >= 50:
                return True
        #if refName == oQueryName:
        #    if percentOverlap(minRefStart,maxRefEnd,oMinQueryStart,oMaxQueryEnd) >= 50:
        #        print "filtering out ", refName, queryName
        #        return True
    return False

def checkOverlap(start1,end1,start2,end2):
    if start1 <= start2 and start2 <= end1:
        return True
    if start1 <= end2 and end2 <= end1:
        return True
    return False

def percentOverlap(start1,end1,start2,end2):
    shorterLength = min(end1-start1+1,end2-start2+1)
    overlapLength = 0
    # case I
    if start1 <= start2 and start2 <= end1:
        if end1 < end2:
            # case A    
            overlapLength = end1 - start2 + 1
        else:
            # case B
            overlapLength = end2 - start2 + 1
    # case II
    if start1 <= end2 and end2 <= end1:
        if start2 < start1:
            # case A
            overlapLength = end2 - start1 + 1
        else:
            # case B
            overlapLength = end2 - start2 + 1
    # case III
    if start2 <= start1 and start1 <= end2:
        if end2 < end1:
            # case A    
            overlapLength = end2 - start1 + 1
        else:
            # case B
            overlapLength = end1 - start1 + 1
    # case IV
    if start2 <= end1 and end1 <= end2:
        if start1 < start2:
            # case A
            overlapLength = end1 - start2 + 1
        else:
            # case B
            overlapLength = end1 - start1 + 1
    
    overlapFraction = float(overlapLength)/float(shorterLength)
    return 100.0*overlapFraction

def summarizeAlignments(alignmentList):
    scoreSum = 0.0
    coverageSum = 0
    identitySum = 0
    continuitySum = 0
    minQueryStart = 10e10
    maxQueryEnd = 0    
    minRefStart = 10e10
    maxRefEnd = 0   
    parseFraction = re.compile('(\d+)\/(\d+)')
    for alignment in alignmentList: 
        score,refName,refStrand,refLength,refStart,refEnd,queryName,queryStrand,queryLength,queryStart,queryEnd,identity,percentIdentity,continuity,percentContinuity,coverage,percentCoverage = alignment
        idNumDen = parseFraction.search(identity)
        conNumDen = parseFraction.search(continuity)
        covNumDen = parseFraction.search(coverage)
        idNum = idNumDen.group(1)
        conNum = conNumDen.group(1)
        covNum = covNumDen.group(1)
        scoreSum += score
        coverageSum += int(covNum)
        identitySum += int(idNum)
        continuitySum += int(conNum)
        #print refName, queryName, queryStrand, percentCoverage, percentIdentity, identitySum
        if queryStart < minQueryStart:
            minQueryStart = queryStart
        if queryEnd > maxQueryEnd:
            maxQueryEnd = queryEnd
        if refStart < minRefStart:
            minRefStart = refStart
        if refEnd > maxRefEnd:
            maxRefEnd = refEnd
    if scoreSum > 0:
        scoreDensity = float(scoreSum)/refLength
        idDensity = 100.0*float(identitySum)/refLength
        conDensity = 100.0*float(continuitySum)/refLength
        covDensity = 100.0*float(coverageSum)/refLength
        #print "->", refName,refStrand,minRefStart,maxRefEnd,refLength,queryName,queryStrand,minQueryStart,maxQueryEnd,queryLength,scoreSum,scoreDensity,covDensity,idDensity,conDensity
        return (refName,refStrand,minRefStart,maxRefEnd,refLength,queryName,queryStrand,minQueryStart,maxQueryEnd,queryLength,scoreSum,scoreDensity,covDensity,idDensity,conDensity)
    else:
        return None

def getSide(start,end,length):
    beginPerc = round(100.0*float(start)/length,1)
    endPerc = round(100.0*float(length-end)/length,1)
    if beginPerc > endPerc:
        return beginPerc
    else:
        return endPerc

############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <LASTZ Output File List>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

minPercentIdentity = 90
minOverlapPercent = 50
lastzOutputFileList = sys.argv[1]

topHits = []
with open(lastzOutputFileList,'r') as F:
    for line in F:
        lastzOutputFile = line.strip()
        # LASTZ alignment for a pair of contigs
        alignmentList = readLastzOutputFile(lastzOutputFile)
        # remove overlapping alignments
        filteredAlignmentList = filterAlignmentList(alignmentList)        
        hit = summarizeAlignments(filteredAlignmentList)
        if hit:
            topHits.append(hit)

# sort by coverage
topHits.sort(key=lambda x: x[12], reverse=True)

# Filter the top hits to not overlap same ref/target
filteredTopHits = []
for hit in topHits:
    #print hit
    if not checkTopHitsOverlap(hit,filteredTopHits):
        filteredTopHits.append(hit)

for topHit in filteredTopHits:
    refName,refStrand,minRefStart,maxRefEnd,refLength,queryName,queryStrand,minQueryStart,maxQueryEnd,queryLength,scoreSum,scoreDensity,coverageSum,identitySum,continuitySum=topHit
    biggerPerc = getSide(minRefStart,maxRefEnd,refLength)
    #querySide,queryOffset = getSide(minQueryStart,maxQueryEnd,queryLength)
    #barcode = refStrand + refSide + str(refOffset) + queryStrand + querySide + str(queryOffset)
    #refName,queryName,minQueryStart,maxQueryEnd,queryStrand,scoreSum,scoreDensity,coverageSum,identitySum = topHit
    print("%s\t%s\t%d\t%d\t%s\t%.2f\t%.2f\t%f\t%f\t%f\t%.1f" % (refName,queryName,minRefStart,maxRefEnd,queryStrand,scoreSum,scoreDensity,coverageSum,identitySum,continuitySum,biggerPerc))
    #print refName, queryName, minRefStart, maxRefEnd, minQueryStart, maxQueryEnd
    
