import sys, os, re
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

def processContigs(clusterFile,primaryIDs,contigMap,coords,scores,contigLengths,purgeHaplotigsFile):
    homologousPrimaryContigs = []
    haplotigCoords = {}
    LOG = open(outBase + '.log.txt','w')
    clusterCoords = readClusterFileWithCoords(clusterFile,coords)
    # load the coordinates
    for primaryID in contigMap:
        # for each haplotig
        for haplotigID in contigMap[primaryID]:
            hapRefCovDensity,refStart,refEnd = getHapMumsCoverageDensity(primaryID,haplotigID)
            if primaryID not in haplotigCoords:
                haplotigCoords[primaryID] = []
            # store haplotig coordinates
            haplotigCoords[primaryID].append((hapRefCovDensity,refStart,refEnd,haplotigID,"FALCON-unzip"))

    homologousCoords = {}
    # read all purge-haplotig IDs
    purgeHap = readPurgeHaplotigsFile(purgeHaplotigsFile,contigLengths)
    for primaryID in purgeHap:
        if primaryID not in homologousCoords:
            homologousCoords[primaryID] = []
        for purgeHapHomologID in purgeHap[primaryID]:
            queryCovDensity,refCovDensity,refStart,refEnd = getMumsCoverageDensity(primaryID,purgeHapHomologID,"purgehap")
            homologousPrimaryContigs.append(purgeHapHomologID)
            # store to a seperate coordsList
            homologousCoords[primaryID].append((refCovDensity,refStart,refEnd,purgeHapHomologID,"purgehaplotigs"))
 
    # read cluster file
    for primaryID in clusterCoords:
        # this if-statement was added in the update from 1->2
        if primaryID not in homologousPrimaryContigs:
            # sort by highest coverageDensity
            if primaryID not in homologousCoords:
                homologousCoords[primaryID] = []
            if primaryID not in haplotigCoords:
                haplotigCoords[primaryID] = []
            clusterCoords[primaryID].sort(key=lambda x:x[0], reverse=True)
            for putativeHomolog in clusterCoords[primaryID]:
                coverageDensity,refStart,refEnd,associateID,source = putativeHomolog
                if not overlapsCoordList(haplotigCoords[primaryID],putativeHomolog):
                    if not overlapsCoordList(homologousCoords[primaryID],putativeHomolog):
                        # if doesn't overlap either, then keep it
                        homologousPrimaryContigs.append(associateID)
                        # store to a seperate coordsList
                        homologousCoords[primaryID].append(putativeHomolog)
                    else:
                        LOG.write("%s\t%s\t%s\t%.1f\t%d\t%d\n" % ("filtered by homolog:",primaryID,associateID,coverageDensity,refStart,refEnd))
                else:
                    LOG.write("%s\t%s\t%s\t%.1f\t%d\t%d\n" % ("filtered by haplotig:",primaryID,associateID,coverageDensity,refStart,refEnd))
    # return list of new homologous contigs
    return(homologousPrimaryContigs,haplotigCoords,homologousCoords)

def readClusterFileWithCoords(clusterFile,coords):
    clusterCoords = {}

    with open(clusterFile,'r') as F:
        for line in F:
            ##Associate     960202656 bp
            if '##Associate' in line:
                contigType,totalLength = line.strip().split('\t')
            elif "Associate" in line:
                # adding to existing cluster
                # 000000F<-010708F        Associate (55040bp)
                contigLinks,contigType = line.strip().split('\t')
                contigs = re.split('<-',contigLinks)
                if len(contigs) == 2:
                    ##########
                    # Here is where we read in a pair of contigs
                    ##########
                    primaryID,associateID = contigs
                    start,end,strand,source,score = coords[(associateID,primaryID)]
                    if source != "purgehap":
                        queryCovDensity,refCovDensity,refStart,refEnd= getMumsCoverageDensity(primaryID,associateID,source)
                        if primaryID not in clusterCoords:
                            clusterCoords[primaryID] = []
                        if queryCovDensity >= mummerCovDensLimit:
                            clusterCoords[primaryID].append((refCovDensity,refStart,refEnd,associateID,"LASTZ"))
                else:
                    print "Error parsing clusters:",  line
                    print contigs
                    sys.exit()
    return clusterCoords

def readPurgeHaplotigsFile(purgeHaplotigsFile,contigLengths):
    p0 = re.compile('(\d+F),PRIMARY')
    pHap = re.compile('-> (\d+F),HAPLOTIG')
    purgeHap = {}
    primaryID = ""
    with open(purgeHaplotigsFile,'r') as F:
        for line in F:
            line = line.strip()
            #print line
            if p0.search(line):
                # store the primary IDs, but this isn't used downstream yet
                match = p0.search(line)
                primaryID = match.group(1)
            if pHap.search(line):
                # this should only store the first HAPLOTIG in the line
                # so haplotigs of haplotigs are not included
                match = pHap.search(line)
                thisID = match.group(1)
                #if contigLengths[primaryID] < contigLengths[thisID]:
                #    primaryID,thisID = thisID,primaryID
                if primaryID not in purgeHap:
                    purgeHap[primaryID] = []
                purgeHap[primaryID].append(thisID)
    return purgeHap

def readClusterFile(clusterFile):
    clusterDict = {}
    with open(clusterFile,'r') as F:
        for line in F:
            if '##Associate' in line:
                contigType,totalLength = line.strip().split('\t')
            elif "Cluster" in line:
                # in a new cluster
                clusterID,sizeText = line.strip().split('\t')
                label,clusterSize = sizeText.split('=')
            elif "Associate" in line:
                # adding to existing cluster
                # 000000F<-010708F        Associate (55040bp)
                contigLinks,contigType = line.strip().split('\t')
                contigs = re.split('<-',contigLinks)
                if len(contigs) == 2:
                    primaryID,associateID = contigs
                    clusterDict[(associateID,primaryID)] = True
                    #print "Stored: ", associateID, primaryID
                else:
                    print "Error parsing clusters:",  line
                    print contigs
                    sys.exit()
    return clusterDict

def readLASTZFiles(lastzOutputFileList,clusterDict,coords,scores,label):
    #/nfs0/Hendrix_Lab/Hops/LASTZ/version6/above1_hsp20000/LASTZFiles/000028F_vs_000011F_lastzMinusStrand.txt
    lastzPattern = re.compile('LASTZFiles/(.*?)_vs_(.*?)_lastz')
    with open(lastzOutputFileList,'r') as F:
        for line in F:
            lastzOutputFile = line.strip()
            match = lastzPattern.search(lastzOutputFile)
            refID = match.group(1)
            queryID = match.group(2)
            if (refID,queryID) in clusterDict:
                #print refID, queryID
                # LASTZ alignment for a pair of contigs
                alignmentList = readLastzOutputFile(lastzOutputFile)
                # remove overlapping alignments
                filteredAlignmentList = filterAlignmentList(alignmentList)
                hit = summarizeAlignments(filteredAlignmentList)
                if hit:
                    # ref is the shorter contig
                    refName,refStrand,minRefStart,maxRefEnd,refLength,queryName,queryStrand,minQueryStart,maxQueryEnd,queryLength,scoreSum,scoreDensity,coverageSum,identitySum,continuitySum = hit
                    # now resStart and refStop are in the longer sequence
                    (refStart,refEnd) = (minRefStart,maxRefEnd)
                    (queryStart,queryEnd) = (minQueryStart,maxQueryEnd)
                    if queryStrand == '-': 
                        (queryStart,queryEnd) = (queryLength-maxQueryEnd,queryLength-minQueryStart)
                    if ((refName,queryName)) not in scores:
                        scores[(refName,queryName)] = 0.0
                    if scoreDensity > scores[(refName,queryName)]:
                        # store queryStart,queryEnd because that is the longer (primary) contig.
                        coords[(refName,queryName)] = (queryStart,queryEnd,queryStrand,label,scoreSum/(queryEnd-queryStart+1))
                        scores[(refName,queryName)] = scoreDensity
    return (coords,scores)

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
    # sort by score in descending order
    alignmentList.sort(key=lambda x: x[0], reverse=True)
    return alignmentList

def filterAlignmentList(alignmentList):    
    filteredAlignmentList = []
    for alignment in alignmentList:
        if not checkListOverlap(alignment,filteredAlignmentList):
            if passesFilter(alignment):
                filteredAlignmentList.append(alignment)
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
            #print refName, queryName, oQueryName, minRefStart,maxRefEnd,oMinRefStart,oMaxRefEnd, percentOverlap(minRefStart,maxRefEnd,oMinRefStart,oMaxRefEnd)
            if percentOverlap(minRefStart,maxRefEnd,oMinRefStart,oMaxRefEnd) >= 50:
                return True
        if refName == oQueryName:
            #print "filtering out ", refName, queryName
            if percentOverlap(minRefStart,maxRefEnd,oMinQueryStart,oMaxQueryEnd) >= 50:
                return True
    return False

def overlapsCoordList(coordList,coord):
    # alignment is from alignmentList
    coverageDensity,refStart,refEnd,associateID,source = coord
    for otherCoord in coordList:
        oCoverageDensity,oRefStart,oRefEnd,oAssociateID,source = otherCoord
        if percentOverlap(refStart,refEnd,oRefStart,oRefEnd) >= minQueryOverlap:
            return True
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
        return (refName,refStrand,minRefStart,maxRefEnd,refLength,queryName,queryStrand,minQueryStart,maxQueryEnd,queryLength,scoreSum,scoreDensity,covDensity,idDensity,conDensity)
    else:
        return None

def getContigLengths(contigLengthsFile):
    contigLengths = {}
    for line in open(contigLengthsFile):
        contigID,Len = line.strip().split(" ")
        contigLengths[contigID] = int(Len)
    return contigLengths

#in mums file, forward is first, reverse is second
def readMumsFile(mumsFile):
    startPosFwd = []
    startPosRev = []
    REVERSE = False
    queryContigID = ""
    with open(mumsFile,'r') as F:
        for line in F:
            if line.startswith('>'):
                foundQueryContigID = re.search('(\d+F)',line)
                if foundQueryContigID:
                    queryContigID = foundQueryContigID.group(1)
                if 'Reverse' in line:
                    REVERSE = True
            else:
                refStart,queryStart,mumLen = line.strip().split()
                if REVERSE:
                    startPosRev.append((int(refStart),int(queryStart)))
                else:
                    startPosFwd.append((int(refStart),int(queryStart)))
    return (queryContigID,startPosFwd,startPosRev)

def getDistance(contig2Pos1,contig2Pos2):
    # renamed from getQueryDistance from v3
    return abs(contig2Pos2 - contig2Pos1)

def getClusterDistance(contig1Pos,contig2Pos,cluster):
    minDist1 = distThreshold
    minDist2 = distThreshold
    for oContig1Pos,oContig2Pos in cluster:
        dist1 = getDistance(contig1Pos,oContig1Pos)
        dist2 = getDistance(contig2Pos,oContig2Pos) 
        if dist1 < minDist1:
            minDist1 = dist1
        if dist2 < minDist2:
            minDist2 = dist2
    return (minDist1,minDist2)

def positionPairsNearby(contig1Pos1,contig2Pos1,contig1Pos2,contig2Pos2):
    if abs(contig1Pos2 - contig1Pos1) <= distThreshold and abs(contig2Pos2 - contig2Pos1) <= distThreshold:
        return True
    else:
        return False

def positionPairNearCluster(contig1Pos,contig2Pos,cluster):
    for oContig1Pos,oContig2Pos in cluster:
        if positionPairsNearby(contig1Pos,contig2Pos,oContig1Pos,oContig2Pos):
            return True
    return False 

def clusterPositions(posList):
    clusterList = []
    # check all pairs of positions in posList
    for contig1Pos,contig2Pos in posList:
        NO_CLUSTER = True
        # compare to all existing clusters
        minDist = distThreshold
        for i in range(len(clusterList)):
            cluster = clusterList[i]
            # if new, add to cluster
            if positionPairNearCluster(contig1Pos,contig2Pos,cluster):
                dist1,dist2 = getClusterDistance(contig1Pos,contig2Pos,cluster)
                #if dist1 < minDist and dist2 < minDist:
                if dist2 < minDist:
                    # this has always been dist2 for query distace, but renamed from "dist" in v3.
                    minDist = dist2
                    minClusterID = i
        # at this point, we've compared to all cluster, minDist is specific to this pair of positions.
        if minDist < distThreshold:
            clusterList[minClusterID].append((contig1Pos,contig2Pos))
            # we found a cluster, so NO_CLUSTER becomes False
            NO_CLUSTER = False
        # couldn't find a near cluster, so create new one
        if NO_CLUSTER:
            newCluster = []
            newCluster.append((contig1Pos,contig2Pos))
            clusterList.append(newCluster)
    return clusterList

def filterClusterList(clusterList):
    clusterTupleList = []
    for cluster in clusterList:
        if len(cluster) >= minClusterCount:
            clusterQueryPos = []
            for pos1,pos2 in cluster:
                clusterQueryPos.append(pos2)
            clusterQueryPos.sort()
            clusterStart = clusterQueryPos[0]
            clusterEnd = clusterQueryPos[-1]
            clusterLength = abs(clusterEnd - clusterStart + 1)
            clusterTupleList.append((clusterLength,clusterStart,clusterEnd,cluster))
    # sort by length in descending order
    clusterTupleList.sort(key=lambda x: x[0], reverse=True)
    # filter by overlaps going in descending order
    filteredClusterTupleList = []
    for clusterTuple in clusterTupleList:
        clusterLen,clusterStart,clusterEnd,cluster = clusterTuple
        #print "CLUSTER: ", clusterStart, clusterEnd, clusterLen
        if not clusterOverlapsList(clusterTuple,filteredClusterTupleList):
            filteredClusterTupleList.append(clusterTuple)
            #print "stored ", clusterStart, clusterEnd, clusterLen;
    # now have a temp filtered cluster list
    # addition 4/28/19 sort by position:
    filteredClusterTupleList.sort(key=lambda x: x[1])
    filteredClusterList = []
    for clusterTuple in filteredClusterTupleList:
        clusterLen,clusterStart,clusterEnd,cluster = clusterTuple
        if clusterLen >= 20000:
            filteredClusterList.append(cluster)

    return filteredClusterList
      
def clusterOverlapsList(clusterTuple,clusterTupleList):
    clusterLength,clusterStart,clusterEnd,cluster = clusterTuple
    # loop through current cluster list
    for oClusterLength,oClusterStart,oClusterEnd,oCluster in clusterTupleList:
        if clusterOverlap(clusterStart,clusterEnd,oClusterStart,oClusterEnd) >= minMumClusterOverlap:
            return True
    return False

def checkOverlap(start1,end1,start2,end2):
    if start1 <= start2 and start2 <= end1:
        return True
    if start1 <= end2 and end2 <= end1:
        return True
    return False

def clusterOverlap(start1,end1,start2,end2):
    overlapLength = 0
    shorterLength = min(end1-start1+1,end2-start2+1)
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

def computeCoverageDensity(clusterList,contigLengths,queryContigID,refContigID):
    queryContigLen = contigLengths[queryContigID]
    refContigLen = contigLengths[refContigID]
    ######################
    ## filtered clusters #    
    ######################
    filteredClusterList = filterClusterList(clusterList)
    totalFilteredQueryClusterCoverage = 0
    totalFilteredRefClusterCoverage = 0
    refPositions = []
    for cluster in filteredClusterList:
        queryClusterPositions = []
        refClusterPositions = []
        for pos1,pos2 in cluster:
            #print pos2
            queryClusterPositions.append(pos2)
            refClusterPositions.append(pos1)
            # pos1 should be x-axis, reference sequence
            refPositions.append(pos1)
        queryClusterPositions.sort()
        filteredQueryClusterCoverage = abs(queryClusterPositions[-1] - queryClusterPositions[0])
        filteredRefClusterCoverage = abs(refClusterPositions[-1] - refClusterPositions[0])
        #print 'first filtered cluster start pos,second filtered cluster start pos'
        #print queryPositions[0],queryPositions[-1]
        totalFilteredQueryClusterCoverage += filteredQueryClusterCoverage
        totalFilteredRefClusterCoverage += filteredRefClusterCoverage
    totalFilteredQueryClusterCoverageDensity = 100.0*totalFilteredQueryClusterCoverage / float(queryContigLen)
    totalFilteredRefClusterCoverageDensity = 100.0*totalFilteredRefClusterCoverage / float(refContigLen)
    refStart = -1
    refEnd = -1
    if len(refPositions) > 0:
        refStart = min(refPositions)
        refEnd = max(refPositions)    
    return(totalFilteredQueryClusterCoverageDensity,totalFilteredRefClusterCoverageDensity,refStart,refEnd)

def computeHapCoverageDensity(clusterList,contigLengths,queryContigID,refContigID):
    refContigLen = contigLengths[refContigID]
    ######################
    ## filtered clusters # 
    ######################
    filteredClusterList = filterClusterList(clusterList)
    totalFilteredQueryClusterCoverage = 0
    totalFilteredRefClusterCoverage = 0
    refPositions = []
    for cluster in filteredClusterList:
        queryClusterPositions = []
        refClusterPositions = []
        for pos1,pos2 in cluster:
            #print pos2
            queryClusterPositions.append(pos2)
            refClusterPositions.append(pos1)
            # pos1 should be x-axis, reference sequence
            refPositions.append(pos1)
        queryClusterPositions.sort()
        filteredQueryClusterCoverage = abs(queryClusterPositions[-1] - queryClusterPositions[0])
        filteredRefClusterCoverage = abs(refClusterPositions[-1] - refClusterPositions[0])
        #print 'first filtered cluster start pos,second filtered cluster start pos'
        #print queryPositions[0],queryPositions[-1]
        totalFilteredQueryClusterCoverage += filteredQueryClusterCoverage
        totalFilteredRefClusterCoverage += filteredRefClusterCoverage
    totalFilteredRefClusterCoverageDensity = 100.0*totalFilteredRefClusterCoverage / float(refContigLen)
    refStart = -1
    refEnd = -1
    if len(refPositions) > 0:
        refStart = min(refPositions)
        refEnd = max(refPositions)    
    return(totalFilteredRefClusterCoverageDensity,refStart,refEnd)

def printClusters(clusterList):
    COUNT = 1
    for cluster in clusterList:
        print "Cluster ", COUNT
        for pos1,pos2 in cluster:
            print "\t",pos1,pos2
        print "\n"
        COUNT += 1

def plotMummerFigure(queryID,refID):
    path = "/nfs0/Hendrix_Lab/Hops/Genome/Cascade/lastzPrimarySymLinks"
    Q = path + "/" + queryID  + ".fa"; 
    R = path + "/" + refID + ".fa"; 
    prefix = queryID + "_vs_" + refID
    os.system("mummer -maxmatch -b -c -l 50 " + R + " " + Q + " > Mummerplot/" + prefix + ".mums")
    os.system("mummerplot -postscript -p Mummerplot/" + prefix + " Mummerplot/" + prefix + ".mums")

def getMumsCoverageDensity(primaryID,associateID,source):
    # build mums file name
    if source == 'blast':
        mumsFile = "/nfs0/Hendrix_Lab/Hops/LASTZ/version6/above1_hsp20000/Mummerplot/"+primaryID+"_vs_"+associateID+".mums"
    elif source == 'purgehap':
        mumsFile = "/nfs0/Hendrix_Lab/Hops/LASTZ/versionPurgeHaps/Mummerplot/"+primaryID+"_vs_"+associateID+".mums"
    else:
        print "unknown file"
        sys.exit()
        #plotMummerFigure(associateID,primaryID)
        #mumsFile = "Mummerplot/"+primaryID+"_vs_"+associateID+".mums"
    # read the mums file
    queryContigID,startPosFwd,startPosRev = readMumsFile(mumsFile)
    if queryContigID != associateID:
        print queryContigID, " and ", associateID, " don't match"
        sys.exit()
    # cluster positions by minDist, on query contig
    fwdClusters = clusterPositions(startPosFwd)
    revClusters = clusterPositions(startPosRev)
    # forward
    fwdQueryCovDensity,fwdRefCovDensity,fwdRefStart,fwdRefEnd = computeCoverageDensity(fwdClusters,contigLengths,associateID,primaryID)
    # reverse
    revQueryCovDensity,revRefCovDensity,revRefStart,revRefEnd = computeCoverageDensity(revClusters,contigLengths,associateID,primaryID)
    if fwdQueryCovDensity > revQueryCovDensity:
        return(fwdQueryCovDensity,fwdRefCovDensity,fwdRefStart,fwdRefEnd)
    else:
        return(revQueryCovDensity,revRefCovDensity,revRefStart,revRefEnd)

def getHapMumsCoverageDensity(primaryID,haplotigID):
    # build mums file name

    mumsFile = "MummerplotHaplotig/"+primaryID+"_vs_"+haplotigID+".mums"
    # read the mums file
    queryContigID,startPosFwd,startPosRev = readMumsFile(mumsFile)
    # cluster positions by minDist
    fwdClusters = clusterPositions(startPosFwd)
    revClusters = clusterPositions(startPosRev)
    # forward
    fwdRefCovDensity,fwdRefStart,fwdRefEnd = computeHapCoverageDensity(fwdClusters,contigLengths,haplotigID,primaryID)
    # reverse
    revRefCovDensity,revRefStart,revRefEnd = computeHapCoverageDensity(revClusters,contigLengths,haplotigID,primaryID)
    if fwdRefCovDensity > revRefCovDensity:
        return(fwdRefCovDensity,fwdRefStart,fwdRefEnd)
    else:
        return(revRefCovDensity,revRefStart,revRefEnd)

def readContigFileList(contigFileList):
    contigList = []
    with open(contigFileList,'r') as CFL:
        for line in CFL:
            fileName = line.strip()
            contigID = os.path.splitext(os.path.basename(fileName))[0]
            contigList.append(contigID)
    return contigList

def createContigMap(primaryIDs,haplotigIDs):
    contigMap = {}
    # make sure each primaryID is a key
    for primaryID in primaryIDs:
        contigMap[primaryID] = []
    for haplotigID in haplotigIDs:
        primaryID,count = haplotigID.split('_')
        if primaryID in contigMap:
            contigMap[primaryID].append(haplotigID)
    return contigMap
            
###########
# MAIN ####
###########

usage = "Usage: " + sys.argv[0] + " <Contig cluster file> <contig lengths file> <lastz file list> <purge hap lastz file list> <primary file list> <haplotig file list> <purge haps file> <outbase>"
if len(sys.argv) != 9:
    print usage
    sys.exit()

minQueryOverlap = 50 # min overlap of two associates on the primary 
minMumClusterOverlap = 50
minPercentIdentity = 90
minOverlapPercent = 50
minClusterCount = 5
distThreshold = 10000
mummerCovDensLimit = 25

clusterFile = sys.argv[1]
contigLengthsFile = sys.argv[2]
lastzOutputFileList = sys.argv[3]
purgeHapLastzFileList = sys.argv[4]
primaryContigFileList = sys.argv[5]
haplotigContigFileList = sys.argv[6]
purgeHapsFile = sys.argv[7]
outBase = sys.argv[8]

contigLengths = getContigLengths(contigLengthsFile)

clusterDict = readClusterFile(clusterFile)
coords = {}
scores = {}
coords,scores = readLASTZFiles(lastzOutputFileList,clusterDict,coords,scores,"blast")
coords,scores = readLASTZFiles(purgeHapLastzFileList,clusterDict,coords,scores,"purgehap")

primaryIDs = readContigFileList(primaryContigFileList)
haplotigIDs = readContigFileList(haplotigContigFileList)
contigMap = createContigMap(primaryIDs,haplotigIDs)

homologousPrimaryContigs,haplotigCoords,homologousCoords = processContigs(clusterFile,primaryIDs,contigMap,coords,scores,contigLengths,purgeHapsFile)

outputGenomeFile = outBase + ".fasta"
os.system("rm -f " + outputGenomeFile)
with open(primaryContigFileList,'r') as CFL:
    for line in CFL:
        fileName = line.strip()
        primaryID = os.path.splitext(os.path.basename(fileName))[0]
        if primaryID not in homologousPrimaryContigs:
            command = "cat " + fileName + " >> " + outputGenomeFile
            os.system(command)


outputFile = outBase + "_contigMap.txt"
F = open(outputFile,'w')
for primaryID in primaryIDs:
    if primaryID not in homologousPrimaryContigs:
        if primaryID in haplotigCoords:
            for refCovDensity,refStart,refEnd,associateID,source in haplotigCoords[primaryID]:
                F.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (primaryID, associateID, refStart, refEnd, refCovDensity, source))
        if primaryID in homologousCoords:
            for refCovDensity,refStart,refEnd,associateID,source in homologousCoords[primaryID]:
                F.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (primaryID, associateID, refStart, refEnd, refCovDensity, source))
