import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines

###############
# SUBROUTINES #
###############


def readBacteriaHomologyFile(bacteriaHomologyFile):
    bacteriaHomologyDict = {}
    with open(bacteriaHomologyFile,'r') as B:
        for line in B:
            hopGeneID,uniprotID,eValue,bitScore,queryCov,uniprotDescription = line.strip().split('\t')
            bacteriaHomologyDict[hopGeneID] = 1
    return(bacteriaHomologyDict)


def getUniprotIDDescription(fastaFile):
    uniprotDict = {}
    for record in SeqIO.parse(fastaFile,"fasta"):
        getUniprotID = re.search('sp\|(.+)\|',record.id)
        uniprotID = getUniprotID.group(1)
        if uniprotID not in uniprotDict:
            uniprotDict[uniprotID] = record.description
    return(uniprotDict)


def readHopGFF(gffFile1,bacteriaHomologyDict):
    hopCoordDict = {}
    hopTotalGeneCount = 0
    with open(gffFile1,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'transcript':
                    getTranscriptID = re.search('ID=(\d+F\.g\d+\.t\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in bacteriaHomologyDict:
                        if transcriptID not in hopCoordDict:
                            hopCoordDict[transcriptID] = (contigID,int(start),int(end))
                            hopTotalGeneCount += 1
    return(hopCoordDict,hopTotalGeneCount)


def readCannabisGFF(gffFile2):
    canCoordDict = {}
    canTotalGeneCount = 0
    with open(gffFile2,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                if not line.isspace():
                    canChromID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                    if feature == 'mRNA':
                        getGeneID = re.search('ID=(.+);Parent',attribute)
                        geneID = getGeneID.group(1)
                        if geneID not in canCoordDict:
                            canCoordDict[geneID] = (canChromID,int(start),int(end))
                            canTotalGeneCount += 1
    return(canCoordDict,canTotalGeneCount)


def readHitFile(hitFile):
    topHit = {}
    topScore = {}
    with open(hitFile,'r') as HF:
        for line in HF:
            # hopGene2,hUniprotID2,hPIdent2,hEValue2,hBitScore2,hQueryCov2
            geneID,uniprotID,pIden,eValue,bitScore,queryCov = line.strip().split('\t')
            eValue = float(eValue)
            bitScore = float(bitScore)
            if geneID in topHit:
                if bitScore > topScore[geneID]:
                    topScore[geneID] = bitScore
                    topHit[geneID] = (uniprotID,pIden,eValue,bitScore,queryCov)
            else:
                topHit[geneID] = (uniprotID,pIden,eValue,bitScore,queryCov)
                topScore[geneID] = bitScore
    return(topHit)


def countGenes(coordDict,uniprotDict,topHit):
    totalGeneWithUniprotCount = 0
    totalUniqueGeneWithUniprotCountDict = {}
    totalUniqueGeneWithUniprotCount = 0
    uniprotGeneCount = {} 
    countList = []

    for geneID in coordDict:
        primaryID,geneStart,geneStop = coordDict[geneID]
        if geneID in topHit:
            uniprotID,pIden,eValue,bitScore,queryCov = topHit[geneID]
            totalGeneWithUniprotCount += 1
            if uniprotID not in uniprotGeneCount:
                uniprotGeneCount[uniprotID] = 0
            uniprotGeneCount[uniprotID] += 1
            if uniprotID not in totalUniqueGeneWithUniprotCountDict:
                totalUniqueGeneWithUniprotCountDict[uniprotID] = 1
                totalUniqueGeneWithUniprotCount += 1
                #print tID,uniprotID,eValue,bitScore,queryCov
    for uniprotGeneID in uniprotGeneCount:
        countList.append(uniprotGeneCount[uniprotGeneID])
    return(totalGeneWithUniprotCount,totalUniqueGeneWithUniprotCount,totalUniqueGeneWithUniprotCountDict,uniprotGeneCount,countList)


def createHist(hopUniprotGeneCount,hopCountList,canUniprotGeneCount,canCountList):

    binWidth = 1
    allValues = hopCountList + canCountList

    SMALL = 8
    MEDIUM = 10
    LARGE = 12

    legendHop = mlines.Line2D([], [], color='#812581', linestyle='solid',label='Hop',linewidth=6)
    legendCan = mlines.Line2D([], [], color='#27ad81', linestyle='solid',label='Cannabis',linewidth=6, alpha=0.5)
    legendEnrich = mlines.Line2D([], [], color='black', linestyle='solid',label='Enrichment Value',linewidth=2)

    fig, ax1 = plt.subplots(figsize=(5,4))

    plt.rc('axes', titlesize=LARGE)  
    plt.rc('axes', labelsize=LARGE)  
    plt.rc('xtick', labelsize=LARGE) 
    plt.rc('ytick', labelsize=LARGE) 
    plt.rc('legend', fontsize=LARGE) 
    plt.rc('figure', titlesize=LARGE)

    bins = np.append(np.arange(min(allValues),max(allValues),binWidth), [np.inf])

    ax1.spines['left'].set_color('black')
    ax1.spines['right'].set_color('black')

    hopBlastCounts, hopBlastBins, hopBlastBars = ax1.hist(hopCountList, bins = bins, histtype = 'stepfilled', color ='#812581', density=True, label ='Hop', linestyle = 'solid')
    canBlastCounts, canBlastBins, canBlastBars = ax1.hist(canCountList, bins = bins, histtype = 'stepfilled', color ='#27ad81', density=True, label ='Cannabis', linestyle = 'solid',alpha=0.5)
    plt.xlim(1,20)
    plt.xlabel('UniProt Hit Copy Number',size=16)
    plt.ylabel('Density',size=16)

    ax2 = ax1.twinx()

    enrichmentValueList = []
    
    hopCountList = np.asarray(hopCountList)
    canCountList = np.asarray(canCountList)

    for hopCount,canCount in zip(hopBlastCounts,canBlastCounts):
        enrichmentValue = float(hopCount - canCount) / (canCount + 1)
        enrichmentValueList.append(enrichmentValue)
        #print hopCount,canCount,enrichmentValue
    xRangeList = []
    for i in range(len(enrichmentValueList)):
        shifted_i = i + 1.5
        xRangeList.append(shifted_i)

    ax2.plot(xRangeList,enrichmentValueList,color='black')
    ax2.set_ylim(ymin=0, ymax=None)

    ax2.spines['left'].set_color('black')
    ax2.spines['right'].set_color('black')

    plt.legend(handles=[legendHop,legendCan,legendEnrich], framealpha=None, edgecolor='inherit', frameon=False,loc='upper right', fontsize = '16')
    plt.savefig('blastHitCopyNumberHist_zoom_v2.pdf')
    plt.savefig('blastHitCopyNumberHist_zoom_v2.svg')

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hop hit file> <can hit file> <hop GFF> <can GFF> <UniProt fasta> <bacteria homology file>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()

hopHitFile = sys.argv[1]
canHitFile = sys.argv[2]
gffFile1 = sys.argv[3]
gffFile2 = sys.argv[4]
fastaFile = sys.argv[5]
bacteriaHomologyFile = sys.argv[6]

eValueThreshold = 1e-5

bacteriaHomologyDict = readBacteriaHomologyFile(bacteriaHomologyFile)
uniprotDict = getUniprotIDDescription(fastaFile)
hopCoordDict,hopTotalGeneCount = readHopGFF(gffFile1,bacteriaHomologyDict)
canCoordDict,canTotalGeneCount = readCannabisGFF(gffFile2)

hopTopHit = readHitFile(hopHitFile)
canTopHit = readHitFile(canHitFile)

hopTotalGeneWithUniprotCount,hopTotalUniqueGeneWithUniprotCount,hopTotalUniqueGeneWithUniprotCountDict,hopUniprotGeneCount,hopCountList = countGenes(hopCoordDict,uniprotDict,hopTopHit)

canTotalGeneWithUniprotCount,canTotalUniqueGeneWithUniprotCount,canTotalUniqueGeneWithUniprotCountDict,canUniprotGeneCount,canCountList = countGenes(canCoordDict,uniprotDict,canTopHit)

createHist(hopUniprotGeneCount,hopCountList,canUniprotGeneCount,canCountList)


