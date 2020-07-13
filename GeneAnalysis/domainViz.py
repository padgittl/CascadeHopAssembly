import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
import svgwrite
from svgwrite import cm, mm
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

def readHMMFile(hmmFile,eValueThresh,scaleFactor):
    hmmCoordList = []
    with open(hmmFile,'r') as HMM:
        for line in HMM:
            if not line.startswith('#'):
                if not line.startswith('-'):
                    if not line.isspace():
                        domainInfo = line.strip().split()[0:22]
                        accessionID = domainInfo[1]
                        hopGeneID = domainInfo[3]
                        hopGeneLen = domainInfo[5]
                        eValue = domainInfo[6]
                        score = domainInfo[7]
                        domainStart = domainInfo[17]
                        domainStop = domainInfo[18]
                        alignmentProbability = domainInfo[21]
                        accessionDesc = line.strip().split()[22:]
                        accessionDesc = ' '.join(accessionDesc)
                        eValue = float(eValue)
                        hopGeneLen = int(hopGeneLen)
                        scaledHopGeneLen = float(hopGeneLen) / scaleFactor
                        domainStart = int(domainStart)
                        domainStop = int(domainStop)
                        scaledStart = float(domainStart) / scaleFactor
                        scaledStop = float(domainStop) / scaleFactor
                        if eValue < eValueThresh:
                            hmmCoordList.append((hopGeneID,accessionID,accessionDesc,int(domainStart),int(domainStop),int(hopGeneLen),int(scaledStart),int(scaledStop),int(scaledHopGeneLen)))
                            #print accessionID,accessionDesc
    return(hmmCoordList,hopGeneLen)

def readDomainMutationFile(domainMutationFile):
    mutationDict = {}
    # print priTranscriptID,hapTranscriptID,cStart,cEnd,primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,aaPos,domainStart,domainStop,status,dNdS
    # 000000F.g109.t1 000000F_051.g109.t1 1132 1134 TGC C TGT C 378.0 366 459 withinDomain -0.0
    with open(domainMutationFile,'r') as DMF:
        for line in DMF:
            priTranscriptID,hapTranscriptID,cStart,cEnd,primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,aaPos,domainStart,domainStop,status,dNdS = line.strip().split()
            if priTranscriptID not in mutationDict:
                mutationDict[priTranscriptID] = []
            mutationDict[priTranscriptID].append((hapTranscriptID,int(cStart),int(cEnd),primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,float(aaPos),int(domainStart),int(domainStop),status,dNdS))
    return(mutationDict)

def createCountArray(hopGeneLen,scaleFactor):
    scaledGeneLen = hopGeneLen / scaleFactor
    countArray = [0]*hopGeneLen
    return(countArray)

def getPosCount(countArray,start,stop,updatedMaxCount):
    for i in range(start,stop):
        countArray[i] = updatedMaxCount
    return(countArray)

def getMaxCountValue(countArray,start,stop):
    maxValue = 0
    for i in range(start,stop):
        if maxValue < countArray[i]:
            maxValue = countArray[i]
    return(maxValue)

def drawDomains(countArray,hmmCoordList,hopGeneLen,mutationDict):
    figWidth = '20' # '50in'
    figHeight = '10' # '30in'
    fig = svgwrite.Drawing(filename="domainViz_000829F.g21.t2.svg", size=(figWidth, figHeight), profile='full')
    fig.add(fig.rect((1,0), (hopGeneLen,30),fill='#2166ac'))
    for hopGeneID,accessionID,accessionDesc,domainStart,domainStop,hopGeneLen,scaledStart,scaledStop,scaledHopGeneLen in hmmCoordList:
        alignLen = abs(domainStop-domainStart)
        maxCount = getMaxCountValue(countArray,domainStart,domainStop)
        shapes = fig.add(fig.g(id='shapes', fill='darkmagenta'))
        for hapTranscriptID,cStart,cEnd,primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,aaPos,domainStartPos,domainStopPos,status,dNdS in mutationDict[hopGeneID]:
            if 'Mediator' in accessionDesc:
                fig.add(fig.rect((domainStart,-2.5), (alignLen,35),fill='darkmagenta'))
            if hopGeneID in mutationDict:
                for hapTranscriptID,cStart,cEnd,primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,aaPos,domainStartPos,domainStopPos,status,dNdS in mutationDict[hopGeneID]:
                    fig.add(fig.rect((aaPos,-2.5), (1,35),fill='gold'))

        updatedMaxCount = maxCount + 1
        countArray = getPosCount(countArray,domainStart,domainStop,updatedMaxCount)
        maxCount = getMaxCountValue(countArray,domainStart,domainStop)
    fig.save()


domainInfo = mlines.Line2D([], [], color='darkmagenta', label='Mediator complex subunit 25 von Willebrand factor type A', linewidth=6)
mutationInfo = mlines.Line2D([], [], color='gold', label='Mutation', linewidth=6)

plt.legend(handles=[domainInfo,mutationInfo], frameon=False)
plt.axis('off')
plt.savefig('000829F.g21.t2_legend.svg')
plt.savefig('000829F.g21.t2_legend.pdf')


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hmmscan domtblout file> <domain mutation file>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

hmmFile = sys.argv[1]
domainMutationFile = sys.argv[2]

eValueThresh = 1e-5
scaleFactor = 10

mutationDict = readDomainMutationFile(domainMutationFile)
hmmCoordList,hopGeneLen = readHMMFile(hmmFile,eValueThresh,scaleFactor)
countArray = createCountArray(hopGeneLen,scaleFactor)

drawDomains(countArray,hmmCoordList,hopGeneLen,mutationDict)

