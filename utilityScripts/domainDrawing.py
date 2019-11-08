import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

sys.path.append('/nfs0/BB/Hendrix_Lab/Hops/svgwrite/svgwrite-1.2.1/')
import svgwrite
from svgwrite import cm, mm

def readHMMFile(hmmFile,eValueThresh,scaleFactor):
    hmmCoordDict = {}
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
                            if (hopGeneID,hopGeneLen) not in hmmCoordDict:
                                hmmCoordDict[(hopGeneID,hopGeneLen)] = []
                            hmmCoordDict[(hopGeneID,hopGeneLen)].append((accessionID,accessionDesc,int(domainStart),int(domainStop),int(scaledStart),int(scaledStop),int(scaledHopGeneLen)))
    return(hmmCoordDict)

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

def drawGene(hmmCoordDict):
    for hopGeneID,hopGeneLen in hmmCoordDict:
        hopGeneLen = int(hopGeneLen)
        figWidth = '20'
        figHeight = '10'
        fig = svgwrite.Drawing(filename=hopGeneID + ".svg", size=(figWidth, figHeight), profile='full')

        countArray = [0]*hopGeneLen

        fig.add(fig.rect((1,0), (hopGeneLen,30),fill='#2166ac'))
        
        for accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop,scaledHopGeneLen in hmmCoordDict[(hopGeneID,hopGeneLen)]:
            drawDomains(hopGeneID,countArray,domainStart,domainStop,hopGeneLen,accessionID,accessionDesc,fig)
            

def drawDomains(geneID,countArray,domStart,domStop,geneLen,accID,accDesc,fig):
    alignLen = abs(domStop-domStart)
    maxCount = getMaxCountValue(countArray,domStart,domStop)
    shapes = fig.add(fig.g(id='shapes', fill='darkmagenta'))
    
    fig.add(fig.rect((domStart,(maxCount*(-1))-2.5), (alignLen,35),fill='#b2182b'))

    fig.add(fig.text(accID,insert=(domStart,maxCount-4.5)))
    fig.add(fig.text(accDesc,insert=(domStart,maxCount+45)))

    updatedMaxCount = maxCount + 1
    countArray = getPosCount(countArray,domStart,domStop,updatedMaxCount)
    maxCount = getMaxCountValue(countArray,domStart,domStop)

    fig.save()


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hmmscan domtblout file> \n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

hmmFile = sys.argv[1]

eValueThresh = 1e-5
scaleFactor = 10

hmmCoordDict = readHMMFile(hmmFile,eValueThresh,scaleFactor)

drawGene(hmmCoordDict)
