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
                            hmmCoordList.append((accessionID,accessionDesc,int(domainStart),int(domainStop),int(hopGeneLen),int(scaledStart),int(scaledStop),int(scaledHopGeneLen)))
                            print accessionID,accessionDesc
    return(hmmCoordList,hopGeneLen)

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

def drawDomains(countArray,hmmCoordList,hopGeneLen):
    figWidth = '50in'
    figHeight = '30in'
    fig = svgwrite.Drawing(filename="domainViz_CBDAS.svg", size=(figWidth, figHeight), profile='full')
    fig.add(fig.rect((1,0), (hopGeneLen,30),fill='#2166ac'))
    for accessionID,accessionDesc,domainStart,domainStop,hopGeneLen,scaledStart,scaledStop,scaledHopGeneLen in hmmCoordList:
        alignLen = abs(domainStop-domainStart)
        maxCount = getMaxCountValue(countArray,domainStart,domainStop)

        shapes = fig.add(fig.g(id='shapes', fill='darkmagenta'))

        if 'FAD' in accessionDesc:
            accessionDesc = 'FAD-binding domain'
            fig.add(fig.rect((domainStart,-30), (alignLen,30),fill='#b2182b'))
        elif 'BBE' in accessionDesc or 'Berberine' in accessionDesc:
            accessionDesc = 'BBE'
            fig.add(fig.rect((domainStart,-30), (alignLen,30),fill='#fdae61'))

        fig.add(fig.text(accessionDesc,
                         insert=(domainStart, -70),
                         fill='black', font_size='18px')
            )

        updatedMaxCount = maxCount + 1
        countArray = getPosCount(countArray,domainStart,domainStop,updatedMaxCount)
        maxCount = getMaxCountValue(countArray,domainStart,domainStop)
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

hmmCoordList,hopGeneLen = readHMMFile(hmmFile,eValueThresh,scaleFactor)
countArray = createCountArray(hopGeneLen,scaleFactor)

drawDomains(countArray,hmmCoordList,hopGeneLen)

