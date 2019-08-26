#!/usr/bin/python
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os, re

def readFilteredContigMap(filteredContigMap):
    filteredContigMapDict = {}
    with open(filteredContigMap,'r') as FCM:
        for line in FCM:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in filteredContigMapDict:
                filteredContigMapDict[hapID] = True
                #print hapID
    return(filteredContigMapDict)

def readMafFileList(mafFileList):
    mafFiles = []
    with open(mafFileList,'r') as M:
        for mafFile in M:
            mafFiles.append(mafFile.strip())
    return mafFiles

def readMafFile(mafFile):
    purines = ['A','G']
    pyrimidines = ['C','T']
    fullPath = mafFile.strip()
    fileName = os.path.basename(fullPath)
    baseName,fileExt = os.path.splitext(fileName)
    getContigIDs = re.search('(\d+F)_vs_(\d+F_\d+|\d+F)',baseName)
    primaryID = getContigIDs.group(1)
    hapID = getContigIDs.group(2)
    
    output = []
    with open(mafFile,'r') as M:
        for primaryLine in M:
            if primaryLine.startswith('s'):
                sP,seqNameP,startP,blockLenP,strandP,sourceLenP,seqP = primaryLine.strip().split()
                seqP = seqP.upper()
                #print seqP
                primaryBlockStart = int(startP)
                primaryBlockStop = int(primaryBlockStart) + int(blockLenP)
                haplotigLine = M.next()
                sH,seqNameH,startH,blockLenH,strandH,sourceLenH,seqH = haplotigLine.strip().split()
                seqH = seqH.upper()
                hapBlockStart = int(startH)
                hapBlockStop = int(hapBlockStart) + int(blockLenH)
                if seqNameH != hapID:
                    print "error: " + seqNameH + "!=" + hapID
                    sys.exit()
                if seqNameP != primaryID:
                    print "error: " + seqNameP + "!=" + primaryID
                    sys.exit()
                # at this point, we are ready to lok at seqP and seqH
                totalLen = 0
                transitions = 0
                transversions = 0
                for i in range(len(seqP)):
                    if seqP[i] != '-' and seqH[i] != '-':
                        if seqP[i] != seqH[i]:
                            if seqP[i] in purines and seqH[i] in pyrimidines:
                                transversions += 1
                            elif seqP[i] in pyrimidines and seqH[i] in purines:
                                transversions += 1
                            else:
                                transitions += 1
                        # indentation here counts all "sites" (not gapped)
                        totalLen += 1
                p = float(transitions)/float(totalLen)
                q = float(transversions)/float(totalLen)
                d = -(0.5)*np.log((1-2*p-q)*np.sqrt(1-2*q))    
                output.append((seqNameP,primaryBlockStart,primaryBlockStop,seqNameH,hapBlockStart,hapBlockStop,p,q,d))
    return output

usage = "Usage: " + sys.argv[0] + " <filtered contig map> <maf file list> \n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

filteredContigMap = sys.argv[1]
mafFileList = sys.argv[2]

fullPath = filteredContigMap.strip()
fileName = os.path.basename(fullPath)
baseName,baseExt = os.path.splitext(fileName)

filteredContigMapDict = readFilteredContigMap(filteredContigMap)

mafFiles = readMafFileList(mafFileList)

haplotigDistances = []
homolotigDistances = []
for mafFile in mafFiles:
    result = readMafFile(mafFile)
    for primaryID,priBlockStart,priBlockStop,haplotigID,hapBlockStart,hapBlockStop,p,q,d in result:
        if haplotigID in filteredContigMapDict:
            #print primaryID,priBlockStart,priBlockStop,haplotigID,hapBlockStart,hapBlockStop,p,q,d
            if "_" in haplotigID:
                haplotigDistances.append(d)
            else:
                homolotigDistances.append(d)

binWidth = 0.01
allRates = homolotigDistances + haplotigDistances
bins = np.arange(min(allRates),max(allRates),binWidth)
#lCounts, lBins, lBars = plt.hist(haplotigDistances, bins = bins, histtype = 'stepfilled', color ='r', normed=True, label= 'Alignment blocks to haplotigs')
#mCounts, mBins, mBars = plt.hist(homolotigDistances, bins = bins, histtype = 'stepfilled', color ='b', normed=True, label ='Alignment blocks to HPCs', alpha=0.5)
lCounts, lBins, lBars = plt.hist(haplotigDistances, bins = bins, histtype = 'stepfilled', color ='r', normed=True, label= 'Haplotigs')
mCounts, mBins, mBars = plt.hist(homolotigDistances, bins = bins, histtype = 'stepfilled', color ='b', normed=True, label ='HPCs', alpha=0.5) 

plt.xlabel("Kimura distance", size=20)
plt.xlim(0, 0.6)
plt.ylabel("Density", size=20)
plt.legend(fancybox=False, framealpha=None, edgecolor='inherit',frameon=False,fontsize=20)
plt.savefig('kimura_blocks_' + baseName + '.pdf')
plt.savefig('kimura_blocks_' + baseName + '.svg')
