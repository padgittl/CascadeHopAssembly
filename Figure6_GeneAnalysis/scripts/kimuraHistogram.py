#!/usr/bin/python
from Bio import SeqIO
import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def readFilteredContigMap(filteredContigMap):
    filteredContigMapDict = {}
    with open(filteredContigMap,'r') as FCM:
        for line in FCM:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in filteredContigMapDict:
                filteredContigMapDict[hapID] = True
    return(filteredContigMapDict)

def readFileListFile(fileListFile):
    fileList = []
    with open(fileListFile,'r') as F:
        for line in F:
            fileList.append(line.strip())
    return fileList

def computeKimuraDistance(fastaFile,filteredContigMapDict):
    purines = ['A','G']
    pyrimidines = ['C','T']
    sequences = SeqIO.parse(fastaFile,'fasta')
    pair = []
    for record in sequences:
        pair.append(record)
    primaryCDS = pair[0]
    haplotigCDS = pair[1]
    totalLen = 0
    transitions = 0
    transversions = 0
    for i in range(len(primaryCDS)):
        if primaryCDS[i] != '-' and haplotigCDS[i] != '-':
            if primaryCDS[i] != haplotigCDS[i]:
                if primaryCDS[i] in purines and haplotigCDS[i] in pyrimidines:
                    transversions += 1
                elif  primaryCDS[i] in pyrimidines and haplotigCDS[i] in purines:
                    transversions += 1
                else:
                    transitions += 1
            # count all ungaped positions as 
            totalLen += 1
    print fastaFile,transitions,totalLen
    p = float(transitions)/float(totalLen)
    q = float(transversions)/float(totalLen)
    d = -(0.5)*np.log((1-2*p-q)*np.sqrt(1-2*q))    
    return (primaryCDS.id, haplotigCDS.id, p,q,d)

usage = "Usage: " + sys.argv[0] + " <filtered contig map> <exportAlignment file list> \n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

filteredContigMap = sys.argv[1]
exportAlignmentFileList = sys.argv[2]

fullPath = filteredContigMap.strip()
fileName = os.path.basename(fullPath)
baseName,baseExt = os.path.splitext(fileName)

filteredContigMapDict = readFilteredContigMap(filteredContigMap)

fileList = readFileListFile(exportAlignmentFileList)
haplotigDistances = []
homolotigDistances = []
for fastaFile in fileList:
    primaryProteinID,haplotigProteinID,p,q,distance = computeKimuraDistance(fastaFile,filteredContigMapDict)
    haplotigID = haplotigProteinID.split('.')[0]
    getHapID = re.search('(.+)_g',haplotigProteinID)
    hapID = getHapID.group(1)
    if hapID in filteredContigMapDict:
        if len(haplotigID.split('_')) == 3:
            haplotigDistances.append(distance)
        elif len(haplotigID.split('_')) == 2:
            homolotigDistances.append(distance)
        else:
            print "error parsing haplotigID " + haplotigID
            sys.exit()

binWidth = 0.01
allDistances = haplotigDistances + homolotigDistances
bins = np.arange(min(allDistances),max(allDistances),binWidth)
mCounts, mBins, mBars = plt.hist(haplotigDistances, bins = bins, histtype = 'stepfilled', color ='r', normed=True, label ='primary vs haplotigs' )
lCounts, lBins, lBars = plt.hist(homolotigDistances, bins = bins, histtype = 'stepfilled', color ='b', normed=True, alpha= 0.5, label= 'primary vs HPCs')

SMALL = 8
MEDIUM = 10
LARGE = 12

#plt.rc('font', size=LARGE)      
#plt.rc('axes', titlesize=LARGE)  
#plt.rc('axes', labelsize=20)  
#plt.rc('xtick', labelsize=20) 
#plt.rc('ytick', labelsize=20) 
plt.rc('legend', fontsize=LARGE) 
plt.rc('figure', titlesize=LARGE)

plt.xlabel("Kimura distance", size = 16)
plt.xlim(0, 0.6)
plt.ylabel("Density", size = 16)
plt.legend(fancybox=False, framealpha=None, edgecolor='inherit',frameon=False)
plt.savefig('kimura_CDS_' + baseName + '_test.pdf')
plt.savefig('kimura_CDS_' + baseName + '_test.svg')
