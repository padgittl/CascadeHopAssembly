#!/usr/bin/python
from Bio import SeqIO
import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def readGFF(gffFile):
    coordDict = {}
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'transcript':
                    getTranscriptID = re.search('ID=(\d+F\.g\d+\.t\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in coordDict:
                        coordDict[transcriptID] = (contigID,int(start),int(end))
                        # print transcriptID
    return(coordDict)

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

def computeTransitionsTransversions(fastaFile,filteredContigMapDict,coordDict):
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
    primaryCDS.id = primaryCDS.id.replace('_','.')
    if primaryCDS.id in coordDict:
        # print primaryCDS.id
        for i in range(len(primaryCDS)):
            if primaryCDS[i] != '-' and haplotigCDS[i] != '-':
                # ungapped positions
                totalLen += 1
                if primaryCDS[i] != haplotigCDS[i]:
                    if primaryCDS[i] in purines and haplotigCDS[i] in pyrimidines:
                        transversions += 1
                    elif  primaryCDS[i] in pyrimidines and haplotigCDS[i] in purines:
                        transversions += 1
                    else:
                        transitions += 1
    # print primaryCDS.id,haplotigCDS.id,totalLen
    return(primaryCDS.id,haplotigCDS.id,transitions,transversions,totalLen)

usage = "Usage: " + sys.argv[0] + " <filtered contig map> <exportAlignment file list> <gene gff file>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

filteredContigMap = sys.argv[1]
exportAlignmentFileList = sys.argv[2]
gffFile = sys.argv[3]

fullPath = filteredContigMap.strip()
fileName = os.path.basename(fullPath)
baseName,baseExt = os.path.splitext(fileName)

coordDict = readGFF(gffFile)

filteredContigMapDict = readFilteredContigMap(filteredContigMap)

fileList = readFileListFile(exportAlignmentFileList)

hapDistanceFile = open('hapDistances.txt','w')
hpcDistanceFile = open('hpcDistances.txt','w')

for fastaFile in fileList:
    primaryProteinID,haplotigProteinID,transitions,transversions,totalLen = computeTransitionsTransversions(fastaFile,filteredContigMapDict,coordDict)
    # 000000F_011_g100.t1
    haplotigID = haplotigProteinID.split('.')[0]
    getHapID = re.search('(.+)_g',haplotigProteinID)
    hapID = getHapID.group(1)
    if primaryProteinID in coordDict:
        if hapID in filteredContigMapDict:
            if len(haplotigID.split('_')) == 3:
                hapDistanceFile.write("%s\t%s\t%s\t%s\t%s\n" % (primaryProteinID,haplotigProteinID,transitions,transversions,totalLen))  
            elif len(haplotigID.split('_')) == 2:
                hpcDistanceFile.write("%s\t%s\t%s\t%s\t%s\n" % (primaryProteinID,haplotigProteinID,transitions,transversions,totalLen))
            else:
                print "error parsing haplotigID " + haplotigID
                sys.exit()
            
