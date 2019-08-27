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
            if primaryID not in filteredContigMapDict:
                filteredContigMapDict[primaryID] = {}
            if hapID not in filteredContigMapDict[primaryID]:
                filteredContigMapDict[primaryID][hapID] = source
    return(filteredContigMapDict)

def read_dNdS(dNdSFile,filteredContigMapDict,bestHits,longestAssociateGenes):
    dNdSDict = {}
    
    with open(dNdSFile,'r') as F:
        for line in F:
            primaryGeneID,hapGeneID,strand,nonSynSubs,nonSynSites,pN,synSubs,synSites,pS,dN,dS,dNdS = line.strip().split('\t')
            # at this step, undefined values for dn/ds are removed
            primaryGeneID = primaryGeneID.replace('_g','.g')
            hapGeneID = hapGeneID.replace('_g','.g')
            
            getPrimaryContigID = re.search('(\d+F)',primaryGeneID)
            primaryContigID = getPrimaryContigID.group(1)
            
            getHaplotigID = re.search('(\d+F_\d+|\d+F)',hapGeneID)
            haplotigID = getHaplotigID.group(1)

            if primaryGeneID in bestHits:
                uniprotID,eValue,bitScore,queryCov = bestHits[primaryGeneID]
                if primaryContigID in filteredContigMapDict:
                    if haplotigID in filteredContigMapDict[primaryContigID]:
                        source = filteredContigMapDict[primaryContigID][haplotigID]
                        if primaryGeneID in longestAssociateGenes:
                            if hapGeneID in longestAssociateGenes[primaryGeneID]:
                                if not dNdS == '--' :
                                    dN = float(dN)
                                    dN = round(dN,5)
                                    dS = float(dS)
                                    dS = round(dS,5)
                                    # this is to remove negative sign
                                    dNdS = dNdS.replace('-','')
                                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (primaryGeneID,hapGeneID,source,strand,uniprotID,eValue,bitScore,queryCov,dN,dS,dNdS))

def getLongestGeneLen(mutationFile):
    associateContigGenes = {}
    longestAssociateGenes = {}

    with open(mutationFile,'r') as F:
        for line in F:
            primaryProteinID,haplotigProteinID,transitions,transversions,totalLen = line.strip().split('\t')
            primaryProteinID = primaryProteinID.replace('_g','.g')
            haplotigProteinID = haplotigProteinID.replace('_g','.g')

            totalLen = int(totalLen)

            haplotigID = haplotigProteinID.split('.')[0]
            getHapID = re.search('(.+)\.g',haplotigProteinID)
            hapID = getHapID.group(1)

            if primaryProteinID not in associateContigGenes:
                associateContigGenes[primaryProteinID] = []
            associateContigGenes[primaryProteinID].append((haplotigProteinID,totalLen))

        for priGeneID in associateContigGenes:
            if len(associateContigGenes[priGeneID]) > 1:
                associateContigGenes[priGeneID].sort(key=lambda x:x[1], reverse=True)
                hapGeneID,longestLen = associateContigGenes[priGeneID][0]
                if priGeneID not in longestAssociateGenes:
                    longestAssociateGenes[priGeneID] = {}
                if hapGeneID not in longestAssociateGenes[priGeneID]:
                    longestAssociateGenes[priGeneID][hapGeneID] = longestLen
            else:
                for haplotigProteinID,ungappedLen in associateContigGenes[priGeneID]:
                    if priGeneID not in longestAssociateGenes:
                        longestAssociateGenes[priGeneID] = {}
                    if haplotigProteinID not in longestAssociateGenes[priGeneID]:
                        longestAssociateGenes[priGeneID][haplotigProteinID] = ungappedLen
    return(longestAssociateGenes)

def readTopHitsFile(topHitsFile,eValueThreshold):
    bestHits = {}
    lowestEValue = {}
    bestScore = {}
    bestQueryCov = {}
    bestPercIden = {}
    with open(topHitsFile,'r') as HIT:
        for line in HIT:
            geneID,uniprotID,pIdent,eValue,bitScore,queryCov = line.strip().split("\t")
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            if eValue < eValueThreshold:
                if geneID in bestHits:
                    if bestScore[geneID] < bitScore:
                        bestScore[geneID] = bitScore
                        bestHits[geneID] = (uniprotID,eValue,bitScore,queryCov)
                else:
                    bestScore[geneID] = bitScore
                    bestHits[geneID] = (uniprotID,eValue,bitScore,queryCov)
    return(bestHits)

usage = "Usage: " + sys.argv[0] + " <filtered contig map> <dnds file> <top hits file> <mutation & lengths files>\n"
if len(sys.argv) != 5:
    print usage
    sys.exit()

filteredContigMap = sys.argv[1]
dNdSFile = sys.argv[2]
topHitsFile = sys.argv[3]
mutationFile = sys.argv[4]

eValueThreshold = 1e-5

fullPath = filteredContigMap.strip()
fileName = os.path.basename(fullPath)
baseName,baseExt = os.path.splitext(fileName)

filteredContigMapDict = readFilteredContigMap(filteredContigMap)
bestHits = readTopHitsFile(topHitsFile,eValueThreshold)
longestAssociateGenes = getLongestGeneLen(mutationFile)

read_dNdS(dNdSFile,filteredContigMapDict,bestHits,longestAssociateGenes)


