#!/bin/python
from __future__ import division
import sys, os
from scipy.stats import hypergeom

# v3.2: From given list of Fbgns:
# Find associated pathways, calculate expected instances, p-value and q-value
# From the above, find highest rank with q-value <= p-value
# Reject null hypothesis for all ranks up to that highest rank; make sure q-values ascend with rank
# Output to file:
# Pathway ID, pathway name, FDR, expected instances
#   list of genes associated with this pathway: gene ID, gene symbol, youngMean, oldMean, logFoldChange

###############
# SUBROUTINES #
###############

#def readTotalGeneFile(totalGeneFile):
#    with open(totalGeneFile,'r') as F:
#        for line in F:

def readGeneToPathwayFile(geneToPathwayFile):
    pathwayDesc = {}      # dictionary of all pathwayID:pathway desc
    geneToPathways = {}   # dictionary of all Fbgns:pathway info
    populationCount = 0   # dictionary of pathways:number of genes annotated with that pathway
    populationTotal = 0   # total number of known pathway mappings to known Fbgns

    sampleCount = {}
    sampleTotal = {}
    below_dNdS_thresh = []
    pathwayToGenes = {}
    for line in open(geneToPathwayFile,'r'):
        # Q949P1  GO:0046345      abscisic acid catabolic process         1       000645F.g81.t1  007393F.g81.t1  0.0847950964688
        uniprotID,goTerm,goDesc,goCount,primaryGeneID,hapGeneID,dNdS = line.strip().split('\t')
        dNdS = float(dNdS)
        populationTotal += 1
        pathwayDesc[goTerm] = goDesc
        if goTerm not in sampleTotal:
            sampleTotal[goTerm] = 0
        sampleTotal[goTerm] += 1
        if dNdS >= dNdS_thresh:
            populationCount += 1
            if goTerm not in sampleCount:
                sampleCount[goTerm] = 0
            sampleCount[goTerm] += 1
            #if sampleCount[goTerm] >= sampleCountThresh:
        if goTerm not in pathwayToGenes:
            pathwayToGenes[goTerm] = []
        pathwayToGenes[goTerm].append((uniprotID,primaryGeneID))

        #else:
        #below_dNdS_thresh.append((goTerm,goDesc,uniprotID,dNdS))
        #print goTerm,goDesc,uniprotID,dNdS
        #return (pathwayDesc,geneToPathways,populationCount,populationTotal)
    return (pathwayDesc,geneToPathways,populationCount,populationTotal,sampleTotal,sampleCount,pathwayToGenes)

# calculate statistics for each pathway
def computePathwaySignificance(populationCount,populationTotal,sampleCount,sampleTotal,pathwayToGenes):
    pathwaysUnsorted = []
    i = 0

    for pathwayID in sampleCount:
        if sampleCount[pathwayID] >= sampleCountThresh:
            N = sampleTotal[pathwayID]      # sample size, n in wikipedia
            k = sampleCount[pathwayID]     # sample successes, k in wikipedia
            M = populationTotal            # population size, N in wikipedia
            #n = populationCount[pathwayID] # population successes, K in wikipedia
            n = populationCount            # population successes, K in wikipedia
            p = hypergeom.sf(k, M, n, N)   # p value
            expK = sampleTotal[pathwayID]*n/M         # expected successes
            i += 1
            pathwaysUnsorted.append((pathwayID,p,k,expK))
        
    # sort the list obtained above by p value:
    pathwaysSorted = sorted(pathwaysUnsorted, key=lambda p : p[1])
        
    # calculate q values from sorted list, going backwards so q values can be adjusted if needed
    qPrev = 1.0
    pathwayStats = []
    for j in reversed(range(len(pathwaysSorted))):
        pathwayID,p,k,expK = pathwaysSorted[j]
        r = j+1
        qTemp = p*i/r
        qValue = min(qTemp,qPrev)
        qPrev = qValue
        pathwayStats.append((pathwayID,k,expK,p,qValue))
        
    # sort pathway list with stats, ascending by q value
    pathwayStats = sorted(pathwayStats, key=lambda q: q[4])

    return (pathwayStats,pathwayToGenes)

########
# MAIN #
########

usage = "Usage: python " + sys.argv[0] + " <gene 2 GO term file>"
if len(sys.argv) != 2:
    print usage
    print "\n"
    sys.exit()

geneToPathwayFile = sys.argv[1]

dNdS_thresh = 1
sampleCountThresh = 2

pathwayDesc,geneToPathways,populationCount,populationTotal,sampleTotal,sampleCount,pathwayToGenes = readGeneToPathwayFile(geneToPathwayFile)
pathwayStats,pathwayToGenes = computePathwaySignificance(populationCount,populationTotal,sampleCount,sampleTotal,pathwayToGenes)

fileName = os.path.basename(geneToPathwayFile)
baseName,fileExt = os.path.splitext(fileName)
#print baseName
outputFile = baseName + "_hypergeometricAnalysis.txt"   # output file name
OUT = open(outputFile,'w')

for item in pathwayStats:
    pathwayID,k,expK,p,qValue = item
    #OUT.write("%s\t%s\tFDR=%.3f\t%d genes, %.2e expected\n" % (pathwayID,pathwayDesc[pathwayID],qValue,k,expK))
    # alt printing style with commas; easier for parsing
    OUT.write("%s,%s,FDR=%.3f,%d genes, %.2e expected\n" % (pathwayID,pathwayDesc[pathwayID],qValue,k,expK))
         
#print sampleTotal
#for goTerm in sampleCount:
#    print goTerm,sampleCount[goTerm]

#for goTerm in pathwayToGenes:
#    print goTerm,pathwayToGenes[goTerm]

