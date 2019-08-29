import re,os,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

def readFilteredContigMap(filteredContigMap):
    filteredContigMapDict = {}
    with open(filteredContigMap,'r') as FCM:
        for line in FCM:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in filteredContigMapDict:
                filteredContigMapDict[hapID] = True
    return(filteredContigMapDict)

def readAlignmentFastaListFile(alignmentFastaListFile):
    alignmentFastaList = []
    with open(alignmentFastaListFile,'r') as F:
        for line in F:
            alignmentFastaList.append(line.strip())
    return alignmentFastaList

def readAlignmentFastaFile(alignmentFastaFile,filteredContigMapDict):
    with open(alignmentFastaFile) as F:
        for line in F:
            if line.startswith('>'):
                primaryDef = line.strip().replace('>','')
                primarySeq = F.next().strip()
                haplotigDef = F.next().strip().replace('>','')
                haplotigSeq = F.next().strip()
                primaryBase,primaryExt = primaryDef.split('_')
                getHaplotigInfo = re.search('(.+)_(g.+)',haplotigDef)
                haplotigBase = getHaplotigInfo.group(1)
                haplotigExt = getHaplotigInfo.group(2)
                # 000000F_vs_000000F_001_lastzPlusStrand.fasta
                if haplotigBase in filteredContigMapDict:
                    if '_' not in haplotigBase:
                        if primaryExt == haplotigExt:
                            if len(primarySeq) == len(haplotigSeq):
                                countBases(primarySeq,haplotigSeq)

def countBases(seq1,seq2):
    if(len(seq1) == len(seq2)):
        for i in range(len(seq1)):
            if(seq1[i] != '-') and (seq2[i] != '-'):
                priNucCount[seq1[i]] += 1
                hapNucCount[seq2[i]] += 1
                pairs[seq1[i],seq2[i]] += 1
                

    #print("%s\t%s\t%s\t%s\t" % (pairs['A','A'],pairs['A','C'],pairs['A','G'],pairs['A','T']))
    #print("%s\t%s\t%s\t%s\t" % (pairs['C','A'],pairs['C','C'],pairs['C','G'],pairs['C','T']))
    #print("%s\t%s\t%s\t%s\t" % (pairs['G','A'],pairs['G','C'],pairs['G','G'],pairs['G','T']))
    #print("%s\t%s\t%s\t%s\t" % (pairs['T','A'],pairs['T','C'],pairs['T','G'],pairs['T','T']))


##########
## MAIN ##
##########

if len(sys.argv) != 3 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <alignment file> <filtered contig map file>"
    sys.exit()

# read in the input file 
alignmentFastaListFile = sys.argv[1]
filteredContigMap = sys.argv[2]

filteredContigMapDict =readFilteredContigMap(filteredContigMap)

heatMapList = []

pairs = {}
pairs['A','A'] = 0
pairs['A','C'] = 0
pairs['A','G'] = 0
pairs['A','T'] = 0

pairs['C','A'] = 0
pairs['C','C'] = 0
pairs['C','G'] = 0
pairs['C','T'] = 0

pairs['G','A'] = 0
pairs['G','C'] = 0
pairs['G','G'] = 0
pairs['G','T'] = 0

pairs['T','A'] = 0
pairs['T','C'] = 0
pairs['T','G'] = 0
pairs['T','T'] = 0

hapNucCount = {}
hapNucCount['A'] = 0
hapNucCount['C'] = 0
hapNucCount['G'] = 0
hapNucCount['T'] = 0

priNucCount = {}
priNucCount['A'] = 0
priNucCount['C'] = 0
priNucCount['G'] = 0
priNucCount['T'] = 0


alignmentFastaList = readAlignmentFastaListFile(alignmentFastaListFile)

for alignmentFastaFile in alignmentFastaList:
    readAlignmentFastaFile(alignmentFastaFile,filteredContigMapDict)

totalPairsCount = sum(pairs.values())
totalPrimaryBaseCount = sum(priNucCount.values())
totalHaplotigBaseCount = sum(hapNucCount.values())
#print totalPrimaryBaseCount,totalHaplotigBaseCount
#print sum(priNucCount['A'])

marginalProb = {}
for base in priNucCount:
    marginalProb[base] = float(priNucCount[base])/totalPrimaryBaseCount
    #print marginalProb,marginalProb[base]

jointProb = {}
for i,j in pairs:
    jointProb[i,j] = float(pairs[i,j])/totalPairsCount
    #print jointProb[i,j]
    
conditionalProb = {}
for x,y in jointProb:
    conditionalProb[x,y] = float(jointProb[x,y])/marginalProb[x]
    print x,y,conditionalProb[x,y]


#row = 4
#column = 4
#countMatrix = [[0 for j in range(row)] for i in range(column)]
#bases = ['A','C','G','T']
#for a in range(len(bases)):
#    for b in range(len(bases)):
#        #countMatrix[a][b]
#        countMatrix[a][b] = conditionalProb[bases[a],bases[b]]
#        #print countMatrix
