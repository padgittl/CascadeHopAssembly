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

def computeUngappedLengthDiff(seq1,seq2):
    seq1NoGaps = seq1.replace('-','')
    seq2NoGaps = seq2.replace('-','')
    diff = len(seq1NoGaps) - len(seq2NoGaps)
    return diff 

def extractCodonList(seq1,seq2):
    codonList = []
    primaryPos = -1
    haplotigPos = -1
    primaryFrame = -1
    haplotigFrame = -1
    primaryCodon = ""
    haplotigCodon = ""
    for i in range(len(seq1)):
        if seq1[i] != "-":
            primaryPos += 1
            primaryFrame = primaryPos % 3
            haplotigFrame = haplotigPos % 3
            primaryCodon += seq1[i]
            haplotigCodon += seq2[i]
        if seq2[i] != "-":
            haplotigPos += 1
        # store data, reset the primary frame:
        if primaryFrame == 2:
            codonList.append((primaryCodon,haplotigCodon))
            primaryFrame = 0
            primaryCodon = ""
            haplotigFrame = 0
            haplotigCodon = ""
    return codonList

def readAlignmentFastaFile(alignmentFastaFile,filteredContigMapDict):
    diffList = []
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
                if haplotigBase in filteredContigMapDict:
                    if primaryExt == haplotigExt:
                        if len(primarySeq) == len(haplotigSeq):
                            diff = computeUngappedLengthDiff(primarySeq,haplotigSeq)
                            codonList = extractCodonList(primarySeq,haplotigSeq)
                            diffList.append(diff)

                        else: 
                            print "error parsing file, lengths don't match!"
                            sys.exit()
                    else:
                        print "error parsing file: ", primaryDef, haplotigDef
                        sys.exit()
            else:
                print "error parsing alignment file: " + alignmentFastaFile
                sys.exit()

    return(diffList)

##########
## MAIN ##
##########

if len(sys.argv) != 3 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <file> <overlapping haplotig file>"
    sys.exit()

# read in the input file 
alignmentFastaListFile = sys.argv[1]
filteredContigMap = sys.argv[2]

filteredContigMapDict = readFilteredContigMap(filteredContigMap)
alignmentFastaList = readAlignmentFastaListFile(alignmentFastaListFile)

fullDiffList = []
for alignmentFastaFile in alignmentFastaList:
    diffList = readAlignmentFastaFile(alignmentFastaFile,filteredContigMapDict)
    fullDiffList.extend(diffList)

binSize = 1
minVal = min(fullDiffList)
maxVal = max(fullDiffList)
bins = np.arange(minVal,maxVal,binSize)
plt.hist(fullDiffList, bins=bins)
plt.xlabel("Difference in CDS length (bp)",size=16)
plt.ylabel("Density",size=16)
plt.xlim(-27,27)
#plt.yscale('symlog')
plt.savefig("cdsLengthDiff_ylog_filteredHaps.pdf")
plt.savefig("cdsLengthDiff_ylog_filteredHaps.svg")
