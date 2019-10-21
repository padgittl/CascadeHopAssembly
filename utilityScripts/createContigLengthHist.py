import sys,os,re
from Bio import SeqIO
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def readLengthsFile(lengthsFile):
    lengthList = []
    with open(lengthsFile,'r') as LF:
        for line in LF:
            contigID,contigLength = line.strip().split()
            lengthList.append(int(contigLength))
    return(lengthList)
            
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + "<contig lengths file> \n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

lengthsFile = sys.argv[1]

lengthList = readLengthsFile(lengthsFile)

stdDev = np.std(lengthList)
# print stdDev

binSize = 75000
minVal = min(lengthList)
maxVal = max(lengthList)
bins = np.arange(minVal,maxVal,binSize)
plt.hist(lengthList, bins=bins, color='mediumseagreen')

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

plt.xlabel("Contig lengths",size=16)
plt.ylabel("Count",size=16)

filePath = lengthsFile.strip()
fileName = os.path.basename(filePath)
baseName,fileExt = os.path.splitext(fileName)
# print baseName
 
plt.savefig(baseName + "Hist.pdf")
plt.savefig(baseName + "Hist.svg")
plt.savefig(baseName + "Hist.png", dpi=600)
