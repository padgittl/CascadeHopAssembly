import sys,re,os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

##############
# SUBROUTINE #
##############

def readFilteredContigMap(filteredContigMap):
    filteredContigMapDict = {}
    with open(filteredContigMap,'r') as FCM:
        for line in FCM:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in filteredContigMapDict:
                filteredContigMapDict[hapID] = True
    return(filteredContigMapDict)

def createHist(haplotigGapLengthsFile,filteredContigMapDict,baseName):
    haplotigGapLengths = []
    homolotigGapLengths = []
    with open(haplotigGapLengthsFile,'r') as HF:
        for line in HF:
            # 000000F_001     1       730     730
            haplotigID,gapLen,gapStart,gapStop = line.strip().split('\t')
            if haplotigID in filteredContigMapDict:
                if '_' in haplotigID:
                    haplotigGapLengths.append(int(gapLen))
                    #print haplotigID
                else:
                    homolotigGapLengths.append(int(gapLen))
                    #print haplotigID

        binWidth = 1.0
        allValues = haplotigGapLengths + homolotigGapLengths

        SMALL = 8
        MEDIUM = 10
        LARGE = 12

        #plt.rc('font', size=LARGE)          # controls default text sizes
        #plt.rc('axes', titlesize=LARGE)     # fontsize of the axes title
        #plt.rc('axes', labelsize=LARGE)    # fontsize of the x and y labels
        #plt.rc('xtick', labelsize=LARGE)    # fontsize of the tick labels
        #plt.rc('ytick', labelsize=LARGE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=LARGE)    # legend fontsize
        plt.rc('figure', titlesize=LARGE)

        bins = np.arange(min(allValues),max(allValues),binWidth)
        #print min(allValues),max(allValues)
        gap1Counts, gap1Bins, gap1Bars = plt.hist(haplotigGapLengths, bins = bins, histtype = 'stepfilled', color ='r', normed=True, label ='Haplotigs')
        gap2Counts, gap2Bins, gap2Bars = plt.hist(homolotigGapLengths, bins = bins, histtype = 'stepfilled', color ='b', normed=True, label= 'HPCs', alpha=0.5)
        plt.legend(loc='upper right')
        #plt.xlim(0,35)
        
        plt.yscale('log')
        plt.xlabel('Gap length',size=20)
        plt.ylabel('Density',size=20)
        plt.legend(fancybox=False, framealpha=None, edgecolor='inherit',frameon=False,fontsize=20)
        plt.savefig("haplotig_vs_homolotig_gapLengthDist_" + baseName + ".svg", format='svg')
        plt.savefig("haplotig_vs_homolotig_gapLengthDist_" + baseName + ".pdf", format='pdf')


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <haplotig gap length file> <filtered contig map>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

haplotigGapLengthsFile = sys.argv[1]
filteredContigMap = sys.argv[2]

fullPath = filteredContigMap.strip()
fileName = os.path.basename(fullPath)
baseName,baseExt = os.path.splitext(fileName)

filteredContigMapDict = readFilteredContigMap(filteredContigMap)
createHist(haplotigGapLengthsFile,filteredContigMapDict,baseName)
