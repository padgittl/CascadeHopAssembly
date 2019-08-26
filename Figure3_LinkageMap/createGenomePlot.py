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

###############
# SUBROUTINES #
###############

def readOverlappingHaplotigList(overlappingHaplotigList):
    overlapHaplotigDict = {}
    with open(overlappingHaplotigList,'r') as OHL:
        for line in OHL:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in overlapHaplotigDict:
                overlapHaplotigDict[hapID] = True
    return(overlapHaplotigDict)

def readContigMapFile_withFiltering(contigMapFile,overlapHaplotigDict):
    contigMapDict = {}
    with open(contigMapFile,'r') as CMF:
        for line in CMF:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if primaryStart != -1:
                if hapID in overlapHaplotigDict:
                    if primaryID not in contigMapDict:
                        contigMapDict[primaryID] = []
                    contigMapDict[primaryID].append((hapID,int(primaryStart),int(primaryStop),source))
    return(contigMapDict)

def readContigMapFile_noFiltering(contigMapFile):
    contigMapDict = {}
    with open(contigMapFile,'r') as CMF:
        for line in CMF:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if primaryStart != -1:
                if primaryID not in contigMapDict:
                    contigMapDict[primaryID] = []
                contigMapDict[primaryID].append((hapID,int(primaryStart),int(primaryStop),source))
    return(contigMapDict)

def readContigLengthsFile(contigLengthsFile):
    lenDict = {}
    with open(contigLengthsFile,'r') as LEN:
        for line in LEN:
            contigID,contigLen = line.strip().split()
            if contigID not in lenDict:
                lenDict[contigID] = int(contigLen)
    return(lenDict)

def readLinkageMapFile(linkageMapFile,lenDict,linkageGroupValue):
    linkageMapList = []
    linkageGroupTotalLen = 0
    with open(linkageMapFile,'r') as QTL:
        for line in QTL:
            # Contig_Position,Lkgrp,pos(cM)
            if not line.startswith('Marker'):
                #if not line.startswith('Contig_Position'):
                line = line.strip().split(',')
                contigPos = line[0]
                # 578_221439F
                getShortenedContigID = re.search('(.+)_',contigPos)
                shortenedContigID = getShortenedContigID.group(1)
                contigID = "0" * (6-len(str(shortenedContigID))) + str(shortenedContigID) + "F"
                linkageGroup = line[1]
                pos_cM = line[2]
                #print contigID
                if linkageGroup == linkageGroupValue:
                    # lnkGrp6
                    #if 842 <= float(pos_cM) and float(pos_cM) <= 995:
                    if 840 <= float(pos_cM) and float(pos_cM) <= 857:
                    #if 842 <= float(pos_cM) and float(pos_cM) <= 859:
                    #if 449 <= float(pos_cM) and float(pos_cM) <= 459:
                        # lnkGrp9
                        #if 40 < float(pos_cM) and float(pos_cM) <= 43.73357591:
                        # if 20 < float(pos_cM) and float(pos_cM) < 26:
                        #if 10 < float(pos_cM) and float(pos_cM) < 50:
                        # lnkGrp10
                        #if 20 < float(pos_cM) and float(pos_cM) < 25:
                        #if 41 < float(pos_cM) and float(pos_cM) < 44:
                        #if 48.5 < float(pos_cM) and float(pos_cM) <= 54:
                        #if float(pos_cM) < 4:
                        #if float(pos_cM) < 10:
                        #if float(pos_cM) < 19.890964 and 11.6930311 < float(pos_cM):
                        if contigID in lenDict:
                            #print contigID,pos_cM
                            contigLen = lenDict[contigID]
                            linkageGroupTotalLen += contigLen
                            linkageMapList.append((contigID,linkageGroup,float(pos_cM)))
    return(linkageMapList,linkageGroupTotalLen)

def drawPrimaryContigs(linkageMapList,lenDict,scaleFactor,linkageGroupTotalLen,megabaseMultiplier,linkageGroupValue,overlappingHaplotigThreshold):
    primaryCoordList = []
    primaryContigStart = 0
    primaryContigStop = 0
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    figWidth = '100in'
    figHeight = '30in'
    fig = svgwrite.Drawing(filename="contigs_lnkGrp" + str(linkageGroupValue) + "_" + overlappingHaplotigThreshold + ".svg", size=(figWidth, figHeight), profile='full')

    scaledMegabase = (float(1000000) / scaleFactor)
    megabaseScaleBar = fig.add(fig.rect((0,30), (scaledMegabase,5),fill='black',stroke='black',stroke_width='0.5'))
    #megabaseScaleBar = fig.add(fig.rect((0,-300), (10,5),fill='black',stroke='black',stroke_width='1'))
    for primaryID,linkageGroup,pos_cM in linkageMapList:
        if primaryID in lenDict:
            primaryContigLen = lenDict[primaryID]
            primaryContigStart = primaryContigStop
            primaryContigStop += primaryContigLen
            primaryScaledStart = primaryContigStart / scaleFactor
            primaryScaledStop = primaryContigStop / scaleFactor
            primaryScaledLen = primaryScaledStop - primaryScaledStart
            fig.add(fig.rect((primaryScaledStart,0), (primaryScaledLen,20),fill='#2166ac',stroke='black',stroke_width='0.5'))
            primaryCoordList.append((primaryID,primaryContigStart,primaryContigStop))
            #print primaryID
    fig.save()
    return(primaryCoordList,fig)

def getHapCoords(primaryCoordList,fig,hapLenDict,contigMapDict,scaleFactor):
    hapCoordList = []
    updatedHapStart = 0
    updatedHapStop = 0
    for primaryID,primaryContigStart,primaryContigStop in primaryCoordList:
        if primaryID in contigMapDict:
            for hapID,hapStart,hapStop,source in contigMapDict[primaryID]:
                hapAlignmentLen = hapStop - hapStart
                updatedHapStart = primaryContigStart + hapStart
                updatedHapStop = updatedHapStart + hapAlignmentLen
                scaledHapStart = updatedHapStart / scaleFactor
                scaledHapStop = updatedHapStop / scaleFactor
                scaledHapLen = scaledHapStop - scaledHapStart
                hapCoordList.append((primaryID,hapID,scaledHapStart,scaledHapStop,scaledHapLen,source))
                print primaryID,hapID,scaledHapStart,scaledHapStop,scaledHapLen,source
                #print primaryID
                #print primaryID,hapID
    hapCoordList.sort(key=lambda x:x[4], reverse=True)
    return(hapCoordList,fig)

def createCountArray(linkageGroupTotalLen,scaleFactor):
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    countArray = [0]*scaledLinkageGroupTotalLen
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

def drawHaplotigs(countArray,hapCoordList,fig,megabaseMultiplier):
    # yScalingFactor = 10
    yScalingFactor = -25
    #yScalingFactor = 50
    for primaryID,hapID,scaledHapStart,scaledHapStop,scaledHapLen,source in hapCoordList:
        maxCount = getMaxCountValue(countArray,scaledHapStart,scaledHapStop)
        #print primaryID,hapID,scaledHapStart,scaledHapStop,scaledHapLen,maxCount*yScalingFactor,source, maxCount
        #print countArray
        if source == "FALCON-unzip":
            #fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-25), (scaledHapLen,20),fill='#b2182b',stroke='black',stroke_width='0.5')) 
            fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-30), (scaledHapLen,20),fill='#b2182b',stroke='black',stroke_width='0.5'))
        if source == "LASTZ":
            #fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-25), (scaledHapLen,20),fill='#f4a582',stroke='black',stroke_width='0.5'))
            fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-30), (scaledHapLen,20),fill='#f4a582',stroke='black',stroke_width='0.5')) 
        if source == "purgehaplotigs":
            #fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-25), (scaledHapLen,20),fill='#f4a582',stroke='black',stroke_width='0.5')) 
            fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-30), (scaledHapLen,20),fill='#f4a582',stroke='black',stroke_width='0.5'))
        updatedMaxCount = maxCount + 1
        countArray = getPosCount(countArray,scaledHapStart,scaledHapStop,updatedMaxCount)
        #print countArray
        maxCount = getMaxCountValue(countArray,scaledHapStart,scaledHapStop)
        #print primaryID,hapID,scaledHapStart,scaledHapStop,scaledHapLen,maxCount*yScalingFactor,source,maxCount
    fig.save()
    return(fig)

def readAugustusGFF(augustusGFF):
    geneData = {}
    with open(augustusGFF,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,stop,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'gene':
                    getGeneID = re.search('ID=(.+)',attribute)
                    geneID = getGeneID.group(1)
                    if contigID not in geneData:
                        geneData[contigID] = []
                    geneData[contigID].append((geneID,strand,int(start),int(stop)))
    return(geneData)

def getHopGeneCoords(primaryCoordList,geneData,scaleFactor):
    genePositionList = []
    for primaryID,primaryContigStart,primaryContigStop in primaryCoordList:
        if primaryID in geneData:
            for geneID,strand,geneStart,geneStop in geneData[primaryID]:
                geneStart = int(geneStart)
                geneStop = int(geneStop)
                geneLen = geneStop - geneStart
                updatedGeneStart = primaryContigStart + geneStart
                updatedGeneStop = updatedGeneStart + geneLen
                scaledGeneStart = updatedGeneStart / scaleFactor
                scaledGeneStop = updatedGeneStop / scaleFactor
                scaledGeneLen = scaledGeneStop - scaledGeneStart
                genePositionList.append((primaryID,geneID,strand,updatedGeneStart,updatedGeneStop,scaledGeneStart,scaledGeneStop,scaledGeneLen))
                #print primaryID,scaledGeneStart,scaledGeneStop,scaledGeneLen
    genePositionList.sort(key=lambda x:x[5], reverse=True)
    return(genePositionList)

def computeGeneDensity(genePositionList,linkageGroupTotalLen,windowSize,scaleFactor,megabaseMultiplier):
    geneDensityList = []
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    scaledWindow = windowSize / scaleFactor
    for i in range(0,scaledLinkageGroupTotalLen,scaledWindow):
        geneCount = 0
        #position = float(i + scaledWindow)
        #position = float(i)
        position = float(i) + (scaledWindow / 2)
        for primaryID,geneID,strand,updatedGeneStart,updatedGeneStop,scaledGeneStart,scaledGeneStop,scaledGeneLen in genePositionList:
            if i <= scaledGeneStart and scaledGeneStart <= i + scaledWindow:
                geneCount += 1
        geneDensity = float(geneCount) * megabaseMultiplier
        geneDensityList.append((position,geneDensity))
        #print primaryID,geneDensity
        #print position,geneCount,geneDensity,scaledMegabaseWindow
    return(geneDensityList)

def read_LTR_GFF_file(LTR_GFF_file):
    ltrData = {}
    # 000000F RepeatMasker    LTR/Gypsy       2768    2924    31.9    +       418     003288F:275585..286585_INT   
    with open(LTR_GFF_file,'r') as LTR:
        for line in LTR:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if contigID not in ltrData:
                    ltrData[contigID] = []
                ltrData[contigID].append((feature,int(start),int(end)))
    return(ltrData)

def getLTRCoords(primaryCoordList,ltrData,scaleFactor):
    ltrPositionList = []
    for primaryID,primaryContigStart,primaryContigStop in primaryCoordList:
        if primaryID in ltrData:
            for LTR_type,ltrStart,ltrStop in ltrData[primaryID]:
                ltrStart = int(ltrStart)
                ltrStop = int(ltrStop)
                ltrLen = ltrStop - ltrStart
                updatedLTRStart = primaryContigStart + ltrStart 
                updatedLTRStop = updatedLTRStart + ltrLen
                scaledLTRStart = updatedLTRStart / scaleFactor
                scaledLTRSTop = updatedLTRStop / scaleFactor
                scaledLTRLen = scaledLTRSTop - scaledLTRStart
                ltrPositionList.append((primaryID,LTR_type,scaledLTRStart,scaledLTRSTop,scaledLTRLen))
    ltrPositionList.sort(key=lambda x:x[4], reverse=True)
    return(ltrPositionList)

def get_LTR_density(ltrPositionList,linkageGroupTotalLen,windowSize,scaleFactor,megabaseMultiplier):
    ltrDensityList = []
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    scaledWindow = windowSize / scaleFactor
    megabaseWindow = windowSize * megabaseMultiplier
    scaledMegabaseWindow = (windowSize * megabaseMultiplier) / scaleFactor
    for i in range(0,scaledLinkageGroupTotalLen,scaledWindow):
        ltrCount = 0
        #position = float(i + scaledWindow)
        position = float(i)
        for primaryID,LTR_type,scaledLTRStart,scaledLTRStop,scaledLTRLen in ltrPositionList:
            if i <= scaledLTRStart and scaledLTRStart <= i + scaledWindow:
                ltrCount += 1
        #ltrDensity = float(ltrCount) / windowSize
        ltrDensity = float(ltrCount) * megabaseMultiplier
        ltrDensityList.append((position,ltrDensity))
    return(ltrDensityList)

def createBarAndLineGraph(geneDensityList,ltrDensityList,windowSize,scaleFactor,megabaseMultiplier,linkageGroupTotalLen,linkageGroupValue,overlappingHaplotigThreshold):
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    scaledWindow = float(windowSize) / scaleFactor
    figWidth = 100
    figHeight = 30
    fig, ax1 = plt.subplots(figsize=(4,1))
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = figWidth
    fig_size[1] = figHeight
    plt.rcParams["figure.figsize"] = fig_size

    #LTR_patch = mpatches.Patch(color='#d6604d', label='LTR density')
    #gene_patch = mpatches.Patch(color='black', label='Gene density')

    geneDensityArray = np.array(geneDensityList)
    ltrDensityArray = np.array(ltrDensityList)
    for position,ltrDensity in ltrDensityArray:
        ax1.bar(position,ltrDensity,scaledWindow, color='#d6604d', align='edge')
        ax1.set_xlim(0,scaledLinkageGroupTotalLen)
        ax1.spines['left'].set_color('#d6604d')
        ax1.spines['right'].set_color('black')
        ax1.tick_params(axis='y', colors='#d6604d')
        #handles1, labels1 = ax1.get_legend_handles_labels()
        #ax1.legend(handles1, labels)
    ax2 = ax1.twinx()
    for j in range(len(geneDensityArray)-1):
        position,geneDensity = geneDensityArray[j]
        nextPosition,nextGeneDensity = geneDensityArray[j+1]
        ax2.plot((position,nextPosition),(geneDensity,nextGeneDensity), linewidth=2, color='black', marker='o', markersize=2)
        ax2.set_xlim(0,scaledLinkageGroupTotalLen)
        ax2.spines['left'].set_color('#d6604d')
        ax2.spines['right'].set_color('black')
        ax2.tick_params(axis='y', colors='black')
        #handles2, labels2 = ax2.get_legend_handles_labels()
        #print handles2,labels2
    #plt.box(False)

    #legend_hetValues = mlines.Line2D([], [], color='#4393c3', marker='o',markersize=6,linestyle="None",label='Heterozygosity values')
    #legend_avgHet = mlines.Line2D([], [], color='black',marker='o',markersize=6,label='Average heterozygosity')
    #plt.legend(handles=[legend_avgHet,legend_hetValues], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))

    legend_genes = mlines.Line2D([], [], color='black', marker='o',markersize=6,label='Gene density')
    legend_LTR = mlines.Line2D([], [], color='#d6604d', marker='s',linestyle="None",markersize=6, label='LTR density')
    #legend_LTR = mpatches.Patch(color='#d6604d', label='LTR density')
    plt.legend(handles=[legend_genes,legend_LTR], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    #plt.legend(handles=[LTR_patch, gene_patch], frameon=False, fontsize='6', ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("gene_repeat_density_lnkGrp" + str(linkageGroupValue) + "_" + overlappingHaplotigThreshold + ".svg", bbox_inches = 'tight', pad_inches=0)

def read_SNP_file(snpFile,primaryCoordList):
    snpDataList = []
    with open(snpFile,'r') as SNP:
        for line in SNP:
            if not line.startswith('Site'):
                Site,SiteName,chrom,pos,numTaxa,Ref,Alt,MajorAllele,MajorAlleleGametes,MajorAlleleProportion,MajorAlleleFreq,MinorAllele,MinorAlleleGametes,MinorAlleleProportion,MinorAlleleFreq = line.strip().split(',')[:15]
                contigID = "0" * (6-len(str(chrom))) + str(chrom)
                contigID += "F"
                snpDataList.append((contigID,int(pos),float(MinorAlleleFreq)))
    return(snpDataList)

def get_SNP_position(snpDataList,primaryCoordList):
    snpPositionList = []
    for primaryContigID,primaryContigStart,primaryContigStop in primaryCoordList:
        for primaryID,snpPos,MinorAlleleFreq in snpDataList:
            if primaryContigID == primaryID:
                updatedSNPPosition = primaryContigStart + snpPos
                snpPositionList.append((primaryID,updatedSNPPosition,MinorAlleleFreq))
    return(snpPositionList)

def get_SNP_density(snpPositionList,linkageGroupTotalLen,windowSize,scaleFactor,megabaseMultiplier):
    hetDensityList = []
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    scaledWindow = windowSize / scaleFactor
    megabaseWindow = windowSize * megabaseMultiplier
    scaledMegabaseWindow = (windowSize * megabaseMultiplier) / scaleFactor
    for i in range(0,scaledLinkageGroupTotalLen,scaledWindow):
        snpCount = 0
        heterozygositySum = 0
        position = float(i) + (scaledWindow / 2)
        #position = float(i)
        for primaryID,updatedSNPPosition,minorAlleleFreq in snpPositionList:
            scaledSNPPosition = float(updatedSNPPosition) / scaleFactor
            p = minorAlleleFreq
            heterozygosityValue = 2*p*(1.0-p)
            if i <= scaledSNPPosition and scaledSNPPosition < i + scaledWindow:
                snpCount += 1
                heterozygositySum += heterozygosityValue
        if snpCount > 1:
            # scale to actual window size!
            #snpDensity = float(snpCount) / windowSize
            snpDensity = float(snpCount) * megabaseMultiplier
            heterozygosityDensity = float(heterozygositySum) / snpCount
            hetDensityList.append((position,heterozygosityDensity,snpDensity))
            #print position,heterozygositySum,snpCount,heterozygosityDensity,snpCount,windowSize,snpDensity
    return(hetDensityList)

def drawSNPs(snpPositionList,hetDensityList,fig,scaleFactor,megabaseMultiplier,linkageGroupTotalLen,linkageGroupValue,overlappingHaplotigThreshold):
    figWidth = 100
    figHeight = 30
    fig, ax1 = plt.subplots(figsize=(4,1))
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = figWidth
    fig_size[1] = figHeight
    plt.rcParams["figure.figsize"] = fig_size
    
    primaryPatch = mpatches.Patch(color='#2166ac', label='Primary contigs')
    HPC_patch = mpatches.Patch(color='#f4a582', label='Homologous primary contigs')
    haplotigPatch = mpatches.Patch(color='#b2182b', label='FALCON-Unzip haplotigs')

    #hetPatch = mpatches.Patch(color='#4393c3', label='Heterozygosity values')
    #avgHetPatch = mpatches.Patch(color='black', label='Average heterozygosity')

    scaledWindow = float(windowSize) / scaleFactor
    scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
    hetDensityArray = np.array(hetDensityList)

    for primaryID,updatedSNPPosition,minorAlleleFreq in snpPositionList:
        scaledSNPPosition = float(updatedSNPPosition) / scaleFactor
        p = minorAlleleFreq
        heterozygosityValue = 2*p*(1.0-p)
        axis1 = ax1.scatter(scaledSNPPosition,heterozygosityValue, color='#4393c3', s=2)
        ax1.set_xlim(0,scaledLinkageGroupTotalLen)
        ax1.spines['left'].set_color('#4393c3')
        ax1.spines['right'].set_color('black')
        ax1.tick_params(axis='y', colors='#4393c3')
        #ax1.margins(0,0.6)
    ax2 = ax1.twinx()
    for i in range(len(hetDensityArray)-1):
        position,hetDensity,snpDensity = hetDensityArray[i]
        nextPosition,nextHetDensity,next_SNP_density = hetDensityArray[i+1]
        axis2 = ax2.plot((position,nextPosition),(hetDensity,nextHetDensity), linewidth=2, color='black', marker='o', markersize=2) 
        ax2.set_xlim(0,scaledLinkageGroupTotalLen)
        ax2.spines['left'].set_color('#4393c3')
        ax2.spines['right'].set_color('black')
        ax2.tick_params(axis='y', colors='black')
        #ax2.margins(0,0.6)
    legend_hetValues = mlines.Line2D([], [], color='#4393c3', marker='o',markersize=6,linestyle="None",label='Heterozygosity values')
    #legend_hetValues = mpatches.Patch(color='#4393c3', label='Heterozygosity values')
    legend_avgHet = mlines.Line2D([], [], color='black',marker='o',markersize=6,label='Average heterozygosity')
    # fontsize='6'
    plt.legend(handles=[legend_avgHet,legend_hetValues], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("snp_density_lnkGrp" + str(linkageGroupValue) + "_" + overlappingHaplotigThreshold + ".svg", bbox_inches = 'tight', pad_inches=0)

    plt.legend(handles=[primaryPatch,HPC_patch,haplotigPatch], frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("contigLegend.svg", bbox_inches = 'tight', pad_inches=0)


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <SNP file> <augustus gene gff file> <map file> <contig lengths file> <contig map file> <haplotig lengths file> <linkage group number> <LTR GFF file> <overlapping haplotig list> <filter overlaps Boolean: yes or no"
if len(sys.argv) != 11:
    print usage
    sys.exit()

snpFile = sys.argv[1]
augustusGFF = sys.argv[2]
linkageMapFile = sys.argv[3]
contigLengthsFile = sys.argv[4]
contigMapFile = sys.argv[5]
haplotigLengthsFile = sys.argv[6]
linkageGroupValue = sys.argv[7]
LTR_GFF_file = sys.argv[8]
overlappingHaplotigList = sys.argv[9]
filterOverlapsBoolean = sys.argv[10]

#overlappingHaplotigFileName = overlappingHaplotigList.strip()
#fileName = os.path.basename(overlappingHaplotigFileName)
#overlappingHaplotigFileBaseName,fileExt = os.path.splitext(fileName)
#baseName,overlappingHaplotigThreshold = overlappingHaplotigFileBaseName.split('_')

scaleFactor = 10000
windowSize = 100000
megabaseMultiplier = 10

FILTER_OVERLAPS = None

if filterOverlapsBoolean == 'yes':
    FILTER_OVERLAPS = True
    overlapHaplotigDict = readOverlappingHaplotigList(overlappingHaplotigList)
    contigMapDict = readContigMapFile_withFiltering(contigMapFile,overlapHaplotigDict)
    overlappingHaplotigFileName = overlappingHaplotigList.strip()
    fileName = os.path.basename(overlappingHaplotigFileName)
    overlappingHaplotigFileBaseName,fileExt = os.path.splitext(fileName)
    baseName,overlappingHaplotigThreshold = overlappingHaplotigFileBaseName.split('_')
else:
    FILTER_OVERLAPS = False
    #overlapHaplotigDict = {}
    contigMapDict = readContigMapFile_noFiltering(contigMapFile)
    overlappingHaplotigThreshold = "noHapFilter"


# contigMapDict = readContigMapFile(contigMapFile,overlapHaplotigDict)
primaryLenDict = readContigLengthsFile(contigLengthsFile)
hapLenDict = readContigLengthsFile(haplotigLengthsFile)
linkageMapList,linkageGroupTotalLen = readLinkageMapFile(linkageMapFile,primaryLenDict,linkageGroupValue)

primaryCoordList,fig = drawPrimaryContigs(linkageMapList,primaryLenDict,scaleFactor,linkageGroupTotalLen,megabaseMultiplier,linkageGroupValue,overlappingHaplotigThreshold)
hapCoordList,fig = getHapCoords(primaryCoordList,fig,hapLenDict,contigMapDict,scaleFactor)

countArray = createCountArray(linkageGroupTotalLen,scaleFactor)
fig = drawHaplotigs(countArray,hapCoordList,fig,megabaseMultiplier)

geneData = readAugustusGFF(augustusGFF)
genePositionList = getHopGeneCoords(primaryCoordList,geneData,scaleFactor)
geneDensityList = computeGeneDensity(genePositionList,linkageGroupTotalLen,windowSize,scaleFactor,megabaseMultiplier)

ltrData = read_LTR_GFF_file(LTR_GFF_file)
ltrPositionList = getLTRCoords(primaryCoordList,ltrData,scaleFactor)
ltrDensityList = get_LTR_density(ltrPositionList,linkageGroupTotalLen,windowSize,scaleFactor,megabaseMultiplier)

createBarAndLineGraph(geneDensityList,ltrDensityList,windowSize,scaleFactor,megabaseMultiplier,linkageGroupTotalLen,linkageGroupValue,overlappingHaplotigThreshold)

snpDataList = read_SNP_file(snpFile,primaryCoordList)
snpPositionList = get_SNP_position(snpDataList,primaryCoordList)
hetDensityList = get_SNP_density(snpPositionList,linkageGroupTotalLen,windowSize,scaleFactor,megabaseMultiplier)
drawSNPs(snpPositionList,hetDensityList,fig,scaleFactor,megabaseMultiplier,linkageGroupTotalLen,linkageGroupValue,overlappingHaplotigThreshold)
