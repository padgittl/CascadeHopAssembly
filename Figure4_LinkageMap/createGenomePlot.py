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

sys.path.append('/svgwrite/svgwrite-1.2.1/')
import svgwrite
from svgwrite import cm, mm

###############
# SUBROUTINES #
###############

def readContigMapFile(contigMapFile):
    contigMapDict = {}
    with open(contigMapFile,'r') as CMF:
        for line in CMF:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            #print primaryID,hapID,primaryStart,primaryStop,scoreDens,source
            if primaryStart != '-1':
                if primaryID not in contigMapDict:
                    contigMapDict[primaryID] = []
                contigMapDict[primaryID].append((hapID,int(primaryStart),int(primaryStop),source))
    return(contigMapDict)


def readContigLengthsFile(contigLengthsFile):
    lenDict = {}
    with open(contigLengthsFile,'r') as LEN:
        for line in LEN:
            contigID,contigLen = line.strip().split()
            #print contigID
            if contigID not in lenDict:
                lenDict[contigID] = int(contigLen)
    return(lenDict)


def readLinkageMapFileAndGetOrderedPrimaryContigs(linkageMapFile,lenDict,linkageGroupValue):
    linkageMapList = []
    linkageGroupTotalLen = 0
    with open(linkageMapFile,'r') as QTL:
        for line in QTL:
            # Contig_Position,Lkgrp,pos(cM)
            if not line.startswith('Marker'):
                line = line.strip().split(',')
                contigPosStuff = line[0]
                # 578_221439F
                shortContigID,contigPosExtra = contigPosStuff.split('_')
                getContigPos = re.search('(\d+)',contigPosExtra)
                contigPos = getContigPos.group(1)
                contigPos = int(contigPos)  
                getShortenedContigID = re.search('(.+)_',contigPosStuff)
                shortenedContigID = getShortenedContigID.group(1)
                contigID = "0" * (6-len(str(shortenedContigID))) + str(shortenedContigID) + "F"
                linkageGroup = line[1]
                pos_cM = line[2]
                linkageGroup = int(linkageGroup)
                if linkageGroup == linkageGroupValue:
                    # lnkGrp6
                    #if 840 <= float(pos_cM) and float(pos_cM) <= 857:
                    if contigID in lenDict:
                        linkageMapList.append((contigID,linkageGroup,float(pos_cM)))
                        #print contigID
                    else:
                        print "problem with contigLenDict"
                        sys.exit()

    # this is to make sure that contigs that occur two in a row are not added to the primary contig list - this is to prevent gene density from being counted twice erroneously
    primaryContigList = []
    for i in range(len(linkageMapList)-1):
        contigID,linkageGroup,cMPos = linkageMapList[i]
        nextContigID,nextLinkageGroup,next_cMPos = linkageMapList[i+1]
        if contigID != nextContigID:
            primaryContigList.append(contigID)

    lastIndex = len(linkageMapList) - 1
    secondLastContigID,secondLastLinkageGroup,secondLast_cMPos = linkageMapList[-2]
    lastContigID,lastLinkageGroup,last_cMPos = linkageMapList[-1]

    if secondLastContigID != lastContigID:
        primaryContigList.append(lastContigID)

    # check occurrence of pattern with script -->
    # /nfs0/BB/Hendrix_Lab/Hops/version5/primary/GenomeDensityFigure/scripts/detectPatterns.py
    # pairs occur maximally twice in a row, at least in this case. for other linkage groups, will need further evaluation

    '''
    for primaryID in primaryContigList:
        print primaryID
    '''

    posDict = {}
    for i in range(len(primaryContigList)-3):
        contigID1 = primaryContigList[i]
        contigID2 = primaryContigList[i+1]
        contigID3 = primaryContigList[i+2]
        contigID4 = primaryContigList[i+3]
        if contigID1 == contigID3:
            # store position index as key to dictionary - storing contigID won't work because its presence in linkage group is not necessarily unique
            posDict[i+2] = contigID3
        if contigID2 == contigID4:
            posDict[i+3] = contigID4

    filteredPrimaryList = []
    for i in range(len(primaryContigList)):
        contigID = primaryContigList[i]
        if i not in posDict:
            #print primaryContigList[i]
            filteredPrimaryList.append(contigID)
            if contigID in lenDict:
                contigLen = lenDict[contigID] 
                linkageGroupTotalLen += contigLen
            else:
                print "contig not in length dictionary, subroutine 'parseMapFile'"
                sys.exit()
    
    '''
    for contigID in filteredPrimaryList:
        print contigID
    '''
    return(filteredPrimaryList,linkageGroupTotalLen)


def drawPrimaryContigs(filteredPrimaryList,lenDict,scaleFactor,linkageGroupTotalLen,linkageGroupValue):
    primaryContigSubset = {}
    #windowSubsetMin = 325000000
    #windowSubsetMax = 335000000
    windowSubsetMin = 255000000
    windowSubsetMax = 275000000
    primaryCoordList = []
    primaryContigStart = 0
    primaryContigStop = 0
    figWidth = '100in'
    figHeight = '30in'
    fig = svgwrite.Drawing(filename="contigDrawingLnkGrp" + str(linkageGroupValue) + ".svg", size=(figWidth, figHeight), profile='full')

    # draw scale bar for legend
    scaleBar = (float(10000000) / scaleFactor)
    scaleBar = fig.add(fig.rect((0,30), (scaleBar,5),fill='black',stroke='black',stroke_width='0.5'))

    for primaryID in filteredPrimaryList:
        if primaryID in lenDict:
            primaryContigLen = lenDict[primaryID]
            primaryContigStart = primaryContigStop
            primaryContigStop += primaryContigLen
            primaryCoordList.append((primaryID,primaryContigStart,primaryContigStop))
            if windowSubsetMin <= primaryContigStart and primaryContigStart < windowSubsetMax:
                if primaryID not in primaryContigSubset:
                    primaryContigSubset[primaryID] = []
                primaryContigSubset[primaryID].append((primaryContigStart,primaryContigStop,primaryContigLen))
                #print primaryID,primaryContigStart,primaryContigStart
                # primary contig coordinates scaled for visualization
                primaryScaledStart = primaryContigStart / scaleFactor
                primaryScaledStop = primaryContigStop / scaleFactor
                primaryScaledLen = primaryScaledStop - primaryScaledStart
                # fig.add(fig.rect((primaryScaledStart,0), (primaryScaledLen,20),fill='#2166ac',stroke='black',stroke_width='0.5'))
                fig.add(fig.rect((primaryScaledStart,0), (primaryScaledLen,75),fill='#2166ac',stroke='black',stroke_width='0.5'))
        else:
            print "problem with primary contig length dict"
            sys.exit()
    fig.save()
    return(primaryCoordList,primaryContigSubset,fig)

def getHapCoords(primaryCoordList,fig,hapLenDict,contigMapDict,scaleFactor):
    hapCoordList = []
    updatedHapStart = 0
    updatedHapStop = 0
    for primaryID,primaryContigStart,primaryContigStop in primaryCoordList:
        if primaryID in contigMapDict:
            for hapID,hapStart,hapStop,source in contigMapDict[primaryID]:
                #print primaryID,primaryContigStart,primaryContigStop,hapID,hapStart,hapStop,source
                # coordinates are 0-based 
                hapAlignmentLen = hapStop - hapStart
                updatedHapStart = primaryContigStart + hapStart
                updatedHapStop = updatedHapStart + hapAlignmentLen
                scaledHapStart = updatedHapStart / scaleFactor
                scaledHapStop = updatedHapStop / scaleFactor
                scaledHapLen = scaledHapStop - scaledHapStart
                hapCoordList.append((primaryID,hapID,updatedHapStart,updatedHapStop,scaledHapStart,scaledHapStop,scaledHapLen,source))
                #print primaryID,hapID,primaryContigStart,primaryContigStop
                #print primaryID,hapID,scaledHapStart,scaledHapStop,scaledHapLen,source
    #hapCoordList.sort(key=lambda x:x[4], reverse=True)
    return(hapCoordList,fig)


def calculateAssociateDensity(hapCoordList,linkageGroupTotalLen,windowSize,multiplier):
    hapDensityList = []
    hpcDensityList = []
    for i in range(0,linkageGroupTotalLen,windowSize):
        hapCount = 0
        hpcCount = 0
        position = float(i) + (windowSize / 2)
        windowInterval = i + windowSize
        for primaryID,hapID,updatedHapStart,updatedHapStop,scaledHapStart,scaledHapStop,scaledHapLen,source in hapCoordList:
            #print primaryID,hapID,position
            if source == "FALCON-unzip":
                if i <= updatedHapStart and updatedHapStart < windowInterval:
                    hapCount += 1
            elif source != "FALCON-unzip":
                if i <= updatedHapStart and updatedHapStart < windowInterval:
                    hpcCount += 1
                    #print primaryID,hapID,updatedHapStart,updatedHapStop
            else:
                print "problem with source"
                sys.exit()

        if hapCount > 0:
            megabaseHapCount = hapCount * multiplier
            hapDensityList.append((position,megabaseHapCount))
        if hapCount == 0:
            megabaseHapCount = 0
            hapDensityList.append((position,megabaseHapCount))
            
        if hpcCount > 0:
            megabaseHPCCount = hpcCount * multiplier
            hpcDensityList.append((position,megabaseHPCCount))
        if hpcCount == 0:
            megabaseHPCCount = 0
            hpcDensityList.append((position,megabaseHPCCount))

    return(hapDensityList,hpcDensityList)


def plotAssociateDensity(hapDensityList,hpcDensityList,windowSize,scaleFactor,linkageGroupTotalLen,linkageGroupValue,scaledWindow,scaledLinkageGroupTotalLen):
    fig, ax1 = plt.subplots(figsize=(4,0.5))
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0

    primaryPatch = mpatches.Patch(color='#2166ac', label='Primary Contigs')
    haplotigPatch = mpatches.Patch(color='#b2182b', label='FALCON-Unzip Haplotigs')
    HPCPatch = mpatches.Patch(color='#f4a582', label='Homologous Primary Contigs')

    hapDensityArray = np.array(hapDensityList)
    hpcDensityArray = np.array(hpcDensityList)

    for i in range(len(hapDensityArray)-2):
        position,hapDensity = hapDensityArray[i]
        nextPosition,nextHapDensity = hapDensityArray[i+1]
        axis1 = ax1.plot((position,nextPosition),(hapDensity,nextHapDensity), linewidth=1, color='#b2182b', marker='o', markersize=2)
        ax1.set_xlim(0,linkageGroupTotalLen)
        ax1.spines['left'].set_color('#b2182b')
        ax1.spines['right'].set_color('black')
        ax1.tick_params(axis='y', colors='#b2182b')
        # changes units from scientific notation to hundred of Mb 
        ax1.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x)/1e6, ',')))
    ax2 = ax1.twinx()
    for i in range(len(hpcDensityArray)-2):
        position,hpcDensity = hpcDensityArray[i]
        nextPosition,nextHPCDensity = hpcDensityArray[i+1]
        #print position,hpcDensity
        axis2 = ax2.plot((position,nextPosition),(hpcDensity,nextHPCDensity), linewidth=1, color='#f4a582', marker='o', markersize=2)
        ax2.set_xlim(0,linkageGroupTotalLen)
        ax2.spines['left'].set_color('black')
        ax2.spines['right'].set_color('#f4a582')
        ax2.tick_params(axis='y', colors='#f4a582')
        ax2.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x)/1e6, ',')))

        legendHapDensity = mlines.Line2D([], [], color='#b2182b',marker='o',markersize=6,label='FALCON-Unzip Haplotigs')
        legendHPCDensity = mlines.Line2D([], [], color='#f4a582',marker='o',markersize=6,label='Homologous Primary Contigs')

    plt.legend(handles=[legendHapDensity,legendHPCDensity], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("associateContigDensityLnkGrp" + str(linkageGroupValue) + ".svg", bbox_inches = 'tight', pad_inches=0)
    plt.close()


def createCountArray(scaledLinkageGroupTotalLen):
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

def drawHaplotigs(countArray,hapCoordList,primaryContigSubset,fig):
    yScalingFactor = -25
    windowSubsetMin = 255000000
    windowSubsetMax = 275000000
    # haplotig coordinates scaled for visualization
    for primaryID,hapID,updatedHapStart,updatedHapStop,scaledHapStart,scaledHapStop,scaledHapLen,source in hapCoordList:
        #print primaryID,hapID,updatedHapStart,updatedHapStop
        if primaryID in primaryContigSubset:
            #print primaryID
            for primaryContigStart,primaryContigStop,primaryContigLen in primaryContigSubset[primaryID]:
                if primaryContigStart <= updatedHapStart and updatedHapStop <= primaryContigStop:
                    #print primaryID,hapID,updatedHapStart,updatedHapStop
                    maxCount = getMaxCountValue(countArray,scaledHapStart,scaledHapStop)
                    # #b2182b
                    if source == "FALCON-unzip":
                        fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-85), (scaledHapLen,75),fill='#f4a582',stroke='black',stroke_width='0.5'))
                    # #f4a582
                    # #833437
                    # #742D33
                    if source == "LASTZ":
                        fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-85), (scaledHapLen,75),fill='#b2182b',stroke='black',stroke_width='0.5')) 
                    # #f4a582
                    # #833437
                    # #742D33
                    if source == "purgehaplotigs":
                        fig.add(fig.rect((scaledHapStart,maxCount*yScalingFactor-85), (scaledHapLen,75),fill='#b2182b',stroke='black',stroke_width='0.5'))
                    updatedMaxCount = maxCount + 1
                    countArray = getPosCount(countArray,scaledHapStart,scaledHapStop,updatedMaxCount)
                    maxCount = getMaxCountValue(countArray,scaledHapStart,scaledHapStop)
                    ##if '_' not in hapID:
                    #print primaryID,hapID,updatedHapStart,updatedHapStop
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


def getHopGeneCoords(primaryCoordList,geneData):
    genePositionList = []
    for primaryID,primaryContigStart,primaryContigStop in primaryCoordList:
        if primaryID in geneData:
            for geneID,strand,geneStart,geneStop in geneData[primaryID]:
                geneStart = int(geneStart)
                geneStop = int(geneStop)
                # coordinates are 1-based
                geneLen = geneStop - geneStart + 1
                updatedGeneStart = primaryContigStart + geneStart
                updatedGeneStop = updatedGeneStart + geneLen
                genePositionList.append((primaryID,geneID,strand,updatedGeneStart,updatedGeneStop,geneLen))
    genePositionList.sort(key=lambda x:x[3], reverse=False)
    return(genePositionList)


def computeGeneDensity(genePositionList,linkageGroupTotalLen,windowSize,multiplier):
    geneDensityList = []
    for i in range(0,linkageGroupTotalLen,windowSize):
        geneCount = 0
        #position = float(i)
        #windowInterval = i + windowSize
        position = float(i) + (windowSize / 2)
        windowInterval = i + windowSize
        for primaryID,geneID,strand,updatedGeneStart,updatedGeneStop,geneLen in genePositionList:
            if i <= updatedGeneStart and updatedGeneStart < windowInterval:
                geneCount += 1
        if geneCount > 0:
            megabaseGeneCount = geneCount * multiplier
            geneDensityList.append((position,megabaseGeneCount))
        else:
            geneCount = 0
            geneDensityList.append((position,geneCount))
    #for position,geneCount in geneDensityList:
    #    print position,geneCount
    return(geneDensityList)


def readLTRGFF(ltrGFF):
    ltrData = {}
    # 000000F RepeatMasker    LTR/Gypsy       2768    2924    31.9    +       418     003288F:275585..286585_INT   
    with open(ltrGFF,'r') as LTR:
        for line in LTR:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if contigID not in ltrData:
                    ltrData[contigID] = []
                ltrData[contigID].append((feature,int(start),int(end)))
    return(ltrData)


def getLTRCoords(primaryCoordList,ltrData):
    ltrPositionList = []
    for primaryID,primaryContigStart,primaryContigStop in primaryCoordList:
        if primaryID in ltrData:
            for LTR_type,ltrStart,ltrStop in ltrData[primaryID]:
                # coordinates are 1-based   
                ltrStart = int(ltrStart)
                ltrStop = int(ltrStop)
                ltrLen = ltrStop - ltrStart + 1
                updatedLTRStart = primaryContigStart + ltrStart 
                updatedLTRStop = updatedLTRStart + ltrLen
                ltrPositionList.append((primaryID,LTR_type,updatedLTRStart,updatedLTRStop,ltrLen))
    ltrPositionList.sort(key=lambda x:x[2], reverse=False)
    return(ltrPositionList)


def computeLTRDensity(ltrPositionList,linkageGroupTotalLen,windowSize,multiplier):
    ltrDensityList = []
    for i in range(0,linkageGroupTotalLen,windowSize):
        ltrCount = 0
        position = float(i)
        windowInterval = i + windowSize
        for primaryID,ltrType,updatedLTRStart,updatedLTRStop,ltrLen in ltrPositionList:
            if i <= updatedLTRStart and updatedLTRStart < windowInterval:
                ltrCount += 1
                #print primaryID,ltrType,scaledUpdatedLTRStart,scaledUpdatedLTRStop,scaledLTRLen,ltrCount
        if ltrCount > 0:
            megabaseLTRCount = ltrCount * multiplier
            ltrDensityList.append((position,megabaseLTRCount))
        else:
            ltrCount = 0
            ltrDensityList.append((position,ltrCount))
    #for position,ltrCount in ltrDensityList:
    #    print position,ltrCount
    return(ltrDensityList)


def plotGeneAndLTRDensity(geneDensityList,ltrDensityList,windowSize,scaleFactor,linkageGroupTotalLen,linkageGroupValue,scaledWindow,scaledLinkageGroupTotalLen):
    figWidth = 100
    figHeight = 30
    fig, ax1 = plt.subplots(figsize=(4,1))
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = figWidth
    fig_size[1] = figHeight
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0

    geneDensityArray = np.array(geneDensityList)
    ltrDensityArray = np.array(ltrDensityList)

    #print "ltr coords for plotting"
    #for position,ltrDensity in ltrDensityArray:
    for i in range(len(ltrDensityArray)-1):
        position,ltrDensity = ltrDensityArray[i]
        ax1.bar(position,ltrDensity,windowSize, color='#d6604d', align='edge')
        #ax1.set_xlim(0,linkageGroupTotalLen)
        ax1.set_ylim(5000,7200)
        ax1.spines['left'].set_color('#d6604d')
        ax1.spines['right'].set_color('black')
        ax1.tick_params(axis='y', colors='#d6604d')
        #print position,ltrDensity
        ax1.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x)/1e6, ',')))

    ax2 = ax1.twinx()
    #print "gene coords for plotting"
    for j in range(len(geneDensityArray)-2):
        position,geneDensity = geneDensityArray[j]
        nextPosition,nextGeneDensity = geneDensityArray[j+1]
        #print position,geneDensity,nextPosition,nextGeneDensity
        ax2.plot((position,nextPosition),(geneDensity,nextGeneDensity), linewidth=1, color='black', marker='o', markersize=2)
        # ax2.set_xlim(0,linkageGroupTotalLen)
        ax2.spines['left'].set_color('#d6604d')
        ax2.spines['right'].set_color('black')
        ax2.tick_params(axis='y', colors='black')
        ax2.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x)/1e6, ',')))

    legendGenes = mlines.Line2D([], [], color='black', marker='o',markersize=6,label='Gene Density')
    legendLTR = mlines.Line2D([], [], color='#d6604d', marker='s',linestyle="None",markersize=6, label='LTR Density')
    plt.legend(handles=[legendGenes,legendLTR], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("geneRepeatDensityLnkGrp" + str(linkageGroupValue) + ".svg", bbox_inches = 'tight', pad_inches=0)


def readSNPFile(snpFile,primaryCoordList):
    snpDataList = []
    with open(snpFile,'r') as SNP:
        for line in SNP:
            if not line.startswith('Site'):
                Site,SiteName,chrom,pos,numTaxa,Ref,Alt,MajorAllele,MajorAlleleGametes,MajorAlleleProportion,MajorAlleleFreq,MinorAllele,MinorAlleleGametes,MinorAlleleProportion,MinorAlleleFreq = line.strip().split(',')[:15]
                contigID = "0" * (6-len(str(chrom))) + str(chrom)
                contigID += "F"
                snpDataList.append((contigID,int(pos),float(MinorAlleleFreq)))
    return(snpDataList)


def getSNPPositions(snpDataList,primaryCoordList):
    snpPositionList = []
    for primaryContigID,primaryContigStart,primaryContigStop in primaryCoordList:
        for primaryID,snpPos,MinorAlleleFreq in snpDataList:
            if primaryContigID == primaryID:
                updatedSNPPosition = primaryContigStart + snpPos
                snpPositionList.append((primaryID,updatedSNPPosition,MinorAlleleFreq))
    return(snpPositionList)

def computeHetAndSNPDensity(snpPositionList,linkageGroupTotalLen,windowSize,multiplier):
    snpAndHetDensityList = []
    for i in range(0,linkageGroupTotalLen,windowSize):
        snpCount = 0
        heterozygositySum = 0
        position = float(i) + (windowSize / 2)
        windowInterval = i + windowSize
        for primaryID,updatedSNPPosition,minorAlleleFreq in snpPositionList:
            p = minorAlleleFreq
            heterozygosityValue = 2*p*(1.0-p)
            if i <= updatedSNPPosition and updatedSNPPosition < windowInterval:
                snpCount += 1
                heterozygositySum += heterozygosityValue
        if snpCount > 1:
            megabaseSNPCount = snpCount * multiplier
            megabaseHeterozygositySum = heterozygositySum * multiplier
            heterozygosityDensity = float(megabaseHeterozygositySum) / megabaseSNPCount
            snpAndHetDensityList.append((position,heterozygosityDensity,megabaseSNPCount))
            #print position,heterozygositySum,snpCount,heterozygosityDensity,snpCount,windowSize,snpDensity
        else:
            snpCount = 0
            heterozygosityDensity = 0
            snpAndHetDensityList.append((position,heterozygosityDensity,snpCount))
    return(snpAndHetDensityList)


def plotSNPs(snpPositionList,snpAndHetDensityList,fig,scaleFactor,linkageGroupValue,linkageGroupTotalLen):
    figWidth = 100
    figHeight = 30
    fig, ax1 = plt.subplots(figsize=(4,0.5))
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = figWidth
    fig_size[1] = figHeight
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    
    primaryPatch = mpatches.Patch(color='#2166ac', label='Primary Contigs')
    HPCPatch = mpatches.Patch(color='#f4a582', label='Homologous Primary Contigs')
    haplotigPatch = mpatches.Patch(color='#b2182b', label='FALCON-Unzip Haplotigs')

    snpAndHetDensityArray = np.array(snpAndHetDensityList)
    
    for i in range(len(snpAndHetDensityArray)-2):
        position,hetDensity,snpDensity = snpAndHetDensityArray[i]
        nextPosition,nextHetDensity,nextSNPDensity = snpAndHetDensityArray[i+1]
        axis1 = ax1.plot((position,nextPosition),(snpDensity,nextSNPDensity), linewidth=1, color='#4393c3', marker='o', markersize=2)
        ax1.set_xlim(0,linkageGroupTotalLen)
        ax1.spines['left'].set_color('#4393c3')
        ax1.spines['right'].set_color('black')
        ax1.tick_params(axis='y', colors='#4393c3')
        # changes units from scientific notation to hundred of Mb
        ax1.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x)/1e6, ',')))

    ax2 = ax1.twinx()
    for i in range(len(snpAndHetDensityArray)-2):
        position,hetDensity,snpCount = snpAndHetDensityArray[i]
        nextPosition,nextHetDensity,nextSNPCount = snpAndHetDensityArray[i+1]
        axis2 = ax2.plot((position,nextPosition),(hetDensity,nextHetDensity), linewidth=1, color='black', marker='o', markersize=2) 
        ax2.set_xlim(0,linkageGroupTotalLen)
        ax2.spines['left'].set_color('#4393c3')
        ax2.spines['right'].set_color('black')
        ax2.tick_params(axis='y', colors='black')
        ax2.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x)/1e6, ',')))

    legendSNPDensity = mlines.Line2D([], [], color='#4393c3',marker='o',markersize=6,label='SNP Density')
    legendAverageHet = mlines.Line2D([], [], color='black',marker='o',markersize=6,label='Average Heterozygosity')
    # fontsize='6'
    plt.legend(handles=[legendSNPDensity,legendAverageHet], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("snpDensityLnkGrp" + str(linkageGroupValue) + ".svg", bbox_inches = 'tight', pad_inches=0)
    plt.close()

    plt.legend(handles=[primaryPatch,HPCPatch,haplotigPatch], frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    plt.savefig("contigLegend.svg", bbox_inches = 'tight', pad_inches=0)


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <overlap-filtered contig map file> <primary contig lengths file> <haplotig lengths file> <snp position file> <linkage map file> <hop gene gff> <ltr gff> <linkage group number> <window size, e.g. 100000>\n"
if len(sys.argv) != 10:
    print usage
    sys.exit()


contigMapFile = sys.argv[1]
primaryContigLengthsFile = sys.argv[2]
haplotigLengthsFile = sys.argv[3]
snpFile = sys.argv[4]
linkageMapFile = sys.argv[5]
augustusGFF = sys.argv[6]
ltrGFF = sys.argv[7]
linkageGroupValue = sys.argv[8]
windowSize = sys.argv[9]


windowSize = int(windowSize)
scaleFactor = 10000
linkageGroupValue = int(linkageGroupValue)
scaledWindow = windowSize / scaleFactor
#megabaseMultiplier = 1000000.0 / windowSize
multiplier = 1

contigMapDict = readContigMapFile(contigMapFile)
primaryLenDict = readContigLengthsFile(primaryContigLengthsFile)
hapLenDict = readContigLengthsFile(haplotigLengthsFile)
filteredPrimaryList,linkageGroupTotalLen = readLinkageMapFileAndGetOrderedPrimaryContigs(linkageMapFile,primaryLenDict,linkageGroupValue)
scaledLinkageGroupTotalLen = linkageGroupTotalLen / scaleFactor
#print linkageGroupTotalLen

primaryCoordList,primaryContigSubset,fig = drawPrimaryContigs(filteredPrimaryList,primaryLenDict,scaleFactor,linkageGroupTotalLen,linkageGroupValue)
hapCoordList,fig = getHapCoords(primaryCoordList,fig,hapLenDict,contigMapDict,scaleFactor)
countArray = createCountArray(scaledLinkageGroupTotalLen)
fig = drawHaplotigs(countArray,hapCoordList,primaryContigSubset,fig)

# plot associate contig density curve
hapDensityList,hpcDensityList = calculateAssociateDensity(hapCoordList,linkageGroupTotalLen,windowSize,multiplier)
plotAssociateDensity(hapDensityList,hpcDensityList,windowSize,scaleFactor,linkageGroupTotalLen,linkageGroupValue,scaledWindow,scaledLinkageGroupTotalLen)

geneData = readAugustusGFF(augustusGFF)
genePositionList = getHopGeneCoords(primaryCoordList,geneData)
geneDensityList = computeGeneDensity(genePositionList,linkageGroupTotalLen,windowSize,multiplier)
ltrData = readLTRGFF(ltrGFF)
ltrPositionList = getLTRCoords(primaryCoordList,ltrData)
ltrDensityList = computeLTRDensity(ltrPositionList,linkageGroupTotalLen,windowSize,multiplier)
plotGeneAndLTRDensity(geneDensityList,ltrDensityList,windowSize,scaleFactor,linkageGroupTotalLen,linkageGroupValue,scaledWindow,scaledLinkageGroupTotalLen)

snpDataList = readSNPFile(snpFile,primaryCoordList)
snpPositionList = getSNPPositions(snpDataList,primaryCoordList)
snpAndHetDensityList = computeHetAndSNPDensity(snpPositionList,linkageGroupTotalLen,windowSize,multiplier)
plotSNPs(snpPositionList,snpAndHetDensityList,fig,scaleFactor,linkageGroupValue,linkageGroupTotalLen)
