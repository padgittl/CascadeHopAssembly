import sys, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from matplotlib.ticker import AutoMinorLocator

################
# SUBROUTINES ##
################

# 000000F 8249941 
# type    count   percentage
# DNA     13671   0.166   
# LTR     5828708 70.652  
# LINE    2268    0.027   
# SINE    0       0.0     
# Satellite       0       0.0     
# Simple  128689  1.56    
# Mobile Element  11290   0.137   
# rRNA    208     0.003   
# Other   370     0.004   
# Total repeat content    5985204 72.549  
# Non-repeat      2264737 27.451

def readFileList(gffFileList,assemblyLen):

    dnaCount = 0
    ltrCount = 0
    lineCount = 0
    sineCount = 0
    satelliteCount = 0
    simpleCount = 0
    mobileElementCount = 0
    rRNACount = 0
    otherCount = 0
    totalRepeatCount = 0
    nonRepeatCount = 0

    with open(gffFileList,'r') as FL:
        for line in FL:
            fileName = line.strip()
            dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount = readRepeatFile(fileName,dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount)

        percDNA = round(float(dnaCount) / assemblyLen * 100,3)
        percLTR = round(float(ltrCount) / assemblyLen * 100,3)
        percLINE = round(float(lineCount) / assemblyLen * 100,3)
        percSINE = round(float(sineCount) / assemblyLen * 100,3)
        percSatellite = round(float(satelliteCount) / assemblyLen * 100,3)
        percSimple = round(float(simpleCount) / assemblyLen * 100,3)
        percMobileElement = round(float(mobileElementCount) / assemblyLen * 100,3)
        percrRNA = round(float(rRNACount) / assemblyLen * 100,3)
        percOTher = round(float(otherCount) / assemblyLen * 100,3)

        totalRepeatPercentage = percDNA + percLTR + percLINE + percSINE + percSatellite + percSimple + percMobileElement + percrRNA + percOTher
        nonRepeatPercentage = 100 - totalRepeatPercentage

        print("PacBioHopCascadeAssembly\t%s\t" % (assemblyLen))
        print("type\tcount\tpercentage")
        print("DNA\t%s\t%s\t" % (dnaCount,percDNA))
        print("LTR\t%s\t%s\t" % (ltrCount,percLTR))
        print("LINE\t%s\t%s\t" % (lineCount,percLINE))
        print("SINE\t%s\t%s\t" % (sineCount,percSINE))
        print("Satellite\t%s\t%s\t" % (satelliteCount,percSatellite))
        print("Simple\t%s\t%s\t" % (simpleCount,percSimple))
        print("Mobile Element\t%s\t%s\t" % (mobileElementCount,percMobileElement))
        print("rRNA\t%s\t%s\t" % (rRNACount,percrRNA))
        print("Other\t%s\t%s\t" % (otherCount,percOTher))
        print("Total Repeat Content\t%s\t%s\t" % (totalRepeatCount,totalRepeatPercentage))
        print("Non-Repeat\t%s\t%s\t" % (nonRepeatCount,nonRepeatPercentage))


def readRepeatFile(repeatFile,dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount):
    with open(repeatFile,'r') as RF:
        for line in RF:
            if 'F' in line:
                contigID,contigLength = line.strip().split('\t')
            if 'DNA' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                dnaCount += repCount
            if 'LTR' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                ltrCount += repCount
            if 'LINE' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                lineCount += repCount
            if 'SINE' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                sineCount += repCount
            if 'Satellite' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                satelliteCount += repCount
            if 'Simple' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                simpleCount += repCount
            if 'Mobile' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                mobileElementCount += repCount
            if 'rRNA' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                rRNACount += repCount
            if 'Other' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                otherCount += repCount
            if 'Total repeat content' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                totalRepeatCount += repCount
            if 'Non-repeat' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                nonRepeatCount += repCount
    return(dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount)


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gff file list> <assembly length>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

gffFileList = sys.argv[1]
assemblyLen = sys.argv[2]

assemblyLen = int(assemblyLen)

readFileList(gffFileList,assemblyLen)
