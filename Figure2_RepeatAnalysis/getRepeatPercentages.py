import sys, re, os
import numpy as np


################
# SUBROUTINES ##
################

def readLengthsFile(lengthsFile):
    lengthsDict = {}
    totalLength = 0
    with open(lengthsFile,'r') as L:
        for line in L:
            scaffoldID,scaffoldLength = line.strip().split()
            scaffoldLength = int(scaffoldLength)
            lengthsDict[scaffoldID] = scaffoldLength
    return(lengthsDict)

def createCountArray(scaffoldLen):
    countArray = [0]*scaffoldLen
    return(countArray)

def updateCount(countArray,start,stop,repState):
    for i in range(start,stop):
        if countArray[i] == 0:
            countArray[i] = repState
    return(countArray)

def readHopGFF(hopGFF):
    coordDict = {}
    with open(hopGFF,'r') as HOP:
        for line in HOP:
            if not line.startswith('#'):
                scaffoldID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                start = int(start)
                end = int(end)
                repLen = abs(end-start)
                #print start,end,repLen
                if scaffoldID not in coordDict:
                    coordDict[scaffoldID] = []
                coordDict[scaffoldID].append((start,end,repLen,float(score),attribute))
    # sort for each scaffoldID
    for scaffoldID in coordDict:
        coordDict[scaffoldID].sort(key=lambda x:x[3], reverse=True)
    return(coordDict)

def getPercentages(coordDict):

    DNA = 1
    LTR = 2
    LINE = 3
    SINE = 4
    SATELLITE = 5
    SIMPLE = 6
    MOBILE_ELEMENT = 7
    rRNA = 8
    OTHER = 9

    dnaCount = 0
    ltrCount = 0
    lineCount = 0
    sineCount = 0
    satelliteCount = 0
    simpleCount = 0
    mobileElementCount = 0
    rRNACount = 0
    otherCount = 0

    for scaffoldID in coordDict:
        scaffoldLen = lengthsDict[scaffoldID]
        countArray = createCountArray(scaffoldLen)
        for repStart,repStop,repeatLen,repScore,repAttribute in coordDict[scaffoldID]:
            if 'Motif' in repAttribute:
                getRepeatType = re.search('Name=Motif:(.+);',repAttribute)
                repType = getRepeatType.group(1)
                if 'DNA' in repType:
                    repState = DNA 
                elif 'LINE' in repType:
                    repState = LINE
                elif 'SINE' in repType:
                    repState = SINE
                elif 'MobileElement' in repType:
                    repState = MOBILE_ELEMENT
                elif 'rRNA' in repType:
                    repState = rRNA
                elif 'Other' in repType:
                    repState = OTHER
                elif repType == "Unknown":
                    repState = OTHER
                elif 'Satellite' in repType:
                    repState = SATELLITE
                elif '(' in repType:
                    repState = SIMPLE
                elif 'Retroelement' in repType:
                    repState = LTR
                elif 'LTR' in repType:
                    repState = LTR
                else:
                    print "Warning: " + repAttribute + " unaccounted for"
                    sys.exit()
            elif 'LTR' in repAttribute:
                getRepeatType = re.search('Name=(LTR.+);',repAttribute)
                repType = getRepeatType.group(1)
                if 'LTR/Gypsy' in repType:
                    repState = LTR
                elif 'LTR/Copia' in repType:
                    repState = LTR
                elif 'LTR/unknown' in repType:
                    repState = LTR
                else:
                    print "Warning: " + repAttribute + " unaccounted for"
                    sys.exit()
            countArray = updateCount(countArray,repStart,repStop,repState)
        dnaCount = countArray.count(DNA)
        ltrCount = countArray.count(LTR)
        lineCount = countArray.count(LINE)
        sineCount = countArray.count(SINE)
        satelliteCount = countArray.count(SATELLITE)
        simpleCount = countArray.count(SIMPLE)
        mobileElementCount = countArray.count(MOBILE_ELEMENT)
        rRNACount = countArray.count(rRNA)
        otherCount = countArray.count(OTHER)
    return(dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gff file> <scaffold lengths file>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

hopGFF = sys.argv[1]
lengthsFile = sys.argv[2]

lengthsDict = readLengthsFile(lengthsFile)
coordDict = readHopGFF(hopGFF)
dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount = getPercentages(coordDict)
#print dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount

baseName = os.path.basename(hopGFF)
scaffoldID,fileExt = os.path.splitext(baseName)

if scaffoldID in lengthsDict:
    scaffoldLen = lengthsDict[scaffoldID]
    # print scaffoldID,scaffoldLen

percDNA = round(float(dnaCount) / scaffoldLen * 100,3)
percLTR = round(float(ltrCount) / scaffoldLen * 100,3)
percLINE = round(float(lineCount) / scaffoldLen * 100,3)
percSINE = round(float(sineCount) / scaffoldLen * 100,3)
percSatellite = round(float(satelliteCount) / scaffoldLen * 100,3)
percSimple = round(float(simpleCount) / scaffoldLen * 100,3)
percMobileElement = round(float(mobileElementCount) / scaffoldLen * 100,3)
percrRNA = round(float(rRNACount) / scaffoldLen * 100,3)
percOTher = round(float(otherCount) / scaffoldLen * 100,3)

totalRepeatCount = dnaCount + ltrCount + lineCount + sineCount + satelliteCount + simpleCount + mobileElementCount + rRNACount + otherCount
totalRepeatPercentage = percDNA + percLTR + percLINE + percSINE + percSatellite + percSimple + percMobileElement + percrRNA + percOTher
nonRepeatCount = scaffoldLen - totalRepeatCount
nonRepeatPercentage = 100 - totalRepeatPercentage


print("%s\t%s\t" % (scaffoldID,scaffoldLen))
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
print("Total repeat content\t%s\t%s\t" % (totalRepeatCount,totalRepeatPercentage))
print("Non-repeat\t%s\t%s\t" % (nonRepeatCount,nonRepeatPercentage))


