import sys, os, re

###############
# SUBROUTINES #
###############

# read in contig-specific error rate file list
def readFileListAndCalculateErrorTotal(errorFileList):
    variantTotal = 0
    featureTotal = 0
    with open(errorFileList,'r') as F:
        for line in F:
            fileName = line.strip()
            variantCount,totalFeatureLen = readErrorFile(fileName)
            variantTotal += variantCount
            featureTotal += totalFeatureLen
    errorPercentage = round((float(variantTotal) / featureTotal), 6) * 100
    print errorPercentage

# parse contig-specific error rate file
def readErrorFile(errorFile):
    with open(errorFile,'r') as E:
        for line in E:
            contigID,percentError,variantCount,totalFeatureLen,contigLen = line.strip().split('\t')
            variantCount = int(variantCount)
            totalFeatureLen = int(totalFeatureLen)
            #print contigID,percentError,variantCount,totalFeatureLen,contigLen
    return(variantCount,totalFeatureLen)


###########
# MAIN ####
###########

usage = "Usage: " + sys.argv[0] + " <error file list>\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

errorFileList = sys.argv[1]

readFileListAndCalculateErrorTotal(errorFileList)
