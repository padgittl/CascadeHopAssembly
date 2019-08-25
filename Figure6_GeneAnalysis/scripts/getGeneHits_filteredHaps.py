import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

# 000000F_g100.t1 000000F_011_g100.t1     Plus    0       928/3   0       0       287/3   0       0.0     0.0     --
# non-synonymous substitutions, non-synonymous sites, non-synonymous substitutions/non-synonymous sites (pn), synonymous substitutions, synonymous sites, synonymous substitutions/synonymous sites (ps), dn, ds, dn/ds

def readFilteredContigMap(filteredContigMap):
    filteredContigMapDict = {}
    with open(filteredContigMap,'r') as FCM:
        for line in FCM:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in filteredContigMapDict:
                filteredContigMapDict[hapID] = True
    return(filteredContigMapDict)

def readDNDSFileList(dndsFileList,pepBlastpOutputFile,filteredContigMapDict,geneHitFile):
    with open(dndsFileList,'r') as FL:
        for dndsFile in FL:
            dndsFile = dndsFile.strip()
            dnds = parseDNDSFile(dndsFile,filteredContigMapDict)
            blastpHits = readGeneHitFile(geneHitFile)
            getHits2Uniprot(dnds,blastpHits)

def parseDNDSFile(fullFilePath,filteredContigMapDict):
    fileName = os.path.basename(fullFilePath)
    baseName,fileExt = os.path.splitext(fileName)
    #getContigIDs = re.search('(.+)_vs_(\d+F_\d+_g\d+.+|\d+F.+)_lastz',baseName)
    #getContigIDs = re.search('(.+)_vs_(.+)',baseName)
    #primaryID = getContigIDs.group(1)
    #hapID = getContigIDs.group(2)
    #primaryID = primaryID.replace('_g','.g')
    #hapID = hapID.replace('_g','.g')
    # 010078F_g5.t1
    # 003057F_004_g17.t1
    # 000000F_011_g100.t1_lastzPlusStrand_NT_dnds
    getStrand = re.search('_(lastz.+Strand)',baseName)
    strand = getStrand.group(1)
    dndsDict = {}
    with open(fullFilePath,'r') as F:
        for line in F:
            primaryGeneID,hapGeneID,strand,nonSynSubs,nonSynSites,pN,synSubs,synSites,pS,dN,dS,dNdS = line.strip().split('\t')
            # at this step, undefined values for dn/ds are removed
            primaryGeneID = primaryGeneID.replace('_g','.g')
            hapGeneID = hapGeneID.replace('_g','.g')
            # 000000F_011.g100.t1
            getHapID = re.search('(.+)\.g',hapGeneID)
            hapID = getHapID.group(1)
            if hapID in filteredContigMapDict:
                if not dNdS == '--' :
                    dN = float(dN)
                    dN = round(dN,5)
                    dS = float(dS)
                    dS = round(dS,5)
                    # this is to remove negative sign
                    dNdS = dNdS.replace('-','')
                    if primaryGeneID not in dndsDict:
                        dndsDict[primaryGeneID] = []
                    dndsDict[primaryGeneID].append((hapGeneID,dNdS))
                    #print primaryGeneID,hapGeneID
    return(dndsDict)

# 005429F.g9.t1 Q04903,Q38920,O80642,Q9LHL5,Q84J75
def readGeneHitFile(geneHitFile):
    blastpHitsDict = {}
    with open(geneHitFile,'r') as GHF:
        for line in GHF:
            geneID,blastHits = line.strip().split('\t')
            if geneID not in blastpHitsDict:
                blastpHitsDict[geneID] = blastHits
    return(blastpHitsDict)

def getHits2Uniprot(dndsDict,blastpHitsDict):
    for primaryGeneID in blastpHitsDict:
        uniprotIDs = blastpHitsDict[primaryGeneID]
        if primaryGeneID in dndsDict:
            for hapID,dNdS in dndsDict[primaryGeneID]:
                print("%s\t%s\t%s\t%s\t" % (primaryGeneID,hapID,dNdS,uniprotIDs))

##########
## MAIN ##
##########

if len(sys.argv) != 4 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <dnds file list> <gene hit file> <filtered contig map file>"
    sys.exit()

# read in the input file 
dndsFileList = sys.argv[1]
geneHitFile = sys.argv[2]
filteredContigMap = sys.argv[3]

filteredContigMapDict = readFilteredContigMap(filteredContigMap)
dnds = readDNDSFileList(dndsFileList,geneHitFile,filteredContigMapDict,geneHitFile)


