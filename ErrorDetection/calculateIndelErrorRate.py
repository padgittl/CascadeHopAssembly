import sys, os, re
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############


# this subroutine reads in the cds fasta from transdecoder
def readFasta(fastaFile):
    cdsRecordDict = {}
    for record in SeqIO.parse(fastaFile,"fasta"):
        getTranscriptID = re.search('::(TCONS_\d+)::',record.id)
        transcriptID = getTranscriptID.group(1)
        getShortenedGeneID = re.search('(.+g\.\d+)::m',record.id)
        shortenedGeneID = getShortenedGeneID.group(1)
        cdsRecordDict[transcriptID] = 1
        #print record.id,transcriptID,shortenedGeneID
    return(cdsRecordDict)


def readLengthsFile(lengthsFile):
    lengthsDict = {}
    totalLength = 0
    with open(lengthsFile,'r') as L:
        for line in L:
            contigID,contigLength = line.strip().split()
            contigLength = int(contigLength)
            lengthsDict[contigID] = contigLength
    return(lengthsDict)


def createCountArray(contigLen):
    countArray = [0]*contigLen
    return(countArray)

def updateCount(countArray,start,stop,status):
    for i in range(start,stop):
        if countArray[i] == 0:
            countArray[i] = status
    return(countArray)


# transcript gff file from cuffmerge
def readGFFFile(gffFile,cdsRecordDict):
    coordDict = {}
    transcriptDict = {}
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                # ID=TCONS_00000073;geneID=XLOC_000052
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                end = int(end)
                start = int(start)
                featureLen = end - start + 1
                if feature == 'exon':
                    getTranscriptID = re.search('Parent=(TCONS_\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID in cdsRecordDict:
                        transcriptDict[transcriptID] = 1
                        #print transcriptID,attribute
                        if contigID not in coordDict:
                            coordDict[contigID] = []
                        coordDict[contigID].append((start,end))
    for contigID in coordDict:
        coordDict[contigID].sort(key=lambda x:x[0], reverse=False)
    print len(transcriptDict)
    return(coordDict)

def createGeneArray(coordDict,lengthsDict):
    totalFeatureLenDict = {}
    for contigID in coordDict:
        contigLen = lengthsDict[contigID]
        countArray = createCountArray(contigLen)
        for start,end in coordDict[contigID]:
            featureStatus = 1
            countArray = updateCount(countArray,start,end,featureStatus)
        featureCount = countArray.count(featureStatus)
        #print contigID,featureCount
        if contigID not in totalFeatureLenDict:
            totalFeatureLenDict[contigID] = featureCount
    return(totalFeatureLenDict)

def readVCF(vcfFile):
    vcfDict = {}
    with open(vcfFile,'r') as VCF:
        for line in VCF:
            if not line.startswith('#'):
                # 000000F 877276  .       A       G       250.77  PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-3.224;ClippingRankSum=0.000;DP
                line = line.strip().split("\t")
                contigID = line[0]
                variantPos = line[1]
                variantPos = int(variantPos)
                idInfo = line[2]
                ref = line[3]
                alt = line[4]
                qual = line[5]
                filterStatus = line[6]
                info = line[7]
                formatValues = line[8] 
                genotypeFields = line[9:]
                genotypeFieldsJoined = '\t'.join(genotypeFields)
                if len(ref) != 1 or len(alt) != 1:
                    #print contigID,variantPos,ref,alt
                    if contigID not in vcfDict:
                        vcfDict[contigID] = []
                    vcfDict[contigID].append(variantPos)
    return(vcfDict)

def assessVariantOverlap(coordDict,vcfDict):
    variantTotalCountDict = {}
    for contigID in vcfDict:
        if contigID in coordDict:
            # every variant for a contig
            for variantPos in vcfDict[contigID]:
                for featureStart,featureEnd in coordDict[contigID]:
                    if featureStart <= variantPos and variantPos < featureEnd:
                        if contigID not in variantTotalCountDict:
                            variantTotalCountDict[contigID] = 0
                        variantTotalCountDict[contigID] += 1
                        break
                    elif variantPos < featureStart:
                        break
    #for contigID in variantTotalCountDict:
    #    variantCount = variantTotalCountDict[contigID]
    #    #print variantCount
    return(variantTotalCountDict)

def calculateErrorInFeature(variantTotalCountDict,totalFeatureLenDict,lengthsDict):
    for contigID in variantTotalCountDict:
        contigLen = lengthsDict[contigID]
        outName = contigID + '.error.txt'
        OUT = open(outName,'w')
        variantCount = variantTotalCountDict[contigID]
        if contigID in totalFeatureLenDict:
            totalFeatureLen = totalFeatureLenDict[contigID]
            percentError = round(float(variantCount)/totalFeatureLen * 100, 3)
            OUT.write("%s\t%s\t%s\t%s\t%s\n" % (contigID,percentError,variantCount,totalFeatureLen,contigLen))
            #print contigID,percentError
        else:
            print "missing contig"
            sys.exit()

###########
# MAIN ####
###########

usage = "Usage: " + sys.argv[0] + " <gene gff file> <vcf file> <lengths file> <fasta file>\n"
if len(sys.argv) != 5:
    print usage
    sys.exit()

gffFile = sys.argv[1]
vcfFile = sys.argv[2]
lengthsFile = sys.argv[3]
fastaFile = sys.argv[4]

cdsRecordDict = readFasta(fastaFile)
lengthsDict = readLengthsFile(lengthsFile)
coordDict = readGFFFile(gffFile,cdsRecordDict)
totalFeatureLenDict = createGeneArray(coordDict,lengthsDict)
vcfDict = readVCF(vcfFile)

variantTotalCountDict = assessVariantOverlap(coordDict,vcfDict)
calculateErrorInFeature(variantTotalCountDict,totalFeatureLenDict,lengthsDict)
