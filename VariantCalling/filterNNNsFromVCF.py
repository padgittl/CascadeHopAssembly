import sys, os, re

###############
# SUBROUTINES #
###############

# contig lengths file
def readLengthsFile(lengthsFile):
    lengthsDict = {}
    totalLength = 0
    with open(lengthsFile,'r') as L:
        for line in L:
            contigID,contigLength = line.strip().split()
            contigLength = int(contigLength)
            lengthsDict[contigID] = contigLength
    return(lengthsDict)


# vcf file created from variant calling pipeline
def readVCF(vcfFile):
    NDict = {}
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
                if 'N' in ref:
                    if contigID not in NDict:
                        NDict[contigID] = []
                    NDict[contigID].append((variantPos,ref,alt))
                qual = line[5]
                filterStatus = line[6]
                info = line[7]
                formatValues = line[8] 
                genotypeFields = line[9:]
                genotypeFieldsJoined = '\t'.join(genotypeFields)
    return(NDict)


# see how many Ns are present per variant
def assessNs(NDict):
    variantsToExclude = {}
    totalNCount = 0
    for contigID in NDict:
        for variantPos,ref,alt in NDict[contigID]:
            nCount = 0
            if len(ref) == 1:
                variantsToExclude[ref] = 1
                totalNCount += 1
            elif len(ref) > 1:
                variantsToExclude[ref] = 1
                totalNCount += 1
                for i in ref:
                    if i == 'N':
                        nCount += 1
                percentN = float(nCount) / len(ref) * 100
                #if percentN <= 50:
                #    print contigID,ref,alt,nCount,len(ref),percentN
            else:
                print "problem"
                sys.exit()
    #print totalNCount
    #print len(variantsToExclude)
    return(variantsToExclude)


# take in dictionary of N-containing variants and read through vcf file again
def filterVCF(vcfFile,variantsToExclude):
    OUT = open('variantsFilteredByN.vcf','w')
    with open(vcfFile,'r') as VCF:
        for line in VCF:
            if not line.startswith('#'):
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
                # if N in variant, exclude from filtered file
                if ref not in variantsToExclude:
                    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (contigID,variantPos,idInfo,ref,alt,qual,filterStatus,info,formatValues,genotypeFieldsJoined))
                    

###########
# MAIN ####
###########


usage = "Usage: " + sys.argv[0] + " <gene gff file> <vcf file> <lengths file>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

gffFile = sys.argv[1]
vcfFile = sys.argv[2]
lengthsFile = sys.argv[3]

lengthsDict = readLengthsFile(lengthsFile)
NDict = readVCF(vcfFile)
variantsToExclude = assessNs(NDict)
filterVCF(vcfFile,variantsToExclude)

