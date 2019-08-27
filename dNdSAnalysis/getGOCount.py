import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

# GO:0072552 000007F.g65.t1 000007F_066.g65.t1 0.272804891733 Q9CA61

def readGeneFile(geneFile):
    goCountDict = {}
    geneDict = {}
    with open(geneFile,'r') as GF:
        for line in GF:
            goID,primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov = line.strip().split('\t')
            # goID,primaryGeneID,hapGeneID,dNdS,uniprotID = line.strip().split('\t')
            if goID not in goCountDict:
                goCountDict[goID] = 0
            goCountDict[goID] += 1
            # print goID,goCountDict[goID]
            if goID not in geneDict:
                geneDict[goID] = []
            geneDict[goID].append((primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov))
            #print goID,primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov
    return(goCountDict,geneDict)

#Entry   Entry name      Gene ontology IDs       Gene ontology (biological process)      Gene ontology (cellular component)      Gene ontology (molecular function)      Gene ontology (GO)
# Q9FY74  CMTA1_ARATH     GO:0001077; GO:0005516; GO:0005634; GO:0009414; GO:0009733; GO:0043565; GO:0045893; GO:0045944; GO:0050826      positive regulation of transcription, DNA-templated

def readGOFile(goFile):
    goInfoDict = {}
    with open(goFile,'r') as GO:
        for line in GO:
            if 'Entry' not in line:
                line = line.strip().split('\t')
                uniprotID = line[0]
                uniprotName = line[1]
                goDescriptions = line[2:]
                #print goDescriptions
                goDescriptions = ''.join(goDescriptions)
                goDescriptions = goDescriptions.split(';')
                for goInfo in goDescriptions:
                    goInfo = goInfo.strip()
                    #print goInfo
                    if 'GO' in goInfo:
                        goInfo = goInfo.split('[GO:')
                        #print goInfo[0]
                        goDesc = goInfo[0]
                        #print goDesc
                        goTerm = goInfo[1]
                        #print goTerm
                        #goDesc,goTerm = goInfo.split('[')
                        goTerm = goTerm.replace(']','')
                        goTerm = goTerm.strip()
                        goTerm = "GO:" + goTerm
                        #print goTerm
                        if goTerm not in goInfoDict:
                            goInfoDict[goTerm] = goDesc
                            #print goTerm,goDesc
                    #else:
                        #print "no GO term"
                        #goInfoDict[goTerm] = "no GO term"
    return(goInfoDict)

def combine(goCountDict,goInfoDict,geneDict):
    for goID in goCountDict:
        goCount = goCountDict[goID]
        if goID in goInfoDict and goID in geneDict:
            goDesc = goInfoDict[goID]
            for primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov in geneDict[goID]:
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (uniprotID,goID,goDesc,goCount,primaryGeneID,hapGeneID,dNdS,strand,eValue,bitScore,queryCov))


##########
## MAIN ##
##########

if len(sys.argv) != 3 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <GO file> <gene file>"
    sys.exit()

# read in the input file 
goFile = sys.argv[1]
geneFile = sys.argv[2]

goCountDict,geneDict = readGeneFile(geneFile)
goInfoDict = readGOFile(goFile)
combine(goCountDict,goInfoDict,geneDict)
