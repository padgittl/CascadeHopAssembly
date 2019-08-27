import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

def readGeneFile(geneFile):
    geneDetailsDict = {}
    with open(geneFile,'r') as GF:
        for line in GF:
            primaryGeneID,hapGeneID,source,strand,uniprotID,eValue,bitScore,queryCov,dN,dS,dNdS = line.strip().split('\t')
            if uniprotID not in geneDetailsDict:
                geneDetailsDict[uniprotID] = []
            geneDetailsDict[uniprotID].append((primaryGeneID,hapGeneID,source,strand,uniprotID,float(eValue),float(bitScore),float(queryCov),dN,dS,float(dNdS)))
    return(geneDetailsDict)

# Entry   Entry name      Gene ontology IDs
# Q9FY74  CMTA1_ARATH     GO:0001077; GO:0005516; GO:0005634; GO:0009414; GO:0009733; GO:0043565; GO:0045893; GO:0045944; GO:0050826

def readGOTermFile(goFile):
    uniprotDict = {}
    with open(goFile,'r') as GO:
        for line in GO:
            if not line.startswith('Entry'):
                #uniprotID,uniprotName,goTerms = line.strip().split('\t')
                #print goTerms
                line = line.strip().split('\t')
                uniprotID = line[0]
                uniprotName = line[1]
                goTerms = line[2:]
                goTerms = ''.join(goTerms)
                goTerms = goTerms.split(';')
                #print goTerms
                for goTerm in goTerms:
                    if uniprotID not in uniprotDict:
                        uniprotDict[uniprotID] = []
                    uniprotDict[uniprotID].append((goTerm))
    return(uniprotDict)

def linkUniprotAndGO(geneDetailsDict,uniprotDict):
    for uniprotID in geneDetailsDict:
        if uniprotID in uniprotDict:
            for primaryGeneID,hapGeneID,source,strand,uniprotID,eValue,bitScore,queryCov,dN,dS,dNdS in geneDetailsDict[uniprotID]:
                goTerms = uniprotDict[uniprotID]
                goTerms = ';'.join(goTerms)
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (primaryGeneID,hapGeneID,source,strand,uniprotID,eValue,bitScore,queryCov,dN,dS,dNdS,goTerms))
        else:
            print "uniprotID has no go term"
            sys.exit()

##########
## MAIN ##
##########

if len(sys.argv) != 3 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <gene file containing dNdS and uniprotIDs> <GO term file> "
    sys.exit()

# read in the input file 
geneFile = sys.argv[1]
goFile = sys.argv[2]

geneDetailsDict = readGeneFile(geneFile)
uniprotDict = readGOTermFile(goFile)
#for uniprotID in uniprotDict:
#    for goTerms in uniprotDict[uniprotID]:
#        print uniprotID,goTerms

linkUniprotAndGO(geneDetailsDict,uniprotDict)
