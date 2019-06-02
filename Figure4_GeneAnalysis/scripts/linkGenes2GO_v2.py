import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

# 000000F.g98.t4  000000F_011.g98.t4      0.0     Q43088,Q9XI84,P94026   

def readGeneFile(geneFile):
    geneDetailsDict = {}
    with open(geneFile,'r') as GF:
        for line in GF:
            primaryGeneID,hapGeneID,dNdS,uniprotIDs = line.strip().split('\t')
            uniprotIDs = uniprotIDs.split(',')
            for uniprotID in uniprotIDs:
                if uniprotID not in geneDetailsDict:
                    geneDetailsDict[uniprotID] = []
                geneDetailsDict[uniprotID].append((primaryGeneID,hapGeneID,dNdS))
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
        for primaryGeneID,hapGeneID,dNdS in geneDetailsDict[uniprotID]:
            if uniprotID in uniprotDict:
                goTerms = uniprotDict[uniprotID]
                goTerms = ';'.join(goTerms)
                print("%s\t%s\t%s\t%s\t%s\t" % (primaryGeneID,hapGeneID,dNdS,uniprotID,goTerms))
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
