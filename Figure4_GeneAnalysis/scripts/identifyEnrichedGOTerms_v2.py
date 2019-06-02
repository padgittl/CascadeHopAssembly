import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

# 000018F.g99.t1  000018F_006.g99.t1      0.0     A8MQY1  GO:0003677; GO:0003700; GO:0005634; GO:0005789; GO:0009736; GO:0009965; GO:0016021; GO:0031965; GO:0045893; GO:0051302

def readGeneFileWithGOTerms(geneFileWithGOTerms):
    goTermDict = {}
    with open(geneFileWithGOTerms,'r') as GF:
        for line in GF:
            line = line.strip().split('\t')
            primaryGeneID = line[0]
            hapGeneID = line[1]
            dNdS = line[2]
            uniprotID = line[3]
            goTermDescriptions = line[4:]
            goTermDescriptions = ''.join(goTermDescriptions)
            if ';' in goTermDescriptions:
                goTermDescriptions = goTermDescriptions.split(';')
            for goTerm in goTermDescriptions:
                goTerm = goTerm.strip()
                #print goTerm
                if goTerm not in goTermDict:
                    goTermDict[goTerm] = []
                goTermDict[goTerm].append((primaryGeneID,hapGeneID,dNdS,uniprotID))
                #print goTerm,primaryGeneID,hapGeneID,dNdS,uniprotID
    return(goTermDict)
                #else:
                #    print "GO term error: " + goTermInfo
                    #sys.exit()

##########
## MAIN ##
##########

if len(sys.argv) != 2 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <gene file containing dNdS, uniprotIDs, and GO terms> "
    sys.exit()

# read in the input file 
geneFileWithGOTerms = sys.argv[1]

goTermDict = readGeneFileWithGOTerms(geneFileWithGOTerms)
for goTerm in goTermDict:
    #print goTerm
    for primaryGeneID,hapGeneID,dNdS,uniprotID in goTermDict[goTerm]:
        print("%s\t%s\t%s\t%s\t%s" % (goTerm,primaryGeneID,hapGeneID,dNdS,uniprotID))

