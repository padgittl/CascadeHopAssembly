import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#################
## SUBROUTINES ##
#################

# primaryGeneID,hapGeneID,source,strand,uniprotID,eValue,bitScore,queryCov,dN,dS,dNdS,goTerms 

def readGeneFileWithGOTerms(geneFileWithGOTerms):
    goTermDict = {}
    with open(geneFileWithGOTerms,'r') as GF:
        for line in GF:
            # primaryGeneID,hapGeneID,source,strand,uniprotID,eValue,bitScore,queryCov,dN,dS,dNdS,goTerms
            line = line.strip().split('\t')
            primaryGeneID = line[0]
            hapGeneID = line[1]
            source = line[2]
            strand = line[3]
            uniprotID = line[4]
            eValue = line[5]
            bitScore = line[6]
            queryCov = line[7]
            dNdS = line[10]
            # print "primaryGeneID,hapGeneID,source,uniprotID,dNdS"
            # print primaryGeneID,hapGeneID,source,uniprotID,dNdS
            goTermDescriptions = line[11:]
            goTermDescriptions = ''.join(goTermDescriptions)
            if ';' in goTermDescriptions:
                goTermDescriptions = goTermDescriptions.split(';')
                for goTerm in goTermDescriptions:
                    goTerm = goTerm.strip()
                    if goTerm not in goTermDict:
                        goTermDict[goTerm] = []
                    goTermDict[goTerm].append((primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov))
                    #print goTerm,primaryGeneID,hapGeneID,dNdS,uniprotID
            else:
                if 'GO' in goTermDescriptions:
                    goTerm = goTermDescriptions.strip()
                    if goTerm not in goTermDict:
                        goTermDict[goTerm] = []
                    goTermDict[goTerm].append((primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov))
                #else:
                #    # print goTermDescriptions
                #    continue
    return(goTermDict)

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
    for primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov in goTermDict[goTerm]:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (goTerm,primaryGeneID,hapGeneID,dNdS,uniprotID,strand,eValue,bitScore,queryCov))

