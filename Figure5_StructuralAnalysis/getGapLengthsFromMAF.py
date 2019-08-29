import sys,re,os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#import numpy as np

def parseMAF(mafFile,outFile1,outFile2):
    list1 = []
    list2 = []
    #000201F_vs_009864F_lastzMinusStrand.txt 
    fullPath = mafFile.strip()
    fileName = os.path.basename(fullPath)
    baseName,fileExt = os.path.splitext(fileName)
    # 000201F_vs_009864F_lastzMinusStrand
    getContigIDs = re.search('(\d+F_\d+)_vs_(\d+F)',baseName)
    contigID1 = getContigIDs.group(1)
    contigID2 = getContigIDs.group(2)    
    contigInfo1 = []
    contigInfo2 = []
    with open(mafFile,'r') as F:
        OUT1 = open(outFile1,'w')
        OUT2 = open(outFile2,'w')
        for line in F:
            if line.startswith('#'):
                continue
            if line.startswith('a'):
                a,score = line.strip().split()
            if line.startswith('s'):
                s,seqName,start,stop,strand,sourceLen,seq = line.strip().split()
                seq = seq.upper()
                if seqName == contigID1:
                    defLine1 = seqName + "(" + strand + ")" + "/" + start + "-" + stop
                    contigInfo1.append((defLine1,seq))
                elif seqName == contigID2:
                    defLine2 = seqName + "(" + strand + ")" + "/" + start + "-" + stop
                    contigInfo2.append((defLine2,seq))
                else:
                    print "unknown contig found!"
                    sys.exit()
            if line.startswith(' '):
                continue
    # after reading the entire file:
    if(len(contigInfo1) == len(contigInfo2)):        
        for (defLine1,seq1),(defLine2,seq2) in zip(contigInfo1,contigInfo2):
            if(len(seq1) == len(seq2)):
                #getGapRates(seq1,seq2,OUT) 
                rateList1,rateList2,totalGapRate1,totalGapRate2 = getGapRates(seq1,seq2,contigID1,contigID2,OUT1,OUT2)
            else:
                print "lengths don't match! \n"
                sys.exit()
    else:
        print "Error reading files. Lines don't match! \n"
        sys.exit()

def getGapRates(seq1,seq2,contigID1,contigID2,OUT1,OUT2):
    gapList1 = []
    gapList2 = []
    GAP1 = False
    GAP2 = False
    BLOCK = True
    if(len(seq1) == len(seq2)):
        for i in range(len(seq1)):
            if(seq1[i] == '-'):
                if BLOCK:
                    # going from BLOCK to GAP1
                    gapLength1 = 1
                    gapPos1 = i
                elif GAP2:
                    # going from GAP2 to GAP1
                    gapLength1 = 1 
                    gapPos1 = i
                else:
                    # already in GAP1, increment length
                    gapLength1 += 1
                    gapPos1 = i
                # update states
                BLOCK = False
                GAP1 = True
                GAP2 = False
            elif(seq2[i] == '-'):
                if BLOCK:
                    # going from BLOCK to GAP2
                    gapLength2 = 1
                    gapPos2 = i
                elif GAP1:
                    # going from GAP1 to GAP2
                    gapLength2 = 1
                    gapPos2 = i
                else:
                    # already in GAP2, increment length
                    gapLength2 += 1
                    gapPos2 = i
                # update states
                BLOCK = False
                GAP1 = False
                GAP2 = True
            else:
                # left a gap, into a block
                if GAP1:
                    gapStop1 = gapPos1
                    gapStart1 = abs(gapStop1 - gapLength1 + 1)
                    gapList1.append((gapLength1,gapStart1,gapStop1))
                if GAP2:
                    gapStop2 = gapPos2
                    gapStart2 = abs(gapStop2 - gapLength2 + 1)
                    gapList2.append((gapLength2,gapStart2,gapStop2))
                # update states
                BLOCK = True
                GAP1 = False
                GAP2 = False
    else:
        print "sequence lengths don't match!"
        sys.exit()

    totalGapRate1 = float(len(gapList1)) / len(seq1)
    totalGapRate1 = round(totalGapRate1,4)
    totalGapRate2 = float(len(gapList2)) / len(seq2)
    totalGapRate2 = round(totalGapRate2,4)
    for gapLength1,gapStart1,gapStop1 in gapList1:
        OUT1.write("%s\t%s\t%s\t%s\n" % (contigID1,gapLength1,gapStart1,gapStop1))
    for gapLength2,gapStart2,gapStop2 in gapList2:
        OUT2.write("%s\t%s\t%s\t%s\n" % (contigID2,gapLength2,gapStart2,gapStop2))
    return(gapList1,gapList2,totalGapRate1,totalGapRate2)


########
# MAIN #
########


usage = "Usage: " + sys.argv[0] + " <maf file> <output gap rates file basename 1> <output gap rates file basename 2>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

mafFile = sys.argv[1]
outFile1 = sys.argv[2]
outFile2 = sys.argv[3]

#print 'gapPos1,gapCount1,gapPos2,gapCount2'
#gapList1,gapList2,totalGapRate1,totalGapRate2 = parseMAF(mafFile)
parseMAF(mafFile,outFile1,outFile2)

#getGapCoords(gapList1,gapList2)
