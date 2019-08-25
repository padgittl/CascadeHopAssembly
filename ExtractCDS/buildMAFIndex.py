import sys,re,os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Bio.AlignIO.MafIO
from Bio import AlignIO

def parseMAF(mafFile):
    # 000000F_001_vs_000000F_lastzPlusStrand.txt
    # 000201F_vs_009864F_lastzMinusStrand.txt 
    # for this analysis, primary will come first! In previous iterations, haplotig came first
    fullPath = mafFile.strip()
    fileName = os.path.basename(fullPath)
    baseName,fileExt = os.path.splitext(fileName)
    outFileName = baseName + ".index"
    # 000201F_vs_009864F_lastzMinusStrand
    getContigIDs = re.search('(\d+F)_vs_(\d+F_\d+|\d+F)',baseName)
    primaryID = getContigIDs.group(1)
    hapID = getContigIDs.group(2)    
    F = open(mafFile,'r')
    OUT = open(outFileName,'w')
    line = F.readline()
    while line:
        position = F.tell() - len(line)
        #print position, line.strip()
        if line.startswith('a'):
            a,score = line.strip().split()
            getScore = re.search('score=(\d+)',score)
            score = getScore.group(1)
            score = int(score)
            #print score
        if line.startswith('s'):
            #print line
            #print position
            # parse the first line
            sP,seqNameP,startP,blockLenP,strandP,sourceLenP,seqP = line.strip().split()
            seqP = seqP.upper()
            startP = int(startP)
            stopP = int(startP) + int(blockLenP)
            # parse the next line
            line = F.readline()            
            sH,seqNameH,startH,blockLenH,strandH,sourceLenH,seqH = line.strip().split()
            seqH = seqH.upper()
            startH = int(startH)
            stopH = int(startH) + int(blockLenH)
            if seqNameH != hapID:
                print "error: " + seqNameH + "!=" + hapID
                sys.exit()
            if seqNameP != primaryID:
                print "error: " + seqNameP + "!=" + primaryID
                sys.exit()
            # at this point, everything is confirmed                
            OUT.write("%s\t%d\t%d\t%s\t%d\t%d\n" % (seqNameP,startP,stopP,strandP,score,position))
        line = F.readline()

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <maf file>\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

mafFile = sys.argv[1]

parseMAF(mafFile)



