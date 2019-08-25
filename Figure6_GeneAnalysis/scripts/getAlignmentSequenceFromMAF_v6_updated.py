import sys,re,os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Bio.AlignIO.MafIO
from Bio import AlignIO

#for a given chrom start stop, grab alignment that corresponds 
#if gap, do not update count
#gap length list

#retrieve maf region
#patch gaps in blocks with original primary sequence
#require everything to be within window for starters 
#subsequence of alignment block

def readAugustusGFF(augustusGFF):
    coordDict = {}
    with open(augustusGFF,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,stop,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'transcript':
                    getTranscriptID = re.search('ID=(.+);',attribute)
                    transcriptID = getTranscriptID.group(1)
                if feature == 'CDS':
                    #if 'g46.t2' in transcriptID:
                    #if 'g21.t2' in transcriptID:  
                    if transcriptID not in coordDict:
                        coordDict[transcriptID] = []
                    coordDict[transcriptID].append((contigID,strand,int(start),int(stop)))
                    #print transcriptID
    return(coordDict)

def readIndexFile(indexFile,coordDict):
    indexDict = {}
    with open(indexFile,'r') as IF:
        for line in IF:
            aSeqName,aStart,aStop,aStrand,aScore,posOffset = line.strip().split('\t')
            aStart = int(aStart)
            aStop = int(aStop)
            if aSeqName not in indexDict:
                indexDict[aSeqName] = []
            indexDict[aSeqName].append((int(aStart),int(aStop),aStrand,int(aScore),int(posOffset)))
        indexDict[aSeqName].sort(key=lambda x:x[0], reverse=True)
    return(indexDict)

def extractExonicMAFCoords(coordDict,indexDict,outFile):
    out = open(outFile,'w')
    for transcriptID in coordDict:
        #print transcriptID
        # initialize list to store exons for this transcriptID
        primaryExonCodingSeqList = []
        hapExonCodingSeqList = []
        transcriptStrand = ''
        # sort list of exons
        exonCount = 0
        ALL_EXONS = True

        coordDict[transcriptID].sort(key=lambda x:x[2], reverse=False)
        # now loop through each exon and extract it's alignment sequence
        #firstExonStartPos = coordDict[transcriptID][0][2]
        #lastExonStopPos = coordDict[transcriptID][-1][3]
        #print firstExonStartPos,lastExonStopPos

        for exonInfo in coordDict[transcriptID]:
            FOUND_EXON = False
            contigID,exonStrand,exonStart,exonStop = exonInfo
            #print exonStart,exonStop
            exonLen = exonStop - exonStart + 1 
            exonCount += exonLen
            #if 'g21.t2' in transcriptID:
            #if 'g46.t2' in transcriptID:
            #print contigID,exonStrand,exonStart,exonStop
            transcriptStrand = exonStrand
            if contigID in indexDict:
                # we know the contigs match up for GFF and the Index
                for seqInfo in indexDict[contigID]:
                    blockStart,blockStop,blockStrand,blockScore,posOffset = seqInfo
                    # the goal is for portion of alignment block to fall entirely within transcript coords, as opposed to simply specifying an overlapping region
                    # firstExonStartPos,lastExonStopPos
                    if blockStart <= exonStart and exonStop <= blockStop:
                        #if blockStart <= exonStart and exonStop <= blockStop:
                        #if 'g46.t2' in transcriptID:
                        FOUND_EXON = True
                        primaryExonSequence,hapExonSequence,primaryContigID,haplotigID,exonLen = extractExonSeqs(mafFile,blockStart,blockStop,blockStrand,exonStart,exonStop,posOffset)
                        #if 'g21.t2' in transcriptID:
                        #    print exonLen
                        #print transcriptID,blockStart,exonStart,exonStop,blockStop
                        if transcriptStrand == "-":
                            primaryExonSequence = revComp(primaryExonSequence)
                            hapExonSequence = revComp(hapExonSequence)
                        primaryExonCodingSeqList.append(primaryExonSequence)
                        hapExonCodingSeqList.append(hapExonSequence)
            else:
                print "Contig not found! " + contigID
                sys.exit()

            # still within specific exon "exonInfo", yet visited all alignment blocks
            if not FOUND_EXON:
                ALL_EXONS = False
            
        if ALL_EXONS:
        #if primaryExonCodingSeqList:
            if transcriptStrand == '-':
                primaryExonCodingSeqList.reverse()
                hapExonCodingSeqList.reverse()
            fullCodingExonSeqPrimary = ''.join(primaryExonCodingSeqList)
            fullCodingExonSeqHap = ''.join(hapExonCodingSeqList)
            out.write(">%s_%s\n%s\n>%s_%s\n%s\n" % (primaryContigID,transcriptID,fullCodingExonSeqPrimary,haplotigID,transcriptID,fullCodingExonSeqHap))

        # we have completed the loop over all exons; now concatenate


def extractExonSeqs(mafFile,alignStart,alignEnd,alignStrand,exonStart,exonEnd,positionOffset):
    fullPath = mafFile.strip()
    fileName = os.path.basename(fullPath)
    baseName,fileExt = os.path.splitext(fileName)
    # primary comes first in this analysis
    getContigIDs = re.search('(\d+F)_vs_(\d+F_\d+|\d+F)',baseName)
    primaryID = getContigIDs.group(1)
    hapID = getContigIDs.group(2)    
    positionOffset = int(positionOffset)
    F = open(mafFile,'r')
    F.seek(positionOffset)
    primaryLine = F.readline()
    exonCount = 0
    
    if primaryLine.startswith('s'):
        # parse the first line
        sP,seqNameP,startP,blockLenP,strandP,sourceLenP,seqP = primaryLine.strip().split()
        #print seqP
        #print primaryLine
        seqP = seqP.upper()
        primaryBlockStart = int(startP)
        primaryBlockStop = int(primaryBlockStart) + int(blockLenP)
        #print primaryBlockStart,primaryBlockStop
        haplotigLine = F.readline()        
        #print haplotigLine
        sH,seqNameH,startH,blockLenH,strandH,sourceLenH,seqH = haplotigLine.strip().split()
        seqH = seqH.upper()
        hapBlockStart = int(startH)
        hapBlockStop = int(hapBlockStart) + int(blockLenH)
        if seqNameH != hapID:
            print "error: " + seqNameH + "!=" + hapID
            sys.exit()
        if seqNameP != primaryID:
            print "error: " + seqNameP + "!=" + primaryID
            sys.exit()
            # at this point, everything is confirmed                

    exonStart = int(exonStart)
    exonEnd = int(exonEnd)
    exonLen = exonEnd - exonStart + 1
    #print exonStart,exonEnd,exonLen

    primaryPos = primaryBlockStart
    primaryExonSeq = ''
    haplotigExonSeq = ''
    if(len(seqP) == len(seqH)):
        for i in range(len(seqP)):
            if seqP[i] != '-':
                primaryPos += 1
            if exonStart <= primaryPos and primaryPos <= exonEnd:
                primaryExonSeq += seqP[i]
                #print primaryExonSeq
                haplotigExonSeq += seqH[i]
    return(primaryExonSeq,haplotigExonSeq,primaryID,hapID,exonLen)

def revComp(sequence):
    revCompSeq = ''
    complementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
    sequence = sequence[::-1]
    for base in sequence:
        base = base.upper()
        complementBase = complementDict[base]
        revCompSeq += complementBase
    #print revCompSeq
    return(revCompSeq)


########
# MAIN #
########


usage = "Usage: " + sys.argv[0] + " <maf file> <out file name> <augustus gff> <lastz maf index file>\n"
if len(sys.argv) != 5:
    print usage
    sys.exit()

mafFile = sys.argv[1]
outFile = sys.argv[2]
augustusGFF = sys.argv[3]
indexFile = sys.argv[4]

coordDict = readAugustusGFF(augustusGFF)
indexDict = readIndexFile(indexFile,coordDict)

extractExonicMAFCoords(coordDict,indexDict,outFile)
