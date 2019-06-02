from __future__ import division
import re,os,sys
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

def computeUngappedLengthDiff(seq1,seq2):
    seq1NoGaps = seq1.replace('-','')
    seq2NoGaps = seq2.replace('-','')
    diff = len(seq1NoGaps) - len(seq2NoGaps)
    return diff 

def extractCodonList(seq1,seq2):
    codonList = []
    primaryCodon = ""
    haplotigCodon = ""
    HAS_GAP = False
    HAS_N = False
    priFrameDisruptingGap = False
    hapFrameDisruptingGap = False
    nCountPri = 0
    nCountHap = 0
    gapCountPri = 0
    gapCountHap = 0
    priFrameDistruptList = []
    hapFrameDistruptList = []
    for j in range(len(seq1)):
        if seq1[j] == 'N':
            nCountPri += 1
        if seq1[j] == '-':
            gapCountPri += 1
        if seq2[j] == 'N':
            nCountHap += 1
        if seq2[j] == '-':
            gapCountHap += 1
    primaryCodonStartEnds = []
    primaryPos = -1
    for i in range(0,len(seq1)):
        if seq1[i] != '-':
            primaryPos += 1
        if primaryPos % 3 == 0:
        # first position of a codon
            primaryCodonStart = i
        if primaryPos % 3 == 2:
            # final position of a codon
            primaryCodonStartEnds.append((primaryCodonStart,i))
    for cStart,cEnd in primaryCodonStartEnds:
        primaryCodon = seq1[cStart:cEnd+1]
        haplotigCodon = seq2[cStart:cEnd+1]
        #print cStart,cEnd,primaryCodon,haplotigCodon
        if "-" not in primaryCodon and "-" not in haplotigCodon:
            if "N" not in primaryCodon and "N" not in haplotigCodon:
                codonList.append((primaryCodon,haplotigCodon))
            else: 
                HAS_N = True
        else:
            HAS_GAP = True
            if '-' in primaryCodon and '---' not in primaryCodon:
                priFrameDistruptList.append(primaryCodon)
                priFrameDisruptingGap = True
            if '-' in haplotigCodon and '---' not in haplotigCodon:
                hapFrameDistruptList.append(haplotigCodon)
                hapFrameDisruptingGap = True
    return(codonList,HAS_GAP,HAS_N,nCountPri,nCountHap,gapCountPri,gapCountHap,priFrameDisruptingGap,hapFrameDisruptingGap,priFrameDistruptList,hapFrameDistruptList)

def readAlignmentFastaFile(alignmentFastaFile):
    #print alignmentFastaFile
    fullFilePath = alignmentFastaFile.strip()
    alignmentFileName = os.path.basename(fullFilePath)
    baseName,fileExt = os.path.splitext(alignmentFileName)
    #print baseName
    # baseName --> 000000F_g101.t10_vs_000000F_011_g101.t10_lastzPlusStrand_NT_NT
    with open(alignmentFastaFile) as F:
        for line in F:
            if line.startswith('>'):
                #OUT = open(baseName + "_processed.fasta",'w')
                primaryDef = line.strip().replace('>','')
                primarySeq = F.next().strip()
                haplotigDef = F.next().strip().replace('>','')
                haplotigSeq = F.next().strip()
                if len(primarySeq) == len(haplotigSeq):
                    diff = computeUngappedLengthDiff(primarySeq,haplotigSeq)
                    codonList,HAS_GAP,HAS_N,nCountPri,nCountHap,gapCountPri,gapCountHap,priFrameDisruptingGap,hapFrameDisruptingGap,priFrameDistruptList,hapFrameDistruptList = extractCodonList(primarySeq,haplotigSeq)
                    primaryCodons,haplotigCodons = zip(*codonList)
                    newPrimarySeq = "".join(primaryCodons)
                    newHaplotigSeq = "".join(haplotigCodons)
                    priLenDiff = len(primarySeq) - len(newPrimarySeq)
                    hapLenDiff = len(haplotigSeq) - len(newHaplotigSeq)
                    seqOUT = open(baseName + ".fasta",'w')
                    seqOUT.write(">%s\n%s\n>%s\n%s\n" % (primaryDef,newPrimarySeq,haplotigDef,newHaplotigSeq))
                    tableOUT = open(baseName + ".txt",'w')
                    if priFrameDistruptList or hapFrameDistruptList:
                        priFrameShiftCodon = '\t'.join(priFrameDistruptList)
                        hapFrameShiftCodon = '\t'.join(hapFrameDistruptList)
                        tableOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (HAS_N,HAS_GAP,nCountPri,nCountHap,gapCountPri,gapCountHap,priLenDiff,hapLenDiff,priFrameDisruptingGap,priFrameShiftCodon,hapFrameDisruptingGap,hapFrameShiftCodon)) 
                        #tableOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (Gap,Ns,Pri_N_Count,Hap_N_Count,PriGapCount,HapGapCount,PriLenDiff,HapLenDiff,priFrameDisruptingGap,hapFrameDisruptingGap,priFrameShiftCodon,hapFrameShiftCodon))
                    else:
                        tableOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (HAS_GAP,HAS_N,nCountPri,nCountHap,gapCountPri,gapCountHap,priLenDiff,hapLenDiff,priFrameDisruptingGap,hapFrameDisruptingGap))
                else: 
                    print "error parsing file, lengths don't match!"
                    sys.exit()
            else:
                print "error parsing alignment file: " + alignmentFastaFile
                sys.exit()

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <exportAlignment output file>\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

alignmentFastaFile = sys.argv[1]

readAlignmentFastaFile(alignmentFastaFile)

