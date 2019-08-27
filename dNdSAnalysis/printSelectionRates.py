from __future__ import division
import re,os,sys
from Bio import SeqIO
from Bio.Seq import Seq
from dnds import dnds, substitutions, dnds_codon, dnds_codon_pair, syn_sum, translate

###############
# SUBROUTINES #
###############

def calcDNDS(alignmentFile):
    alignmentFileFullPath = alignmentFile.strip()
    fileName = os.path.basename(alignmentFileFullPath)
    baseName,fileExt = os.path.splitext(fileName)
    # 000000F_g100.t1_vs_000000F_011_g100.t1_lastzPlusStrand_NT_NT.fasta
    #getContigIDs = re.search('(.+)_vs_(\d+F_\d+|\d+F)',baseName)
    getContigIDs = re.search('(.+)_vs_(.+)_lastz',baseName)
    primaryID = getContigIDs.group(1)
    hapID = getContigIDs.group(2)
    getFullStrandInfo = re.search('_(lastz.+Strand)',baseName)
    fullStrandInfo = getFullStrandInfo.group(1)
    getLimitedStrandInfo = re.search('_lastz(.+)Strand',baseName)
    limitedStrandInfo = getLimitedStrandInfo.group(1)
    getFileNameRemainder = re.search('Strand_(.+)',baseName)
    fileNameRemainder = getFileNameRemainder.group(1)
    outName = primaryID + "_vs_" + hapID + "_" + fullStrandInfo + "_" + fileNameRemainder + "_dnds.txt"
    #print outName
    OUT = open(outName,'w')
    with open(alignmentFile,'r') as F:
        for line in F:
            line = line.strip()
            if line.startswith('>'):
                primaryDef = line.strip().replace('>','')
                primarySeq = F.next().strip() 
                #print primarySeq
                haplotigDef = F.next().strip().replace('>','')
                haplotigSeq = F.next().strip()
                #print haplotigSeq
                #haplotigDef = line.strip().replace('>','')
                #haplotigSeq = F.next().strip()
                #primaryDef = F.next().strip().replace('>','')
                #primarySeq = F.next().strip()
                primaryBase,primaryExt = primaryDef.split('_')
                getHaplotigBase = re.search('(\d+F_\d+|\d+F)_(g\d+\.t\d+)',haplotigDef)
                haplotigBase = getHaplotigBase.group(1)
                haplotigExt = getHaplotigBase.group(2)
                #print haplotigExt
                if primaryExt == haplotigExt:
                    #print primaryExt,haplotigExt 
                    if len(primarySeq) % 3 == 0:
                        if len(primarySeq) == len(haplotigSeq):
                            #print primaryExt,haplotigExt
                            rate = dnds(str(primarySeq),str(haplotigSeq))
                            #print rate
                            rate = '\t'.join(rate)
                            #print rate
                            #OUT.write(">%s\n%s\n>%s\n%s\n" % (primaryDef,primarySeq,haplotigDef,haplotigSeq))
                            OUT.write("%s\t%s\t%s\t%s\n" % (primaryDef,haplotigDef,limitedStrandInfo,rate))
                        else:
                            print "lengths not equal: " + primaryDef + "\t" + haplotigDef
                            sys.exit()
                    else:
                        print "seq not a  multiple of 3"
                        OUT.write("%s\t%s: seq not a multiple of 3" % (primaryDef,haplotigDef))
                        sys.exit()
                            
                else:
                    print "geneIDs not identical: " + primaryDef + "\t"+ haplotigDef
                    sys.exit()


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <processed alignment file from MACSE>\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

alignmentFile = sys.argv[1]

calcDNDS(alignmentFile)

