import sys,os,re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def readFileList(fileList):
    with open(fileList,'r') as FL:
        for fullFilePath in FL:
            fullFilePath = fullFilePath.strip()
            alignmentFileName = os.path.basename(fullFilePath)
            baseName,fileExt = os.path.splitext(alignmentFileName)
            getStrand = re.search('_(lastz.+Strand)',baseName)
            strand = getStrand.group(1)
            # 000000F_vs_000000F_001_lastzPlusStrand
            getContigIDs = re.search('(\d+F)_vs_(\d+F_\d+|\d+F)',baseName)
            primaryID = getContigIDs.group(1)
            hapID = getContigIDs.group(2)
            splitFasta(fullFilePath,strand)

def splitFasta(alignmentFile,strandInfo):
    with open(alignmentFile) as F:
        for line in F:
            if line.startswith('>'):
                primaryDef = line.strip().replace('>','')
                primarySeq = F.next().strip()
                haplotigDef = F.next().strip().replace('>','')
                haplotigSeq = F.next().strip()
                #print haplotigSeq
                primaryBase,primaryExt = primaryDef.split('_')
                getHaplotigBase = re.search('(\d+F_\d+|\d+F)_(g\d+\.t\d+)',haplotigDef)
                haplotigBase = getHaplotigBase.group(1)
                haplotigExt = getHaplotigBase.group(2)
                outName = primaryDef + "_vs_" + haplotigDef + "_" + strandInfo + ".fasta"
                OUT = open(outName,'w')
                if primaryExt == haplotigExt:
                    OUT.write(">%s\n%s\n>%s\n%s" % (primaryDef,primarySeq,haplotigDef,haplotigSeq))
                    #print primaryDef,haplotigDef
                    #print outName

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <file list> \n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

fileList = sys.argv[1]

readFileList(fileList)
