import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############
        
def readFileList(haplotigFileList):
    with open(haplotigFileList,'r') as FL:
        for fileName in FL:
            fileName = fileName.strip()
            prepareHaplotigLASTZCommands(fileName)

def prepareHaplotigLASTZCommands(haplotigFastaFile):
    lastzVersion = "lastz"
    pathToPrimaryContigs = "path2FastaFiles/"
    contigDict = {}
    haplotigFastaBaseName = os.path.basename(haplotigFastaFile)
    haplotigID,fileExt = os.path.splitext(haplotigFastaBaseName)
    getPrimaryID = re.search('(\d+F)_',haplotigID)
    primaryID = getPrimaryID.group(1)
    if haplotigID not in contigDict:
        contigDict[haplotigID] = {}
    if primaryID not in contigDict[haplotigID]:
        contigDict[haplotigID][primaryID] = 1
        
        print lastzVersion + " " + haplotigFastaFile + " " + pathToPrimaryContigs + primaryID + ".fa --gfextend --hspthresh=20000 --chain --gapped --output=LASTZFiles/" + haplotigID + "_vs_" + primaryID + "_lastzPlusStrand.txt --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,continuity,coverage --inner=10000 --identity=80 --strand=plus"
        print lastzVersion + " " + haplotigFastaFile + " " + pathToPrimaryContigs + primaryID + ".fa --gfextend --hspthresh=20000 --chain --gapped --output=LASTZFiles/" + haplotigID + "_vs_" + primaryID + "_lastzMinusStrand.txt --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,continuity,coverage --inner=10000 --identity=80 --strand=minus"

    
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <haplotig file list>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

haplotigFileList = sys.argv[1]

readFileList(haplotigFileList)

