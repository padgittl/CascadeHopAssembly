import sys,os,re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def readRepeatAssociatedDomainFile(repeatAssociatedDomainFile):
    repeatDict = {}
    with open(repeatAssociatedDomainFile,'r') as R:
        for line in R:
            # 000271F.g23.t1  130     193     PF03732.16      Retrotransposon gag protein     4.5e-09 repeat 
            hopGeneID,geneStart,geneStop,accessionID,accessionDesc,eValue,status = line.strip().split('\t')
            if hopGeneID not in repeatDict:
                repeatDict[hopGeneID] = (accessionID,accessionDesc)
    return(repeatDict)

def readFasta(fastaFile,repeatDict):
    fileName = os.path.basename(fastaFile)
    baseName,fileExt = os.path.splitext(fileName)
    recordList = []
    for record in SeqIO.parse(fastaFile, "fasta"):
        if record.id not in repeatDict:
            recordList.append(record)
    SeqIO.write(recordList,baseName + "_repeatDomainsRemoved.fasta", "fasta") 

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <pep fasta file> <repeat-associated domain file>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

fastaFile = sys.argv[1]
repeatAssociatedDomainFile = sys.argv[2]

repeatDict = readRepeatAssociatedDomainFile(repeatAssociatedDomainFile)

readFasta(fastaFile,repeatDict)
