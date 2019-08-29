import sys,os,re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def readUniprot(uniprotFasta):
    uniprotDict = {}
    for record in SeqIO.parse(uniprotFasta, "fasta"):
        uniprotDict[record.id] = record.description
    return(uniprotDict)

def readBlastOutputFile1(blastFile1,uniprotDict,eValueThreshold,pIdenThresh,queryCovThresh):
    bestScore ={}
    bestHits = {}

    blastDict1 = {}
    with open(blastFile1,'r') as B1:
        for line in B1:
            queryID,uniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            pIdent = float(pIdent)
            bitScore = float(bitScore)
            eValue = float(eValue)
            queryCov = int(queryCov)
            if eValue < eValueThreshold and pIdent >= pIdenThresh and queryCov >= queryCovThresh:
                if uniprotID in uniprotDict:
                    uniprotDesc = uniprotDict[uniprotID]
                    if queryID in bestHits:
                        if bitScore > bestScore[queryID]:
                            bestScore[queryID] = bitScore
                            bestHits[queryID] = (uniprotID,uniprotDesc)
                    else:
                        bestScore[queryID] = bitScore
                        bestHits[queryID] = (uniprotID,uniprotDesc)
                else:
                    "uniprotID not in uniprotDict"
                    sys.exit()
        for geneID in bestHits:
            uniprotName,uniprotDescription = bestHits[geneID]
            if 'Cannabidiolic' in uniprotDescription or 'Berberine' in uniprotDescription or 'etrahydrocannabinolic' in uniprotDescription:
                blastDict1[geneID] = uniprotName
                #print geneID,uniprotName
    return(blastDict1)

def readBlastOutputFile2(blastFile2,uniprotDict,eValueThreshold,pIdenThresh,queryCovThresh):
    bestScore ={}
    bestHits = {}

    blastDict2 = {}
    with open(blastFile2,'r') as B2:
        for line in B2:
            uniprotID,subjectID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            pIdent = float(pIdent)
            bitScore = float(bitScore)
            eValue = float(eValue)
            queryCov = int(queryCov)
            if eValue < eValueThreshold and pIdent >= pIdenThresh and queryCov >= queryCovThresh:
                if uniprotID in uniprotDict:
                    uniprotDesc = uniprotDict[uniprotID]
                    if subjectID in bestHits:
                        if bitScore > bestScore[subjectID]:
                            bestScore[subjectID] = bitScore
                            bestHits[subjectID] = (uniprotID,uniprotDesc)
                    else:
                        bestScore[subjectID] = bitScore
                        bestHits[subjectID] = (uniprotID,uniprotDesc)
                else:
                    "uniprotID not in uniprotDict"
                    sys.exit()
        for geneID in bestHits:
            uniprotName,uniprotDescription = bestHits[geneID]
            if 'Cannabidiolic' in uniprotDescription or 'Berberine' in uniprotDescription or 'etrahydrocannabinolic' in uniprotDescription:
                blastDict2[geneID] = uniprotName
                #print geneID,uniprotName
    return(blastDict2)

def readFasta(fastaFile,blastDict1,blastDict2):
    fileName = os.path.basename(fastaFile)
    baseName,fileExt = os.path.splitext(fileName)
    recordList = []
    for record in SeqIO.parse(fastaFile, "fasta"):
        if record.id in blastDict1 or record.id in blastDict2:
            recordList.append(record)
            SeqIO.write(recordList,"cannabisHits_limitedSet.fasta", "fasta") 
            
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <Cannabis pep fasta file> <blast file containing hits 1> <blast file containing hits 2> <uniprot fasta>\n"
if len(sys.argv) != 5:
    print usage
    sys.exit()

fastaFile = sys.argv[1]
blastFile1 = sys.argv[2]
blastFile2 = sys.argv[3]
uniprotFasta = sys.argv[4]

eValueThreshold = 1e-5
pIdenThresh = 50
queryCovThresh = 75

uniprotDict = readUniprot(uniprotFasta)
blastDict1 = readBlastOutputFile1(blastFile1,uniprotDict,eValueThreshold,pIdenThresh,queryCovThresh)
blastDict2 = readBlastOutputFile2(blastFile2,uniprotDict,eValueThreshold,pIdenThresh,queryCovThresh)

readFasta(fastaFile,blastDict1,blastDict2)
