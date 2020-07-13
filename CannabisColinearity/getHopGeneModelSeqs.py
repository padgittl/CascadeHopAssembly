import sys,os,re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#sp|A6P6W0|CASL1_CANSA
#sp|A6P6W1|CASL2_CANSA
#sp|A7IZZ1|TPS1_CANSA
#sp|A7IZZ2|TPS2_CANSA
#sp|B1Q2B6|OLIS_CANSA
#sp|C6KI62|CHSL1_CANSA
#sp|F1LKH5|PKSF3_CANSA
#sp|F1LKH6|PKSG1_CANSA
#sp|F1LKH7|PKSG2_CANSA
#sp|F1LKH8|PKSG4_CANSA
#sp|F1LKH9|PKSG5_CANSA
#sp|Q33DQ2|THCAI_CANSA
#sp|Q8RVK9|CHS_CANSA
#sp|Q95BY0|MATK_CANSA

#>sp|Q33DQ2|THCAI_CANSA Inactive tetrahydrocannabinolic acid synthase OS=Cannabis sativa PE=3 SV=1
#>sp|Q8GTB6|THCAS_CANSA Tetrahydrocannabinolic acid synthase OS=Cannabis sativa PE=1 SV=1

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
                        # store it. 
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
    bestScore = {}
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
    return(blastDict2)

def readFasta(fastaFile,blastDict1,blastDict2):
    fileName = os.path.basename(fastaFile)
    baseName,fileExt = os.path.splitext(fileName)
    recordList = []
    for record in SeqIO.parse(fastaFile, "fasta"):
        if record.id in blastDict1 or record.id in blastDict2:
            recordList.append(record)
            SeqIO.write(recordList,"hopHits_limitedSet.fasta", "fasta") 

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <augustus pep fasta> <blast file containing hits 1> <blast file containing hits 2> <uniprot fasta>\n"
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
