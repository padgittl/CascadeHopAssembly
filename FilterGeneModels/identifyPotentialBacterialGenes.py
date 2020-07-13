import sys,os,re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def readUniprot(uniprotFasta):
    uniprotDict = {}
    for record in SeqIO.parse(uniprotFasta, "fasta"):
        uniprotDict[record.id] = record.description
        extractUniprotID = re.search('sp\|(.+)\|',record.id)
        uniprotID = extractUniprotID.group(1)
        #print record.id
        uniprotDict[uniprotID] = record.description
    return(uniprotDict)

def parse_hop_vs_plants(hop_vs_plants,eValueThreshold):
    hop_vs_plantsDict = {}
    with open(hop_vs_plants,'r') as HP:
        for line in HP:
            queryID,uniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            hop_vs_plantsDict[queryID] = 1
    return(hop_vs_plantsDict)

def parse_plants_vs_hop(plants_vs_hop,eValueThreshold):
    plants_vs_hopDict = {}
    with open(plants_vs_hop,'r') as PH:
        for line in PH:
            uniprotID,queryID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            plants_vs_hopDict[queryID] = 1
    return(plants_vs_hopDict)

def parse_hop_vs_bacteria(hop_vs_bacteria,eValueThreshold):
    bestHits_hop_vs_bacteria = {}
    lowestEValue = {}
    bestScore = {}
    bestQueryCov = {}
    bestPercIden = {}
    with open(hop_vs_bacteria,'r') as HB:
        for line in HB:
            geneID,fullUniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            extractUniprotID = re.search('sp\|(.+)\|',fullUniprotID)
            uniprotID = extractUniprotID.group(1)
            if eValue < eValueThreshold:
                if geneID in bestHits_hop_vs_bacteria:
                    if bestScore[geneID] < bitScore:
                        bestScore[geneID] = bitScore
                        bestHits_hop_vs_bacteria[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
                else:
                    bestScore[geneID] = bitScore
                    bestHits_hop_vs_bacteria[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
    return(bestHits_hop_vs_bacteria)

def parse_bacteria_vs_hop(bacteria_vs_hop,eValueThreshold):
    bestHits_bacteria_vs_hop = {}
    lowestEValue = {}
    bestScore = {}
    bestQueryCov = {}
    bestPercIden = {}
    with open(bacteria_vs_hop,'r') as BH:
        for line in BH:
            fullUniprotID,geneID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            extractUniprotID = re.search('sp\|(.+)\|',fullUniprotID)
            uniprotID = extractUniprotID.group(1)
            if eValue < eValueThreshold:
                if geneID in bestHits_bacteria_vs_hop:
                    if bestScore[geneID] < bitScore:
                        bestScore[geneID] = bitScore
                        bestHits_bacteria_vs_hop[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
                else:
                    bestScore[geneID] = bitScore
                    bestHits_bacteria_vs_hop[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
    return(bestHits_bacteria_vs_hop)

def findBacteriaGenes(bestHits_hop_vs_bacteria,bestHits_bacteria_vs_hop,hop_vs_plantsDict,plants_vs_hopDict,uniprotDict,queryCovThresh):
    uniqueGeneDict = {}
    for hopGeneID1 in bestHits_hop_vs_bacteria:
        uniprotID1,pIdent1,eValue1,bitScore1,queryCov1 = bestHits_hop_vs_bacteria[hopGeneID1]
        if queryCov1 >= queryCovThresh:
            if hopGeneID1 not in hop_vs_plantsDict and hopGeneID1 not in plants_vs_hopDict:
                if uniprotID1 in uniprotDict:
                    uniprotDescription1 = uniprotDict[uniprotID1]
                    if hopGeneID1 not in uniqueGeneDict:
                        uniqueGeneDict[hopGeneID1] = 1
                        print("%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID1,uniprotID1,eValue1,bitScore1,queryCov1,uniprotDescription1))
                    
    for hopGeneID2 in bestHits_bacteria_vs_hop:
        uniprotID2,pIdent2,eValue2,bitScore2,queryCov2 = bestHits_bacteria_vs_hop[hopGeneID2]
        if queryCov2 >= queryCovThresh:
            if hopGeneID2 not in hop_vs_plantsDict and hopGeneID2 not in plants_vs_hopDict:
                if uniprotID2 in uniprotDict:
                    uniprotDescription2 = uniprotDict[uniprotID2]
                    if hopGeneID2 not in uniqueGeneDict:
                        uniqueGeneDict[hopGeneID2] = 1
                        print("%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID2,uniprotID2,eValue2,bitScore2,queryCov2,uniprotDescription2))


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hop vs bacteria> <bacteria vs hop> <hop vs plants> <plants vs hop> <uniprot fasta> <query coverage value>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()


hop_vs_bacteria = sys.argv[1]
bacteria_vs_hop = sys.argv[2]
hop_vs_plants = sys.argv[3]
plants_vs_hop = sys.argv[4]
uniprotFasta = sys.argv[5]
queryCoverageValue = sys.argv[6]

eValueThreshold = 1e-5
pIdenThresh = 50
queryCovThresh = int(queryCoverageValue)

uniprotDict = readUniprot(uniprotFasta)

bestHits_hop_vs_bacteria = parse_hop_vs_bacteria(hop_vs_bacteria,eValueThreshold)
bestHits_bacteria_vs_hop = parse_bacteria_vs_hop(bacteria_vs_hop,eValueThreshold)

hop_vs_plantsDict = parse_hop_vs_plants(hop_vs_plants,eValueThreshold)
plants_vs_hopDict = parse_plants_vs_hop(plants_vs_hop,eValueThreshold)

findBacteriaGenes(bestHits_hop_vs_bacteria,bestHits_bacteria_vs_hop,hop_vs_plantsDict,plants_vs_hopDict,uniprotDict,queryCovThresh)
