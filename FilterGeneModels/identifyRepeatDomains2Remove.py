import os,re,sys

###############
# SUBROUTINES #
###############

def readHMMFile(hmmFile,repeatDict):
    geneDictWithRepeat = {}
    geneDictWithoutRepeat = {}
    with open(hmmFile,'r') as HMM:
        for line in HMM:
            if not line.startswith('#'):
                if not line.startswith('-'):
                    if not line.isspace():
                        domainInfo = line.strip().split()[0:22]
                        accessionID = domainInfo[1]
                        hopGeneID = domainInfo[3]
                        eValue = domainInfo[6]
                        score = domainInfo[7]
                        geneStart = domainInfo[17]
                        geneStop = domainInfo[18]
                        alignmentProbability = domainInfo[21]
                        accessionDesc = line.strip().split()[22:]
                        accessionDesc = ' '.join(accessionDesc)
                        eValue = float(eValue)
                        if eValue < eValueThresh:
                            if accessionID in repeatDict:
                                if hopGeneID not in geneDictWithRepeat:
                                    geneDictWithRepeat[hopGeneID] = []
                                geneDictWithRepeat[hopGeneID].append((accessionID,accessionDesc,float(eValue),float(score),int(geneStart),int(geneStop),float(alignmentProbability)))
                            else:
                                if hopGeneID not in geneDictWithoutRepeat:
                                    geneDictWithoutRepeat[hopGeneID] = []
                                geneDictWithoutRepeat[hopGeneID].append((accessionID,accessionDesc,float(eValue),float(score),int(geneStart),int(geneStop),float(alignmentProbability)))
    return(geneDictWithRepeat,geneDictWithoutRepeat)

def readRepeatFile(repeatFile):
    repeatDict = {}
    with open(repeatFile,'r') as RF:
        for line in RF:
            accessionID,accessionDesc = line.strip().split('\t')
            if accessionID not in repeatDict:
                repeatDict[accessionID] = accessionDesc
    return(repeatDict)

def findGenesContainingOnlyRepeatDomains(geneDictWithRepeat,geneDictWithoutRepeat):
    genesContainingOnlyRepeats_highestScore = {}
    genesContainingOnlyRepeats_bestHits = {}
    status = 'repeat'
    for hopGeneID in geneDictWithRepeat:
        if hopGeneID not in geneDictWithoutRepeat:
            for accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability in geneDictWithRepeat[hopGeneID]:
                eValue = float(eValue)
                score = int(score)
                if hopGeneID in genesContainingOnlyRepeats_bestHits:
                    if genesContainingOnlyRepeats_highestScore[hopGeneID] < score:
                        genesContainingOnlyRepeats_highestScore[hopGeneID] = score
                        genesContainingOnlyRepeats_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,status)
                else:
                    genesContainingOnlyRepeats_highestScore[hopGeneID] = score
                    genesContainingOnlyRepeats_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,status)
    return(genesContainingOnlyRepeats_highestScore,genesContainingOnlyRepeats_bestHits)

def findGenesWithBoth(geneDictWithRepeat,geneDictWithoutRepeat):
    genesWithBoth = {}
    for hopGeneID in geneDictWithoutRepeat:
        if hopGeneID in geneDictWithRepeat:
            # hopGeneID has both
            if hopGeneID not in genesWithBoth:
                genesWithBoth[hopGeneID] = []
            for accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability in geneDictWithoutRepeat[hopGeneID]:
                genesWithBoth[hopGeneID].append((accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,"noRepeat"))
            for rAccessionID,rAccessionDesc,r_eValue,rScore,rGeneStart,rGeneStop,rAlignmentProbability in geneDictWithRepeat[hopGeneID]:
                genesWithBoth[hopGeneID].append((rAccessionID,rAccessionDesc,r_eValue,rScore,rGeneStart,rGeneStop,rAlignmentProbability,"repeat"))
            genesWithBoth[hopGeneID].sort(key=lambda x: x[4], reverse=True)
    return(genesWithBoth)

def getTopHitFromGenesWithBoth(genesWithBoth):
    genesWithBoth_bestHits = {}
    genesWithBoth_highestScore = {}
    for hopGeneID in genesWithBoth:
        for accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status in genesWithBoth[hopGeneID]:
            eValue = float(eValue)
            if hopGeneID in genesWithBoth_bestHits:
                if genesWithBoth_highestScore[hopGeneID] < score:
                    genesWithBoth_highestScore[hopGeneID] = score
                    genesWithBoth_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,status)
            else:
                genesWithBoth_highestScore[hopGeneID] = score
                genesWithBoth_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,status)    
    return(genesWithBoth_highestScore,genesWithBoth_bestHits)
                

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hmmscan domtblout file> <file containing pfam IDs associated with repeats>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

hmmFile = sys.argv[1]
repeatFile = sys.argv[2]

eValueThresh = 1e-5

repeatDict = readRepeatFile(repeatFile)
geneDictWithRepeat,geneDictWithoutRepeat = readHMMFile(hmmFile,repeatDict)
genesWithBoth = findGenesWithBoth(geneDictWithRepeat,geneDictWithoutRepeat)

#for hopGeneID in genesWithBoth:
#    for accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status in genesWithBoth[hopGeneID]:
#        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status))

genesContainingOnlyRepeats_highestScore,genesContainingOnlyRepeats_bestHits = findGenesContainingOnlyRepeatDomains(geneDictWithRepeat,geneDictWithoutRepeat)

for hopGeneID in genesContainingOnlyRepeats_bestHits:
    accessionID,accessionDesc,eValue,score,geneStart,geneStop,status = genesContainingOnlyRepeats_bestHits[hopGeneID]
    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,geneStart,geneStop,accessionID,accessionDesc,eValue,status))
    
genesWithBoth_highestScore,genesWithBoth_bestHits = getTopHitFromGenesWithBoth(genesWithBoth)
for hopGeneID in genesWithBoth_bestHits:
    accessionID,accessionDesc,eValue,score,geneStart,geneStop,status = genesWithBoth_bestHits[hopGeneID]
    #if status == "repeat":
    #    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,geneStart,geneStop,accessionID,accessionDesc,eValue,status))

