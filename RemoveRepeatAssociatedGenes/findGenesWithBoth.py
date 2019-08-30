import os,re,sys

###############
# SUBROUTINES #
###############

def readGFF(gffFile):
    coordDict = {}
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'transcript':
                    getTranscriptID = re.search('ID=(\d+F\.g\d+\.t\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in coordDict:
                        coordDict[transcriptID] = (int(start),int(end))
    return(coordDict)

def readUniProtGeneFile(uniprotGeneFile):
    uniprotDict = {}
    # Entry   Entry name      Gene names      Organism        Organism ID     Protein names
    with open(uniprotGeneFile,'r') as UF:
        for line in UF:
            if not line.startswith('Entry'):
                uniprotID,entryName,geneName,organism,organismID,proteinName = line.strip().split('\t')
                # entryName --> ABC1_ARATH // 6OMT_PAPSO // AAE12_ARATH // AB1I_ARATH // AB1K8_ARATH
                # geneName --> ABC1 AtATH10 At4g01660 T15B16.14  6OMT PSOMT2  AAE12 At1g65890 F12P19.6
                # organism --> Arabidopsis thaliana (Mouse-ear cress) // Papaver somniferum (Opium poppy)
                # organismID --> 3702 // 3469 
                # proteinName --> Protein ABC transporter 1, mitochondrial (ABC1At) (AtABC1) (EC 2.7.-.-)
                if uniprotID not in uniprotDict:
                    uniprotDict[uniprotID] = (entryName,organism,proteinName)
    return(uniprotDict)

def readBlastpOutputFile_hop_vs_uniprot(pepBlastpOutputFile1,eValueThresh):
    bestHits_hop_vs_uniprot = {}
    lowestEValue = {}
    bestScore = {}
    bestQueryCov = {}
    bestPercIden = {}
    with open(pepBlastpOutputFile1,'r') as BF1:
        for line in BF1:
            # 001620F.g33.t1  sp|O80438|MAK3_ARATH    74.444  180     46      0       17      196     11      190     6.62e-100       287     88
            geneID,fullUniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            extractUniprotID = re.search('sp\|(.+)\|',fullUniprotID)
            uniprotID = extractUniprotID.group(1)
            if eValue < eValueThresh:
                if geneID in bestHits_hop_vs_uniprot:
                    if bestScore[geneID] < bitScore:
                        bestScore[geneID] = bitScore
                        bestHits_hop_vs_uniprot[geneID] = (uniprotID,eValue,bitScore,queryCov)
                else:
                    bestScore[geneID] = bitScore
                    bestHits_hop_vs_uniprot[geneID] = (uniprotID,eValue,bitScore,queryCov)
    return(bestHits_hop_vs_uniprot)

def readBlastpOutputFile_uniprot_vs_hop(pepBlastpOutputFile2,eValueThresh):
    bestHits_uniprot_vs_hop = {}
    lowestEValue = {}
    bestScore = {}
    bestQueryCov = {}
    bestPercIden = {}
    with open(pepBlastpOutputFile2,'r') as BF2:
        for line in BF2:
            fullUniprotID,geneID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            extractUniprotID = re.search('sp\|(.+)\|',fullUniprotID)
            uniprotID = extractUniprotID.group(1)
            if eValue < eValueThresh:
                if geneID in bestHits_uniprot_vs_hop:
                    if bestScore[geneID] < bitScore:
                        bestScore[geneID] = bitScore
                        bestHits_uniprot_vs_hop[geneID] = (uniprotID,eValue,bitScore,queryCov)
                else:
                    bestScore[geneID] = bitScore
                    bestHits_uniprot_vs_hop[geneID] = (uniprotID,eValue,bitScore,queryCov)
    return(bestHits_uniprot_vs_hop)

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
                    genesWithBoth_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status)
            else:
                genesWithBoth_highestScore[hopGeneID] = score
                genesWithBoth_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status)
                # genesWithBoth_bestHits[hopGeneID] = (accessionID,accessionDesc,eValue,score,geneStart,geneStop,status)    
    return(genesWithBoth_highestScore,genesWithBoth_bestHits)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hmmscan domtblout file> <file containing pfam IDs associated with repeats> <blastp output file 1> <blastp output file 2> <uniprot geneID file> <gff file>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()

hmmFile = sys.argv[1]
repeatFile = sys.argv[2]
pepBlastpOutputFile1 = sys.argv[3]
pepBlastpOutputFile2 = sys.argv[4]
uniprotGeneFile = sys.argv[5]
gffFile = sys.argv[6]


eValueThresh = 1e-5

coordDict = readGFF(gffFile)
uniprotDict = readUniProtGeneFile(uniprotGeneFile)

bestHits_hop_vs_uniprot = readBlastpOutputFile_hop_vs_uniprot(pepBlastpOutputFile1,eValueThresh)
bestHits_uniprot_vs_hop = readBlastpOutputFile_uniprot_vs_hop(pepBlastpOutputFile2,eValueThresh)

repeatDict = readRepeatFile(repeatFile)
geneDictWithRepeat,geneDictWithoutRepeat = readHMMFile(hmmFile,repeatDict)
genesWithBoth = findGenesWithBoth(geneDictWithRepeat,geneDictWithoutRepeat)

for hopGeneID in genesWithBoth:
    for accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status in genesWithBoth[hopGeneID]:
        print hopGeneID,accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status
#    if hopGeneID in coordDict:
#        geneStart,geneStop = coordDict[hopGeneID]
#        if hopGeneID in bestHits:
#            uniprotID,eValue,bitScore,queryCov = bestHits[hopGeneID]
#            if uniprotID in uniprotDict:
#                entryName,organism,proteinName = uniprotDict[uniprotID]
#                for accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status in genesWithBoth[hopGeneID]:
#                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,geneStart,geneStop,uniprotID,entryName,organism,proteinName,eValue,bitScore,queryCov,accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status))

genesContainingOnlyRepeats_highestScore,genesContainingOnlyRepeats_bestHits = findGenesContainingOnlyRepeatDomains(geneDictWithRepeat,geneDictWithoutRepeat)

#for hopGeneID in genesContainingOnlyRepeats_bestHits:
#    accessionID,accessionDesc,eValue,score,geneStart,geneStop,status = genesContainingOnlyRepeats_bestHits[hopGeneID]
#    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,geneStart,geneStop,accessionID,accessionDesc,eValue,status))

genesWithBoth_highestScore,genesWithBoth_bestHits = getTopHitFromGenesWithBoth(genesWithBoth)
#for hopGeneID in genesWithBoth_bestHits:
#accessionID,accessionDesc,eValue,score,geneStart,geneStop,status = genesWithBoth_bestHits[hopGeneID]
#if status == "repeat":
#print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,geneStart,geneStop,accessionID,accessionDesc,eValue,status))

#for hopGeneID in genesWithBoth_bestHits: 
#    accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status = genesWithBoth_bestHits[hopGeneID]
#    print hopGeneID,accessionID,accessionDesc,eValue,score,geneStart,geneStop,alignmentProbability,status


    #if hopGeneID in coordDict:
    #    geneStart,geneStop = coordDict[hopGeneID]
    #    if hopGeneID in bestHits_hop_vs_uniprot:
    #        uniprotID1,eValue1,bitScore1,queryCov1 = bestHits_hop_vs_uniprot[hopGeneID]
    #        if uniprotID1 in uniprotDict:
    #            entryName1,organism1,proteinName1 = uniprotDict[uniprotID1]
                #for accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status in genesWithBoth[hopGeneID]:
                #print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (hopGeneID,geneStart,geneStop,uniprotID,entryName,organism,proteinName,eValue,bitScore,queryCov,accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status))
    #            continue
           # else:
                
                #print "not in uniprot dict: " + hopGeneID
    #    elif hopGeneID in bestHits_uniprot_vs_hop:
    #        uniprotID2,eValue2,bitScore2,queryCov2 = bestHits_uniprot_vs_hop[hopGeneID]
    #        if uniprotID2 in uniprotDict:
    #            entryName2,organism2,proteinName2 = uniprotDict[uniprotID2]
                #for accessionID,accessionDesc,eValue,score,domainStart,domainStop,alignmentProbability,status in genesWithBoth[hopGeneID]:
            #else:
                #print "not in uniprot dict: " + hopGeneID
                #sys.exit()
        #else:
            #print "not in best hits: " + hopGeneID
            #sys.exit()
    #else:
        #print "not in coord dict: " + hopGeneID
        #sys.exit()
