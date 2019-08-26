from __future__ import division
import re,os,sys
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

codons = {'ATT':'I',   'ATC':'I',  'ATA':'I',  'CTT':'L',  'CTC':'L',  
          'CTA':'L',  'CTG':'L',  'TTA':'L',  'TTG':'L',  'GTT':'V',  'GTC':'V',  
          'GTA':'V',  'GTG':'V',  'TTT':'F',  'TTC':'F',  'ATG':'M',  'TGT':'C',  
          'TGC':'C',  'GCT':'A',  'GCC':'A',  'GCA':'A',  'GCG':'A',  'GGT':'G',  
          'GGC':'G',  'GGA':'G',  'GGG':'G',  'CCT':'P',  'CCC':'P',  'CCA':'P',  
          'CCG':'P',  'ACT':'T',  'ACC':'T',  'ACA':'T',  'ACG':'T',  'TCT':'S',  
          'TCC':'S',  'TCA':'S',  'TCG':'S',  'AGT':'S',  'AGC':'S',  'TAT':'Y',  
          'TAC':'Y',  'TGG':'W',  'CAA':'Q',  'CAG':'Q',  'AAT':'N',  'AAC':'N',  
          'CAT':'H',  'CAC':'H',  'GAA':'E',  'GAG':'E',  'GAT':'D',  'GAC':'D',  
          'AAA':'K',  'AAG':'K',  'CGT':'R',  'CGC':'R',  'CGA':'R',  'CGG':'R',  
          'AGA':'R',  'AGG':'R',  'TAA':'*',  'TAG':'*',  'TGA':'*'}

def readHMMFile(hmmFile,eValueThresh):
    pfamDict = {}
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
                            if hopGeneID not in pfamDict:
                                pfamDict[hopGeneID] = []
                            pfamDict[hopGeneID].append((accessionID,accessionDesc,float(eValue),float(score),int(geneStart),int(geneStop)))
    return(pfamDict)

def read_dNdS_file(dNdS_file):
    dndsDict = {}
    with open(dNdS_file,'r') as D:
        for line in D:
            # 000000F_g100.t1 000000F_011_g100.t1     Plus    0       309.333333333   0.0     0       95.6666666667   0.0     0.0     0.0     --
            # return str(non_subs),str(non_sites),str(pn),str(syn_subs),str(syn_sites),str(ps),str(abs(dn)),str(abs(ds)),str(float(dn) / ds)
            primaryTranscriptID,haplotigTranscriptID,strand,non_subs,non_sites,pn,syn_subs,syn_sites,ps,dn,ds,dNdS = line.strip().split('\t')
            updatedPrimaryTranscriptID = primaryTranscriptID.replace('_g','.g')
            updatedHaplotigTranscriptID = haplotigTranscriptID.replace('_g','.g')
            if dNdS != '--':
                dNdS = float(dNdS)
            if updatedPrimaryTranscriptID not in dndsDict:
                dndsDict[updatedPrimaryTranscriptID] = []
            dndsDict[updatedPrimaryTranscriptID].append((updatedHaplotigTranscriptID,strand,dNdS))
    return(dndsDict)

def computeUngappedLengthDiff(seq1,seq2):
    seq1NoGaps = seq1.replace('-','')
    seq2NoGaps = seq2.replace('-','')
    diff = len(seq1NoGaps) - len(seq2NoGaps)
    return diff 

def identifyDomainMutation(aaPosition,domainStart,domainStop):
    status = ''
    if domainStart <= aaPosition and aaPosition <= domainStop:
        status = 'withinDomain'
    else:
        status = 'outsideDomain'
    return(status)

def extractCodonList(seq1,seq2,pfamDict,dndsDict,priTranscriptID,hapTranscriptID):
    codonList = []
    primaryCodon = ""
    haplotigCodon = ""
    HAS_GAP = False
    HAS_N = False
    priFrameDisruptingGap = False
    hapFrameDisruptingGap = False

    priFrameDistruptList = []
    hapFrameDistruptList = []

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
        # print primaryCodon,haplotigCodon
        priAminoAcid = codons[primaryCodon]
        hapAminoAcid = codons[haplotigCodon]
        #if primaryCodon == 'TAG' or primaryCodon == 'TAA' or primaryCodon == 'TGA':
        #if haplotigCodon == 'TAG' or haplotigCodon == 'TAA' or haplotigCodon == 'TGA':

        aaPos = float(cStart)/3 + 1
        # add 1 to get coordinates 1-based, which matches alignment viz coords
        cStart = cStart + 1
        cEnd = cEnd + 1
        if primaryCodon != haplotigCodon:
            if priTranscriptID in pfamDict:
                if priTranscriptID in dndsDict:
                    for hTranscriptID,strand,dNdS in dndsDict[priTranscriptID]:
                        for accessionID,accessionDesc,eValue,score,domainStart,domainStop in pfamDict[priTranscriptID]:
                            # print domainStart,domainStop
                            
                            status = identifyDomainMutation(aaPos,domainStart,domainStop)
                            if status == 'withinDomain':
                                print priTranscriptID,hapTranscriptID,cStart,cEnd,primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,aaPos,domainStart,domainStop,status,dNdS

            else:
                # print "not in pfam dict: " + priTranscriptID
                continue

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
    return(codonList,HAS_GAP,HAS_N,priFrameDisruptingGap,hapFrameDisruptingGap,priFrameDistruptList,hapFrameDistruptList)

def readAlignmentFastaFile(alignmentFastaFile,pfamDict,dndsDict):
    fullFilePath = alignmentFastaFile.strip()
    alignmentFileName = os.path.basename(fullFilePath)
    baseName,fileExt = os.path.splitext(alignmentFileName)
    getTranscriptID = re.search('(.+)_vs_(.+)_lastz',baseName)
    primaryTranscriptID = getTranscriptID.group(1)
    haplotigTranscriptID = getTranscriptID.group(2)
    updatedPrimaryTranscriptID = primaryTranscriptID.replace('_g','.g')
    updatedHaplotigTranscriptID = haplotigTranscriptID.replace('_g','.g')
    # print updatedPrimaryTranscriptID,updatedHaplotigTranscriptID
    # 000206F_g22.t1
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
                    codonList,HAS_GAP,HAS_N,priFrameDisruptingGap,hapFrameDisruptingGap,priFrameDistruptList,hapFrameDistruptList = extractCodonList(primarySeq,haplotigSeq,pfamDict,dndsDict,updatedPrimaryTranscriptID,updatedHaplotigTranscriptID)
                    primaryCodons,haplotigCodons = zip(*codonList)
                    newPrimarySeq = "".join(primaryCodons)
                    newHaplotigSeq = "".join(haplotigCodons)
                    priLenDiff = len(primarySeq) - len(newPrimarySeq)
                    hapLenDiff = len(haplotigSeq) - len(newHaplotigSeq)
                    if priFrameDistruptList or hapFrameDistruptList:
                        priFrameShiftCodon = '\t'.join(priFrameDistruptList)
                        hapFrameShiftCodon = '\t'.join(hapFrameDistruptList)
                else: 
                    print "error parsing file, lengths don't match!"
                    sys.exit()
            else:
                print "error parsing alignment file: " + alignmentFastaFile
                sys.exit()

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <processedAlignment output file> <pfam dombtlout file> <dNdS file>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

alignmentFastaFile = sys.argv[1]
hmmFile = sys.argv[2]
dNdS_file = sys.argv[3]

eValueThresh = 1e-5

pfamList = readHMMFile(hmmFile,eValueThresh)
dndsDict = read_dNdS_file(dNdS_file)

#print "priTranscriptID,hapTranscriptID,cStart,cEnd,primaryCodon,priAminoAcid,haplotigCodon,hapAminoAcid,aaPos,domainStart,domainStop,status,dNdS"
readAlignmentFastaFile(alignmentFastaFile,pfamList,dndsDict)

