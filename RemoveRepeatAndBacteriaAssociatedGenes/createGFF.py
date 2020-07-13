import os,re,sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#000000F AUGUSTUS        gene    32757   33741   0.09    -       .       ID=000000F.g1
#000000F AUGUSTUS        transcript      32757   33741   0.09    -       .       ID=000000F.g1.t1;Parent=000000F.g1
#000000F AUGUSTUS        exon    32757   32997   .       -       .       Parent=000000F.g1.t1
#000000F AUGUSTUS        intron  32998   33167   0.84    -       .       Parent=000000F.g1.t1
#000000F AUGUSTUS        intron  33357   33430   0.6     -       .       Parent=000000F.g1.t1
#000000F AUGUSTUS        intron  33491   33600   0.66    -       .       Parent=000000F.g1.t1
#000000F AUGUSTUS        CDS     32848   32997   0.84    -       0       ID=000000F.g1.t1.cds1;Parent=000000F.g1.t1

def readFastaFile(fastaFile):
    fastaGeneDict = {}
    for record in SeqIO.parse(fastaFile,"fasta"):
        seqRecord = record.seq
        seqID = record.id
        #getContigID_geneID = re.search('(\d+F_\d+\.g\d+|\d+F\.g\d+)',seqID)
        getGeneID = re.search('(\d+F_\d+\.g\d+.t\d+|\d+F\.g\d+.t\d+)',seqID)
        geneID = getGeneID.group(1)
        if geneID not in fastaGeneDict:
            fastaGeneDict[geneID] = True
    return(fastaGeneDict)

def readGFFFile(gffFile,fastaGeneDict):
    augustusGeneDict = {}
    transcriptDict = {}
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'gene':
                    getGeneID = re.search('ID=(.+)',attribute)
                    geneID = getGeneID.group(1)
                    if geneID not in augustusGeneDict:
                        augustusGeneDict[geneID] = (contigID,source,feature,start,end,score,strand,frame,attribute)
                if feature == 'transcript':
                    getTranscriptID = re.search('ID=(\d+F\.g\d+\.t\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in transcriptDict:
                        transcriptDict[transcriptID] = []
                    transcriptDict[transcriptID].append((contigID,source,feature,start,end,score,strand,frame,attribute))
                if feature == 'exon':
                    #000000F AUGUSTUS        exon    32757   32997   .       -       .       Parent=000000F.g1.t1  
                    getTranscriptID = re.search('Parent=(\d+F\.g\d+\.t\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in transcriptDict:
                        transcriptDict[transcriptID] = []
                    transcriptDict[transcriptID].append((contigID,source,feature,start,end,score,strand,frame,attribute))
                if feature == 'CDS':
                    # ID=000000F.g1.t1.cds1
                    getTranscriptID = re.search('ID=(\d+F\.g\d+\.t\d+).cds',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in transcriptDict:
                        transcriptDict[transcriptID] = []
                    transcriptDict[transcriptID].append((contigID,source,feature,start,end,score,strand,frame,attribute))
                if feature == 'intron':
                    # #000000F AUGUSTUS        intron  32998   33167   0.84    -       .       Parent=000000F.g1.t1 
                    getTranscriptID = re.search('Parent=(\d+F\.g\d+\.t\d+)',attribute)
                    transcriptID = getTranscriptID.group(1)
                    if transcriptID not in transcriptDict:
                        transcriptDict[transcriptID] = []
                    transcriptDict[transcriptID].append((contigID,source,feature,start,end,score,strand,frame,attribute))
    return(augustusGeneDict,transcriptDict)

def filterGFF(fastaGeneDict,augustusGeneDict,transcriptDict):
    filteredGeneList = []
    for transcriptID in fastaGeneDict:
        if transcriptID in transcriptDict:
            # 001620F.g33.t1
            getGeneID = re.search('(\d+F\.g\d+).t\d+',transcriptID)
            geneID = getGeneID.group(1)
            #print otherGFFColumns
            for otherGFFColumns in transcriptDict[transcriptID]:
                contigID,source,feature,start,end,score,strand,frame,attribute = otherGFFColumns
                filteredGeneList.append((contigID,source,feature,start,end,score,strand,frame,attribute))
                #print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (contigID,source,feature,start,end,score,strand,frame,attribute))
                #getTranscriptID = re.search('Parent=(\d+F\.g\d+\.t\d+)',attribute)
                #transcriptID = getTranscriptID.group(1)
                #print transcriptID
                #if feature == 'transcript':
                #getTranscriptID = re.search('ID=(\d+F\.g\d+\.t\d+)',attribute)
                #transcriptID = getTranscriptID.group(1)
                #print transcriptID
            if geneID in augustusGeneDict:
                #geneGFFColumns = '\t'.join(augustusGeneDict[geneID])
                #filteredGeneList.append((geneGFFColumns))
                filteredGeneList.append((augustusGeneDict[geneID]))
    return(filteredGeneList)
                
        #fullPath = fastaFile.strip() 
        #fileName = os.path.basename(fullPath)
        #filePrefix,fileExt = os.path.splitext(fileName)
        #baseName,extra = os.path.splitext(filePrefix)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gff file> <longest pep transcript fasta file>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

gffFile = sys.argv[1]
fastaFile = sys.argv[2]

fastaGeneDict = readFastaFile(fastaFile)
augustusGeneDict,transcriptDict = readGFFFile(gffFile,fastaGeneDict)
filteredGeneList = filterGFF(fastaGeneDict,augustusGeneDict,transcriptDict)
filteredGeneList.sort(key=lambda x: x[0])
for gffInfo in filteredGeneList:
    gffInfo = '\t'.join(gffInfo)
    print gffInfo
#for contigID,source,feature,start,end,score,strand,frame,attribute in filteredGeneList:
#    print contigID,source,feature,start,end,score,strand,frame,attribute
