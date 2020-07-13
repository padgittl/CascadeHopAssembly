import os,re,sys

def readFilteredGFF(filteredGFF):
    geneDict = {}
    with open(filteredGFF,'r') as F:
        for line in F:
            if not line.startswith('#'):
                # 000000F AUGUSTUS        transcript      4099240 4100215 0.08    -       .       ID=000000F.g139.t1;Parent=000000F.g139
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if feature == 'transcript':
                    getTranscriptID = re.search('ID=(.+);Parent',attribute)
                    transcriptID = getTranscriptID.group(1)
                    geneDict[transcriptID] = 1
    return(geneDict)


def getPercentageGenesWithHints(gffFile,percentageThreshold,geneDict):
    fullPath = gffFile.strip() 
    fileName = os.path.basename(fullPath)
    filePrefix,fileExt = os.path.splitext(fileName)
    primaryID,extra = os.path.splitext(filePrefix)
    OUT = open(primaryID + '.genesSupportedByEvidence.txt','w')
    OUT.write("transcriptID\tpercentTranscriptSupportedByHints\n")
    lineList = []
    filteredList = []
    with open(gffFile,'r') as GFF:
        for line in GFF:
            # 000000F AUGUSTUS        transcript      32757   33741   0.09    -       .       ID=g1.t1;Parent=g1 
            line = line.strip()
            lineList.append(line)
    for item in lineList:
        if not item.startswith('#'):
            contigID,source,feature,start,end,score,strand,frame,attribute = item.strip().split("\t")
            if feature == 'transcript' and 'ID=' in attribute:
                getTranscriptSuffix = re.search('ID=(.+);Parent',attribute) 
                transcriptSuffix = getTranscriptSuffix.group(1) 
                transcriptID = contigID + "." + transcriptSuffix 
                filteredList.append(transcriptID)
        if '% of transcript supported by hints (any source):' in item:
            getPercentage = re.search(':\s+(.+)',item)
            percentage = getPercentage.group(1)
            percentage = percentage.strip()
            percentage = float(percentage) 
            filteredList.append(percentage)
    for i in range(1, len(filteredList), 2):
        perc = filteredList[i] 
        transcript_ID = filteredList[i-1]
        getGeneID = re.search('(.+)\.t',transcript_ID)
        geneID = getGeneID.group(1)
        #if geneID in geneDict:
        if transcript_ID in geneDict:
            if perc > percentageThreshold:
                #print perc,transcript_ID
                OUT.write("%s\t%s\n" % (transcript_ID,perc))


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gff file for single contig> <gff for entire assembly, containing longest transcripts per gene>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

gffFile = sys.argv[1]
filteredGFF = sys.argv[2]

percentageThreshold = 0

geneDict = readFilteredGFF(filteredGFF)
getPercentageGenesWithHints(gffFile,percentageThreshold,geneDict)

