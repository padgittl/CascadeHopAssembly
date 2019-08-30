import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

# ACC   PF10417.8
# DESC  C-terminal domain of 1-Cys peroxiredoxin

def readPfamFile(pfamFile):
    repeatDict = {}
    with open(pfamFile,'r') as PFAM:
        for line in PFAM:
            if 'ACC ' in line:
                getAccessionID = re.search('\s(.+)',line)
                accessionID = getAccessionID.group(1)
                accessionID = accessionID.strip()
                fullAccessionDescription = PFAM.next().strip().split('DESC')
                accessionDesc = fullAccessionDescription[1].strip()
                # v1 don't use
                #if 'gag-' in accessionDesc or 'gag ' in accessionDesc or 'GAG-' in accessionDesc or 'GAG ' in accessionDesc or 'ransposase' in accessionDesc or 'LTR' in accessionDesc or 'Ty1' in accessionDesc or 'Ty3' in accessionDesc or 'etrotrans' in accessionDesc or 'epeat of unknown function' in accessionDesc or 'hort repeat of unknown function' in accessionDesc or 'etroviral aspartyl protease' in accessionDesc or 'everse transcriptase' in accessionDesc or 'everse transcriptase' in accessionDesc or 'lant transposon protein' in accessionDesc or 'ntegrase core domain' in accessionDesc or accessionID == "PF13952.5" or accessionID == "PF13650.5":
                # v2 use this
                if 'gag-' in accessionDesc or 'gag ' in accessionDesc or 'GAG-' in accessionDesc or 'GAG ' in accessionDesc or 'ransposase' in accessionDesc or 'LTR' in accessionDesc or 'Ty1' in accessionDesc or 'Ty3' in accessionDesc or 'etrotrans' in accessionDesc  or 'etroviral aspartyl protease' in accessionDesc or 'everse transcriptase' in accessionDesc or 'everse transcriptase' in accessionDesc or 'lant transposon protein' in accessionDesc or 'ntegrase core domain' in accessionDesc or accessionID == "PF13952.5" or "aspartyl" in accessionDesc or accessionID == "PF00075.23" or "PF03184" in accessionID or "PF13976" in accessionID or "PF00665" in accessionID or "PF07727" in accessionID or "PF05699" in accessionID or "PF14372" in accessionID or "PF03004" in accessionID or "PF05699" in accessionID or "PF14372" in accessionID or accessionID == "PF13359.5" or accessionID == "PF13917.5" or accessionID == "PF13837.5" or accessionID == "PF12776.6":
                    if accessionID != "PF00078.26":
                        if accessionID != "PF00026.22":
                            if accessionID != "PF11474.7":
                                if accessionID != "PF06550.10":
                                    if accessionID != "PF12784.6":
                                        if accessionID not in repeatDict:
                                            repeatDict[accessionID] = accessionDesc
    return(repeatDict)
            
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <pfam file> "
if len(sys.argv) != 2:
    print usage
    sys.exit()

pfamFile = sys.argv[1]

repeatDict = readPfamFile(pfamFile)
for accessionID in repeatDict:
    accessionDesc = repeatDict[accessionID]
    print("%s\t%s\t" % (accessionID,accessionDesc))
