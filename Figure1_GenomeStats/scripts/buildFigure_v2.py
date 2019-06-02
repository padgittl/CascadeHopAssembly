#!/local/cluster/bin/python
import sys,os,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

###############
# SUBROUTINES #
###############

def readBuscoOutput(buscoOutput):
    data = {}
    buscoMap = {}
    F = open(buscoOutput,'r')
    for line in F:
        if not line.startswith('#'):
            if "Missing" not in line:
                buscoId,status,contig,start,end,score,length = line.strip().split('\t')
                if status not in data:
                    data[status] = []
                data[status].append((score,length))
                buscoMap[buscoId] = status
            else:
                buscoId,status = line.strip().split('\t')
                buscoMap[buscoId] = status

    return buscoMap,data

def readBuscoFileList(buscoFileList):
    #draftPrimaryContigs,/nfs0/BB/Hendrix_Lab/Hops/BUSCO/Cascade/primary/run_cascadePrimaryBusco_embryophyta/full_table_cascadePrimaryBusco_embryophyta.tsv
    #Shinshuwase,/nfs0/Hendrix_Lab/Hops/BUSCO/PAGCascade/run_shinsuwase_embryophyta_odb9/full_table_shinsuwase_embryophyta_odb9.tsv
    #purgeHapsOnly,/nfs0/BB/Hendrix_Lab/Hops/BUSCO/Cascade/Deduplication_FINAL/run_cascadeDedupPurgeHapsOnly_embryophyta/full_table_cascadeDedupPurgeHapsOnly_embryophyta.tsv
    #deduplicatedGenome,/nfs0/BB/Hendrix_Lab/Hops/BUSCO/Cascade/Deduplication_FINAL/run_sd20_o40_mc25_embryophyta/full_table_sd20_o40_mc25_embryophyta.tsv
    #combinedMaskedCascadePrimary,/nfs0/BB/Hendrix_Lab/Hops/BUSCO/maskedCascade/run_combinedMaskedCascadePrimary_embryophyta/full_table_combinedMaskedCascadePrimary_embryophyta.tsv

    nameList = []
    onlyNames = []
    onlyFiles = []
    with open(buscoFileList,'r') as F:
        for line in F:
            name,fileName = line.strip().split(',')
            onlyFiles.append(fileName)
            if name == 'draftPrimaryContigs':
                label = 'Draft Primary'
                nameList.append((label,fileName))
                onlyNames.append(label)
            if name == 'Teamaker':
                label = 'Teamaker'
                nameList.append((label,fileName))
                onlyNames.append(label)
            if name == 'Shinshuwase':
                label = 'Shinshuwase'
                nameList.append((label,fileName))
                onlyNames.append(label)
            if name == 'purgeHapsOnly':
                label = 'Draft Primary + Purge Haplotigs'
                nameList.append((label,fileName))
                onlyNames.append(label)
            if name == 'deduplicatedGenome':
                label = 'Final Deduplicated primary'
                nameList.append((label,fileName))
                onlyNames.append(label)
            if name == 'combinedMaskedCascadePrimary':
                label = 'Masked Final Deduplicated primary'
                nameList.append((label,fileName))
                onlyNames.append(label)
    return (nameList,onlyNames,onlyFiles)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <full table file list> \n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

buscoFileList = sys.argv[1]
nameList,onlyNames,onlyFiles = readBuscoFileList(buscoFileList)

nameList = np.array(nameList)
onlyNames = np.array(onlyNames)
onlyFiles = np.array(onlyFiles)

red = '#e41a1c'
blue = '#377eb8'
green = '#4daf4a'
purple = '#decbe4'
yellow = '#ffff33'

types = ['Complete','Duplicated','Fragmented','Missing']
complete = []
fragmented = []
duplicated = []
missing = []

for n,f in zip(onlyNames,onlyFiles):
    buscoMap,data = readBuscoOutput(f)
    counts = {}
    for t in types:
        counts[t] = 0
    for buscoId in buscoMap:
        counts[buscoMap[buscoId]] += 1
    print n, counts
    complete.append(counts['Complete'])
    fragmented.append(counts['Fragmented'])
    duplicated.append(counts['Duplicated'])
    missing.append(counts['Missing'])
    dataList = zip(complete,fragmented,duplicated,missing,onlyNames)
    sortedDataList = sorted(dataList, key=lambda x:x[0])
    #print sortedDataList
    comp,frag,dup,miss,names = zip(*sortedDataList)

ind = np.arange(len(names))
width = 0.5 

plt.rcParams['axes.titlesize'] = 12
plt.rcParams['font.size'] = 12

p1 = plt.bar(ind, comp, width, color='green',align='center')
p2 = plt.bar(ind, dup, width, color='blue', bottom=comp)
p3 = plt.bar(ind, frag, width, color='yellow', bottom=[x+y for x,y in zip(comp,dup)])
p4 = plt.bar(ind, miss, width, color='purple', bottom=[x+y+z for x,y,z in zip(comp,dup,frag)])

plt.xticks(ind,names,fontsize=10,ha='right',rotation=30)

plt.title('Assess completeness of hop genome assemblies with BUSCO',size=16)
plt.ylabel('Count')

plt.legend((p1[0],p2[0],p3[0],p4[0]),types,frameon=False,bbox_to_anchor=[1, 1])

plt.savefig('hop_BUSCO_embryophyta.pdf',bbox_inches='tight')
plt.savefig('hop_BUSCO_embryophyta.svg',bbox_inches='tight')
