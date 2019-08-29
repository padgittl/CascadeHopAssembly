import sys,os,re
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#Total number of contigs: 38060 
#Assembly size is: 1346852479
#50 Contig Sum is: 673426239
#Partial Contig Sum is: 673430218, with counter: 7540, and contig length 39429

def readN50SummaryFile(n50SummaryFile):
    fullPath = n50SummaryFile.strip()
    fileName = os.path.basename(fullPath)
    fullBaseName,fileExt = os.path.splitext(fileName)
    baseName,stuff = fullBaseName.split('_N50Summary')
    
    with open(n50SummaryFile,'r') as N50:
        for line in N50:
            if "Total number of contigs" in line:
                stuff,numberContigs = line.strip().split(':')
                numberContigs = numberContigs.strip()
                contigNumberList.append(int(numberContigs))
            if "Assembly size" in line:
                stuff,assemblySize = line.strip().split(':')
                assemblySize = assemblySize.strip()
                #print assemblySize
            if "Partial Contig Sum" in line:
                stuff,maxContigLen = line.strip().split('length')
                maxContigLen = maxContigLen.strip()
                #print maxContigLen
            if "contig length" in line:
                # Partial Contig Sum is: 1881231382, with counter: 1693, and contig length 665508
                stuff,N50Value = line.strip().split('length')
                N50Value = N50Value.strip()
                #print N50Value
                n50List.append((int(N50Value),baseName))
    return(contigNumberList,n50List)
                
def scatterPlot(contigNumberList,n50List):
    for contigNumber,assemblyStuff in zip(contigNumberList,n50List):
        #11705 (2121096227, 'primaryDraft')
        #38060 (673426239, 'draftHap')
        #9214 (1910754940, 'cascadeDedupPurgeHapsOnly')
        #8836 (1880671523, 'sd20_o40_mc25')
        #37443 (870284095, 'cascadeHaplotigs')
        #292698 (906250852, 'shinshuwase')

        N50,assemblyID = assemblyStuff
        print("%s\t%s\t%s\t" % (assemblyID,contigNumber,N50))

        plt.rcParams['axes.titlesize'] = 12
        plt.rcParams['font.size'] = 10
        
        #primaryDraft_N50Summary.txt draftHap_N50Summary.txt cascadeDedupPurgeHapsOnly_N50Summary.txt cascade_v5_N50Summary_withAssemblySize.txt cascadeHaplotig_v2_N50Summary_withAssemblySize.txt shinshuwase_N50Summary.txt filteredCascadeHaplotigs_v5_N50Summary.txt teamaker_N50Summary.txt


        if assemblyID == 'primaryDraft':
            label = 'Draft Primary'
            colors = '#006d2c'
        if assemblyID == 'draftHap':
            label = 'Haplotigs'
            colors = '#41ab5d'
        if assemblyID == 'cascadeDedupPurgeHapsOnly':
            label = 'Draft Primary + Purge Haplotigs'
            colors = '#006d2c'
        #if assemblyID == 'sd20_o40_mc25':
        if assemblyID == 'cascade_v5':
            label = 'Final Deduplicated Primary'
            colors = '#006d2c'
        #if assemblyID == 'cascadeHaplotigs':
        if assemblyID == 'cascadeHaplotig_v2':
            label = 'Haplotigs + Homologous Primary Contigs'
            # label = 'Haplotigs + HPCs'
            colors = '#41ab5d'
        if assemblyID == 'filteredCascadeHaplotigs_v5':
            label = 'Haplotigs + Overlapping Haplotigs Removed'
            colors = '#41ab5d'

        if assemblyID == 'shinshuwase':
            label = 'Shinshuwase'
            colors = '#d95f0e'
        if assemblyID == 'teamaker':
            label = 'Teamaker'
            colors = '#d95f0e'
        
        # #00441b
        # colors = ['darkgreen','lightgreen','darkgreen','darkgreen','lightgreen','orange']

        plt.scatter(contigNumber,N50, color=colors, s=20)
        plt.title('Number of contigs in assembly vs N50', size=16)
        plt.xlabel('Number of contigs in assembly', size=16)
        plt.ylabel('N50',size=16)

        plt.tight_layout()

        if label == 'Shinshuwase':
            horizonAlign = 'right'
        else:
            horizonAlign = 'left'
        if label == 'Deduplicated primary assembly' or label == 'Draft primary assembly':
            vertAlign = 'top'
        else:
            vertAlign = 'bottom'

        plt.annotate(label, (contigNumber+8000,N50),horizontalalignment=horizonAlign, verticalalignment=vertAlign)

        plt.savefig("N50_vs_assemblySize.pdf", format='pdf', dpi=1000)
        plt.savefig("N50_vs_assemblySize.svg", format='svg', dpi=1000)
        

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <draft primary> <draft haps> <dedup Purge Haplotigs> <dedup primary> <dedup haps> <Shinshuwase> <filtered haplotigs> <teamaker>\n"
if len(sys.argv) != 9:
    print usage
    sys.exit()

draftPrimary = sys.argv[1]
draftHaps = sys.argv[2]
dedupPurgeHaps = sys.argv[3]
dedupPrimary = sys.argv[4]
dedupHaps = sys.argv[5]
shinshuwase = sys.argv[6]
filteredHaps = sys.argv[7]
teamaker = sys.argv[8]

contigNumberList = []
n50List = []

contigNumberList,n50List = readN50SummaryFile(draftPrimary)
contigNumberList,n50List = readN50SummaryFile(draftHaps)
contigNumberList,n50List = readN50SummaryFile(dedupPurgeHaps)
contigNumberList,n50List = readN50SummaryFile(dedupPrimary)
contigNumberList,n50List = readN50SummaryFile(dedupHaps)
contigNumberList,n50List = readN50SummaryFile(shinshuwase)
contigNumberList,n50List = readN50SummaryFile(filteredHaps)
contigNumberList,n50List = readN50SummaryFile(teamaker)

scatterPlot(contigNumberList,n50List)
