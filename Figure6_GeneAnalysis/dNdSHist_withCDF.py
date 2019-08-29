import sys,re,os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines

##############
# SUBROUTINE #
##############

def readFilteredContigMap(filteredContigMap):
    filteredContigMapDict = {}
    with open(filteredContigMap,'r') as FCM:
        for line in FCM:
            primaryID,hapID,primaryStart,primaryStop,scoreDens,source = line.strip().split('\t')
            if hapID not in filteredContigMapDict:
                filteredContigMapDict[hapID] = True
    return(filteredContigMapDict)

def readGeneFile(geneFile,filteredContigMapDict):
    HAP = []
    HPC = []
    with open(geneFile,'r') as F:
        for line in F:
            primaryGeneID,hapGeneID,strand,non_subs,non_sites,pn,syn_subs,syn_sites,ps,dN,dS,dNdS = line.strip().split('\t')
            # 000000F_011_g100.t1
            getHapID = re.search('(.+)_g',hapGeneID)
            hapID = getHapID.group(1)
            if hapID in filteredContigMapDict:
                if dNdS != "--":
                    if "-" in dNdS:
                        if dNdS == "-0.0":
                            dNdS = dNdS.replace('-','')
                            dNdS = float(dNdS)
                            if "_" in hapID:
                                HAP.append(float(dNdS))
                                #print primaryGeneID,hapGeneID,strand,non_subs,non_sites,pn,syn_subs,syn_sites,ps,dN,dS,dNdS
                            else:
                                HPC.append(float(dNdS))
                                #print primaryGeneID,hapGeneID,strand,non_subs,non_sites,pn,syn_subs,syn_sites,ps,dN,dS,dNdS
                    else:
                        if "_" in hapID:
                            HAP.append(float(dNdS))
                            #print primaryGeneID,hapGeneID,strand,non_subs,non_sites,pn,syn_subs,syn_sites,ps,dN,dS,dNdS
                        else:
                            HPC.append(float(dNdS))
                            #print primaryGeneID,hapGeneID,strand,non_subs,non_sites,pn,syn_subs,syn_sites,ps,dN,dS,dNdS
    return(HAP,HPC)

def createHist(HAP,HPC,outBaseName):

    #binWidth = 0.05
    binWidth = 0.1
    allValues = HAP + HPC

    SMALL = 8
    MEDIUM = 10
    LARGE = 12

    legendHapCDF = mlines.Line2D([], [], color='r', linestyle='solid',label='Haplotig')
    legendHPCCDF = mlines.Line2D([], [], color='b', linestyle='solid',label='HPC')

    #plt.rc('font', size=LARGE)      
    plt.rc('axes', titlesize=LARGE)  
    plt.rc('axes', labelsize=LARGE)  
    plt.rc('xtick', labelsize=LARGE) 
    plt.rc('ytick', labelsize=LARGE) 
    plt.rc('legend', fontsize=LARGE) 
    plt.rc('figure', titlesize=LARGE)

    bins = np.append(np.arange(min(allValues),max(allValues),binWidth), [np.inf])

    gap1Counts, gap1Bins, gap1Bars = plt.hist(HAP, bins = bins, histtype = 'step', color ='r', density=True, label ='Haplotig', cumulative=True, linestyle = 'solid')
    gap2Counts, gap2Bins, gap2Bars = plt.hist(HPC, bins = bins, histtype = 'step', color ='b', density=True, label= 'Homologous primary contig', alpha=0.5, cumulative=True, linestyle = 'solid')
    
    plt.xscale('log')
    plt.xlabel('dN/dS',size=16)
    plt.ylabel('CDF',size=16)
    # plt.title('CDF of dN/dS values in haplotigs vs homologous primary contigs')
    plt.legend(handles=[legendHapCDF,legendHPCCDF], framealpha=None, edgecolor='inherit', frameon=False,loc='lower right', fontsize = '10')
    plt.savefig(outBaseName + '_CDF.pdf')
    plt.savefig(outBaseName + '_CDF.svg')

########
# MAIN #
########


usage = "Usage: " + sys.argv[0] + " <gene dnds file> <filtered contig map> <out base name>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

geneFile = sys.argv[1]
filteredContigMap = sys.argv[2]
outBaseName = sys.argv[3]

filteredContigMapDict = readFilteredContigMap(filteredContigMap)
HAP,HPC = readGeneFile(geneFile,filteredContigMapDict)
createHist(HAP,HPC,outBaseName)

