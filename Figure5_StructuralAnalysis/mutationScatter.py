import os,re,sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


#################
## SUBROUTINES ##
#################

def readConditionalProbabilityFile(conProbFile):
    conProbDict = {}
    with open(conProbFile,'r') as F:
        for line in F:
            primaryBase,hapBase,conProbValue = line.strip().split()
            conProbDict[primaryBase,hapBase] = float(conProbValue)
    return(conProbDict)

def scatterPlot(conProbDict1,conProbDict2):
    for a,b in conProbDict1:
        print a,b,conProbDict1[a,b],conProbDict2[a,b]
        # transitions
        if a == 'A' and b == 'G' or a == 'G' and b == 'A' or a == 'C' and b == 'T' or a == 'T' and b == 'C':
            plt.scatter(conProbDict1[a,b],conProbDict2[a,b], color='#fb8761')
        # tranversions
        else:
            plt.scatter(conProbDict1[a,b],conProbDict2[a,b], color='#812581')
        #C T 0.0179513113739
        #G A 0.0155610434034
        #A G 0.0127590632678
        #T C 0.0123676704441
        label = "P(" + a + "|" + b + ")"
        #plt.annotate('test', (a,b))
        if 0.01 < conProbDict2[a,b] and conProbDict2[a,b] < 0.05:
            label = "P(" + a + "|" + b + ")"
            plt.annotate(label, (conProbDict1[a,b],conProbDict2[a,b]))


    plt.xlabel('Haplotigs',size=20)
    plt.ylabel('HPCs',size=20)
    # plt.title('Conditional probabilities for haplotigs vs homologous primary contigs',fontsize=14) 
    #plt.xlim(0,1.0)
    #plt.ylim(0,1.0)
    plt.xscale('log')
    plt.yscale('log')
    xline = [10**-3,1]
    yline = [10**-3,1]
    plt.plot(xline,yline,'--')

    legendTransition = mlines.Line2D([], [], color='#fb8761', marker='o',linestyle="None",label='Transitions', markersize='10')
    legendTransversion = mlines.Line2D([], [], color='#812581', marker='o',linestyle="None",label='Transversions', markersize='10')

    plt.legend(handles=[legendTransition,legendTransversion], framealpha=None, edgecolor='inherit', frameon=False, loc='best', fontsize = '18', fancybox=False)

    # plt.legend(fancybox=False, framealpha=None, edgecolor='inherit',frameon=False)
    plt.savefig("conProb_filteredHaps.pdf", format='pdf', dpi=1000)
    plt.savefig("conProb_filteredHaps.svg", format='svg', dpi=1000)


##########
## MAIN ##
##########

if len(sys.argv) != 3 or "-h" in sys.argv or "--help" in sys.argv:
    print >> sys.stderr, "Usage: " + sys.argv[0] + " <haplotig file> <homolotig file>"
    sys.exit()

# read in the input file 
haplotigConditionalProbFile = sys.argv[1]
homolotigConditionalProbFile = sys.argv[2]

haplotigProb = readConditionalProbabilityFile(haplotigConditionalProbFile)
homolotigProb = readConditionalProbabilityFile(homolotigConditionalProbFile)

scatterPlot(haplotigProb,homolotigProb)
