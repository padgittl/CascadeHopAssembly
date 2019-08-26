#!/bin/python
import sys, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO

###############
# SUBROUTINES #
###############

def readLTRFastaFile(fastaFile):
    lengthValues = {}
    sequences = SeqIO.parse(fastaFile,'fasta')    
    for record in sequences:
        desc = record.description
        seq = str(record.seq)
        defLine,ltrType = desc.split('#')
        #print ltrType, len(seq)
        if ltrType not in lengthValues:
            lengthValues[ltrType] = []
        lengthValues[ltrType].append(len(seq))
    return lengthValues['LTR/Gypsy'],lengthValues['LTR/Copia'],lengthValues['LTR/unknown']

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <LTR non-redundant FASTA>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

LTRFastaFile = sys.argv[1]

gypsyLengths,copiaLengths,unknownLengths = readLTRFastaFile(LTRFastaFile)

lengthData = [np.asarray(gypsyLengths), np.asarray(copiaLengths), np.asarray(unknownLengths)]

# Create a figure instance
#plt.rcParams["figure.figsize"] = [8,5]
plt.rcParams["figure.figsize"] = [4,5]
#plt.rcParams["figure.figsize"] = [5,6]
#plt.rcParams["figure.figsize"] = [6,7]

# Create the boxplot

fig = plt.figure()
ax = fig.add_subplot(111)

bp = ax.boxplot(lengthData, labels = ['Gypsy','Copia','Unknown'], patch_artist=True, showfliers=False, zorder=0)
ax.tick_params(axis='x', labelsize=16) 

plt.ylabel('LTR Length (bp)',size=16)
# fill with colors
colors = ['lightpink', 'lightblue','plum']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('black')

y1 = lengthData[0]
x1 = np.random.normal(1, 0.02, len(y1))
plt.plot(x1, y1, 'r.', alpha=0.5, zorder=1)
y2 = lengthData[1]
x2 = np.random.normal(2, 0.02, len(y2))
plt.plot(x2, y2, 'b.', alpha=0.5, zorder=1)
y3 = lengthData[2]
x3 = np.random.normal(3, 0.02, len(y3))
plt.plot(x3, y3, 'm.', alpha=0.5, zorder=1)

# Save the figure
plt.savefig('LTRLengthPlot.pdf', bbox_inches='tight')
plt.savefig('LTRLengthPlot.svg', bbox_inches='tight')
