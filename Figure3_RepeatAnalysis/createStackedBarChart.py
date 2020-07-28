import sys, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np

###############
# SUBROUTINES #
###############

#type    count   percentage
#DNA     7757346 0.292   
#LTR     2574127002      97.04   
#LINE    1457927 0.055   
#SINE    426     0.0     
#Satellite       0       0.0     
#Simple  62951315        2.373   
#Mobile Element  5369783 0.202   
#rRNA    247508  0.009   
#Other   743493  0.028   
#Total repeat content    2652654800      99.999 

def readValuesFile(valuesFile):
    with open(valuesFile,'r') as V:
        for line in V:
            if 'DNA' in line:
                # DNA     7757346 0.292
                repType,dnaCount,dnaPercentage = line.strip().split('\t')
                dnaPercentage = float(dnaPercentage)
            if 'LTR' in line:
                repType,ltrCount,ltrPercentage = line.strip().split('\t')
                ltrPercentage = float(ltrPercentage)
            if 'LINE' in line:
                repType,lineCount,linePercentage = line.strip().split('\t')
                linePercentage = float(linePercentage)
            if 'SINE' in line:
                repType,sineCount,sinePercentage = line.strip().split('\t')
                sinePercentage = float(sinePercentage)
            if 'Satellite' in line:
                repType,satelliteCount,satellitePercentage = line.strip().split('\t')
                satellitePercentage = float(satellitePercentage)
            if 'Simple' in line:
                repType,simpleCount,simplePercentage = line.strip().split('\t')
                simplePercentage = float(simplePercentage)
            if 'Mobile' in line:
                repType,mobileElementCount,mobileElementPercentage = line.strip().split('\t')
                mobileElementPercentage = float(mobileElementPercentage)
            if 'rRNA' in line:
                repType,rRNACount,rRNAPercentage = line.strip().split('\t')
                rRNAPercentage = float(rRNAPercentage)
            if 'Other' in line:
                repType,otherCount,otherPercentage = line.strip().split('\t')
                otherPercentage = float(otherPercentage)
            if 'Non-Repeat' in line:
                repType,nonRepCount,nonRepPercentage = line.strip().split('\t')
                nonRepPercentage = float(nonRepPercentage)
    return(dnaPercentage,ltrPercentage,linePercentage,sinePercentage,satellitePercentage,simplePercentage,mobileElementPercentage,rRNAPercentage,otherPercentage,nonRepPercentage)


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hop repeat file> <hop assembly length> <arabidopsis repeat file> <arabidopsis assembly length> <maize repeat file> <maize assembly length>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()

hopFile = sys.argv[1]
hopAssemblyLength = sys.argv[2]
arabidopsisFile = sys.argv[3]
arabidopsisAssemblyLength = sys.argv[4]
maizeFile = sys.argv[5]
maizeAssemblyLength = sys.argv[6]

hopAssemblyLength = int(hopAssemblyLength)
arabidopsisAssemblyLength = int(arabidopsisAssemblyLength)
maizeAssemblyLength = int(maizeAssemblyLength)

dnaPercentageH,ltrPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH,rRNAPercentageH,otherPercentageH,nonRepPercentageH = readValuesFile(hopFile)
dnaPercentageA,ltrPercentageA,linePercentageA,sinePercentageA,satellitePercentageA,simplePercentageA,mobileElementPercentageA,rRNAPercentageA,otherPercentageA,nonRepPercentageA = readValuesFile(arabidopsisFile)
dnaPercentageM,ltrPercentageM,linePercentageM,sinePercentageM,satellitePercentageM,simplePercentageM,mobileElementPercentageM,rRNAPercentageM,otherPercentageM,nonRepPercentageM = readValuesFile(maizeFile)

pieRadius = 0.8

hopLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element','Non-Repeat']
hopSizes = [ltrPercentageH,dnaPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH,nonRepPercentageH]

arabidopsisLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element','Non-Repeat']
arabidopsisSizes = [ltrPercentageA,dnaPercentageA,linePercentageA,sinePercentageA,satellitePercentageA,simplePercentageA,mobileElementPercentageA,nonRepPercentageA]

maizeLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element','Non-Repeat']
maizeSizes = [ltrPercentageM,dnaPercentageM,linePercentageM,sinePercentageM,satellitePercentageM,simplePercentageM,mobileElementPercentageM,nonRepPercentageM]

print("Hop DNA\t%s\t" % (dnaPercentageH))
print("Hop LTR\t%s\t" % (ltrPercentageH))
print("Hop LINE\t%s\t" % (linePercentageH))
print("Hop SINE\t%s\t" % (sinePercentageH))
print("Hop Satellite\t%s\t" % (satellitePercentageH))
print("Hop Simple\t%s\t" % (simplePercentageH))
print("Hop Mobile Element\t%s\t" % (mobileElementPercentageH))
print("Hop rRNA\t%s\t" % (rRNAPercentageH))
print("Hop Other\t%s\t" % (otherPercentageH))
print("Hop Non-Repeat\t%s\t" % (nonRepPercentageH))

print("Arabidopsis DNA\t%s\t" % (dnaPercentageA))
print("Arabidopsis LTR\t%s\t" % (ltrPercentageA))
print("Arabidopsis LINE\t%s\t" % (linePercentageA))
print("Arabidopsis SINE\t%s\t" % (sinePercentageA))
print("Arabidopsis Satellite\t%s\t" % (satellitePercentageA))
print("Arabidopsis Simple\t%s\t" % (simplePercentageA))
print("Arabidopsis Mobile Element\t%s\t" % (mobileElementPercentageA))
print("Arabidopsis rRNA\t%s\t" % (rRNAPercentageA))
print("Arabidopsis Other\t%s\t" % (otherPercentageA))
print("Arabidopsis Non-Repeat\t%s\t" % (nonRepPercentageA))

print("Maize DNA\t%s\t" % (dnaPercentageM))
print("Maize LTR\t%s\t" % (ltrPercentageM))
print("Maize LINE\t%s\t" % (linePercentageM))
print("Maize SINE\t%s\t" % (sinePercentageM))
print("Maize Satellite\t%s\t" % (satellitePercentageM))
print("Maize Simple\t%s\t" % (simplePercentageM))
print("Maize Mobile Element\t%s\t" % (mobileElementPercentageM))
print("Maize rRNA\t%s\t" % (rRNAPercentageM))
print("Maize Other\t%s\t" % (otherPercentageM))
print("Maize Non-Repeat\t%s\t" % (nonRepPercentageM))

################
################

# ltrPercentageH,dnaPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH
types = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element','Non-Repeat']
LTR = [ltrPercentageH,ltrPercentageA,ltrPercentageM]
DNA = [dnaPercentageH,dnaPercentageA,dnaPercentageM]
LINE = [linePercentageH,linePercentageA,linePercentageM]
SINE = [sinePercentageH,sinePercentageA,sinePercentageM]
Satellite = [satellitePercentageH,satellitePercentageA,satellitePercentageM]
Simple = [simplePercentageH,simplePercentageA,simplePercentageM]
Mobile = [mobileElementPercentageH,mobileElementPercentageA,mobileElementPercentageM]
nonRepeat = [nonRepPercentageH,nonRepPercentageA,nonRepPercentageM]
#print nonRepeat

onlyNames = ['Final Deduplicated Primary','$\it{A.}$ $\it{thaliana}$ TAIR10.1','Maize B73 RefGen_v4']
dataList = zip(LTR,DNA,LINE,SINE,Satellite,Simple,Mobile,nonRepeat,onlyNames)
#sortedDataList = sorted(dataList, key=lambda x:x[0])
ltr,dna,line,sine,satellite,simple,mobile,nonrepeat,names = zip(*dataList)
#ltr,dna,line,sine,satellite,simple,mobile,nonrepeat,names = zip(*sortedDataList)

#print ltr,dna,line,sine,satellite,simple,mobile,nonrepeat,names

plt.rcParams['axes.titlesize'] = 12
plt.rcParams['font.size'] = 10

#figWidth = 20
#figHeight = 10
#figWidth = 100
#figHeight = 30
#fig, ax1 = plt.subplots(figsize=(10,5))
#fig_size = plt.rcParams["figure.figsize"]
#fig_size[0] = figWidth
#fig_size[1] = figHeight
#plt.rcParams["figure.figsize"] = fig_size

ind = np.arange(len(names))
width = 0.5 

# colors = ['#8c510a','#01665e','#d8b365','#5ab4ac','#f6e8c3','#c7eae5','#f5f5f5']

fig, ax = plt.subplots()

p1 = plt.bar(ind, ltr, width, color='#4575b4')
p2 = plt.bar(ind, dna, width, color='#74add1', bottom=ltr)
p3 = plt.bar(ind, line, width, color='#abd9e9', bottom=[x+y for x,y in zip(ltr,dna)])
p4 = plt.bar(ind, sine, width, color='#e0f3f8', bottom=[x+y+z for x,y,z in zip(ltr,dna,line)])
p5 = plt.bar(ind, satellite, width, color='#fee090', bottom=[w+x+y+z for w,x,y,z in zip(ltr,dna,line,sine)])
p6 = plt.bar(ind, simple, width, color='#fdae61', bottom=[v+w+x+y+z for v,w,x,y,z in zip(ltr,dna,line,sine,satellite)])
p7 = plt.bar(ind, mobile, width, color='#f46d43', bottom=[u+v+w+x+y+z for u,v,w,x,y,z in zip(ltr,dna,line,sine,satellite,simple)])
p8 = plt.bar(ind, nonrepeat, width, color='#d73027', bottom=[t+u+v+w+x+y+z for t,u,v,w,x,y,z in zip(ltr,dna,line,sine,satellite,simple,mobile)])

#plt.xlim(-width,len(names)+3*width)
# plt.xlabel('Repeat percentage of assembly')
plt.ylabel('Repeat Percentage of Assembly',size=16)

#plt.set_xticklabels(['Hop primary contig assembly','$\it{Arabidopsis}$ $\it{thaliana}$ TAIR10.1','Maize B73 RefGen_v4 assembly'],minor=True)

#minor_locator = AutoMinorLocator(2)
#ax.set_xticklabels('')
#ax.xaxis.set_minor_locator(minor_locator)
#ax.set_xticklabels(['Hop primary contig assembly','$\it{Arabidopsis}$ $\it{thaliana}$ TAIR10.1','Maize B73 RefGen_v4 assembly'])  

# plt.xticks(ind + width/2.0, names, fontsize=8, ha='center')
# plt.xticks(ind,names,fontsize=8,ha='right')
plt.xticks(ind,names,fontsize=16,ha='right',rotation=30) 


plt.legend((p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p8[0]),types,frameon=False,bbox_to_anchor=[1, 1],fontsize=16)
#plt.tight_layout()
plt.savefig('stackedBarChart_v2.pdf',bbox_inches='tight')
plt.savefig('stackedBarChart_v2.svg',bbox_inches='tight')


