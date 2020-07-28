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
            if 'Total Repeat Content' in line:
                repType,totalRepCount,totalRepPercentage = line.strip().split('\t')
                totalRepPercentage = float(totalRepPercentage)
    return(dnaPercentage,ltrPercentage,linePercentage,sinePercentage,satellitePercentage,simplePercentage,mobileElementPercentage,rRNAPercentage,otherPercentage,totalRepPercentage)
                

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hop repeat file> <arabidopsis repeat file> <maize repeat file>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

hopFile = sys.argv[1]
arabidopsisFile = sys.argv[2]
maizeFile = sys.argv[3]

#fullFileName = arabidopsisGFF.strip()
#reducedFileName = os.path.basename(fullFileName)
#baseName,fileExt = os.path.splitext(reducedFileName)

dnaPercentageH,ltrPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH,rRNAPercentageH,otherPercentageH,totalRepPercentageH = readValuesFile(hopFile)
dnaPercentageA,ltrPercentageA,linePercentageA,sinePercentageA,satellitePercentageA,simplePercentageA,mobileElementPercentageA,rRNAPercentageA,otherPercentageA,totalRepPercentageA = readValuesFile(arabidopsisFile)
dnaPercentageM,ltrPercentageM,linePercentageM,sinePercentageM,satellitePercentageM,simplePercentageM,mobileElementPercentageM,rRNAPercentageM,otherPercentageM,totalRepPercentageM = readValuesFile(maizeFile)

pieRadius = 0.8

hopLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element']
hopSizes = [ltrPercentageH,dnaPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH]

arabidopsisLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element']
arabidopsisSizes = [ltrPercentageA,dnaPercentageA,linePercentageA,sinePercentageA,satellitePercentageA,simplePercentageA,mobileElementPercentageA]

maizeLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element']
maizeSizes = [ltrPercentageM,dnaPercentageM,linePercentageM,sinePercentageM,satellitePercentageM,simplePercentageM,mobileElementPercentageM]

print("Hop DNA\t%s\t" % (dnaPercentageH))
print("Hop LTR\t%s\t" % (ltrPercentageH))
print("Hop LINE\t%s\t" % (linePercentageH))
print("Hop SINE\t%s\t" % (sinePercentageH))
print("Hop Satellite\t%s\t" % (satellitePercentageH))
print("Hop Simple\t%s\t" % (simplePercentageH))
print("Hop Mobile Element\t%s\t" % (mobileElementPercentageH))
print("Hop rRNA\t%s\t" % (rRNAPercentageH))
print("Hop Other\t%s\t" % (otherPercentageH))

print("Arabidopsis DNA\t%s\t" % (dnaPercentageA))
print("Arabidopsis LTR\t%s\t" % (ltrPercentageA))
print("Arabidopsis LINE\t%s\t" % (linePercentageA))
print("Arabidopsis SINE\t%s\t" % (sinePercentageA))
print("Arabidopsis Satellite\t%s\t" % (satellitePercentageA))
print("Arabidopsis Simple\t%s\t" % (simplePercentageA))
print("Arabidopsis Mobile Element\t%s\t" % (mobileElementPercentageA))
print("Arabidopsis rRNA\t%s\t" % (rRNAPercentageA))
print("Arabidopsis Other\t%s\t" % (otherPercentageA))

print("Maize DNA\t%s\t" % (dnaPercentageM))
print("Maize LTR\t%s\t" % (ltrPercentageM))
print("Maize LINE\t%s\t" % (linePercentageM))
print("Maize SINE\t%s\t" % (sinePercentageM))
print("Maize Satellite\t%s\t" % (satellitePercentageM))
print("Maize Simple\t%s\t" % (simplePercentageM))
print("Maize Mobile Element\t%s\t" % (mobileElementPercentageM))
print("Maize rRNA\t%s\t" % (rRNAPercentageM))
print("Maize Other\t%s\t" % (otherPercentageM))

#colors = ['#313695','#4575b4','#74add1','#abd9e9','#ffffbf','#fee090','#fdae61','#f46d43','#d73027']
colors = ['#4575b4','#74add1','#abd9e9','#ffffbf','#fee090','#fdae61','#f46d43']
# colors = ['#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43']
#['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element']

#colors = ['#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43']
 
#plt.rcParams["figure.frameon"] = False
plt.rcParams['axes.titlesize'] = 20
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 16


#figWidth = 100
#figHeight = 30
figWidth = 20
figHeight = 10
fig, ax1 = plt.subplots(figsize=(10,5))
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = figWidth
fig_size[1] = figHeight
plt.rcParams["figure.figsize"] = fig_size

ax1 = plt.subplot(131)
ax1.set_aspect('equal')
        
hopExplode = 0
hopExplodeList = []
filteredHopSizes = []
filteredHopColors = []
hopStuff = []

for percValue,color in zip(hopSizes,colors):
    percValue = round(percValue,1)
    hopStuff.append((percValue,color))
hopStuff.sort(key=lambda x:x[0], reverse=True)
for percValue,color in hopStuff:
    if 0 < percValue and percValue < 50:
        #hopExplode += 0.05
        hopExplode += 0.075
        hopExplodeList.append(hopExplode)
        filteredHopSizes.append(percValue)
        filteredHopColors.append(color)
    else:
        hopExplode = 0.0
        hopExplodeList.append(hopExplode)
        filteredHopSizes.append(percValue)
        filteredHopColors.append(color)

# https://stackoverflow.com/questions/34035427/conditional-removal-of-labels-in-matplotlib-pie-chart
def filterAutopct(perc):
    return ('%1.1f%%' % perc) if perc > 0 else ''

hopWedges, hopTexts = ax1.pie(filteredHopSizes, colors=filteredHopColors, shadow=False, radius = pieRadius,startangle=-45, explode=hopExplodeList)

kw = dict(arrowprops=dict(arrowstyle="-"))

for i, p in enumerate(hopWedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax1.annotate(str(filteredHopSizes[i]) + "%", xy=(x, y), xytext=(0.9*np.sign(x), 1.2*y),
                 horizontalalignment=horizontalalignment,**kw)

plt.title('Hop primary contig assembly', fontsize=12)

ax2 = plt.subplot(132)
ax2.set_aspect('equal')

arabidopsisExplode = 0
arabidopsisExplodeList = []
filteredArabidopsisSizes = []
filteredArabidopsisColors = []
arabidopsisStuff = []

for percValue,color in zip(arabidopsisSizes,colors):
    percValue = round(percValue,1)
    arabidopsisStuff.append((percValue,color))
arabidopsisStuff.sort(key=lambda x:x[0], reverse=True)
for percValue,color in arabidopsisStuff:
    if 0 < percValue and percValue < 25:
        arabidopsisExplode += 0.05
        arabidopsisExplodeList.append(arabidopsisExplode)
        filteredArabidopsisSizes.append(percValue)
        filteredArabidopsisColors.append(color)
    elif 25 <= percValue:
        arabidopsisExplode = 0.0
        arabidopsisExplodeList.append(arabidopsisExplode)
        filteredArabidopsisSizes.append(percValue)
        filteredArabidopsisColors.append(color)


arabidopsisWedges, arabidopsisTexts = ax2.pie(filteredArabidopsisSizes, colors=filteredArabidopsisColors,shadow=False,radius=pieRadius,explode=arabidopsisExplodeList, startangle=-45)

kw = dict(arrowprops=dict(arrowstyle="-"))

for i, p in enumerate(arabidopsisWedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax2.annotate(str(filteredArabidopsisSizes[i]) + "%", xy=(x, y), xytext=(1.1*np.sign(x), 1.1*y),
                 horizontalalignment=horizontalalignment,**kw)

plt.title('$\it{Arabidopsis}$ $\it{thaliana}$ TAIR10.1', fontsize=12)

ax3 = plt.subplot(133)
# wedges, texts, autotexts = ax3.pie(maizeSizes, colors=colors, autopct='%1.1f%%', shadow=False, radius=1.25, startangle=90)
# plt.setp(autotexts, size=7)
ax3.set_aspect('equal')

maizeExplode = 0
maizeExplodeList = []
filteredMaizeSizes = []
filteredMaizeColors = []
maizeStuff = []

for percValue,color in zip(maizeSizes,colors):
    percValue = round(percValue,1)
    maizeStuff.append((percValue,color))
maizeStuff.sort(key=lambda x:x[0], reverse=True)
for percValue,color in maizeStuff:
    if 0 < percValue and percValue < 50:
        maizeExplode += 0.05
        maizeExplodeList.append(maizeExplode)
        filteredMaizeSizes.append(percValue)
        filteredMaizeColors.append(color)
    else:
        maizeExplode = 0.0
        maizeExplodeList.append(maizeExplode)
        filteredMaizeSizes.append(percValue)
        filteredMaizeColors.append(color)


maizeWedges, maizeTexts = ax3.pie(filteredMaizeSizes, colors=filteredMaizeColors, shadow=False, radius = pieRadius,startangle=-30, explode=maizeExplodeList)

kw = dict(arrowprops=dict(arrowstyle="-"))

for i, p in enumerate(maizeWedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax3.annotate(str(filteredMaizeSizes[i]) + "%", xy=(x, y), xytext=(1.1*np.sign(x), 1.2*y),
                 horizontalalignment=horizontalalignment,**kw)

plt.title('Maize B73 RefGen_v4 assembly', fontsize=12)
#plt.title('Total repeat content of assembly', fontsize=16)

#plt.legend([ax1, ax2, ax3],  
#           labels=hopLabels,
#           bbox_to_anchor=(1,1),
#           frameon=False,
#           prop={'size':9}
#       )

#plt.subplots_adjust(left=0.6,right=0.9,wspace=0.2, hspace=0.2)
plt.tight_layout()

#plt.title('Percentage of repeat types', fontsize=12)

plt.savefig("combinedPercentRepeatPieChart_v2.pdf")
plt.savefig("combinedPercentRepeatPieChart_v2.svg")



