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

def readHopGFF(hopGFF):
    DNA = 0
    LTR = 0
    LINE = 0
    SINE = 0
    SATELLITE = 0
    SIMPLE = 0
    MOBILE_ELEMENT = 0
    rRNA = 0
    OTHER = 0
    UNKNOWN = 0
    TOTAL_REPEAT_LEN = 0
    with open(hopGFF,'r') as HOP:
        for line in HOP:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                start = int(start)
                end = int(end)
                repeatLen = end - start
                #contigID = contigID.replace('|arrow','')
                #repType = feature
                if 'Motif' in attribute:
                    getRepeatType = re.search('Name=Motif:(.+);',attribute)
                    repType = getRepeatType.group(1)
                else:
                    getRepeatType = re.search('Name=(LTR.+);',attribute)
                    repType = getRepeatType.group(1)
                if 'DNA' in repType:
                    DNA += repeatLen
                    TOTAL_REPEAT_LEN += repeatLen
                elif 'LTR/Gypsy' in repType:
                    LTR += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'LTR/Copia' in repType:
                    LTR += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'LTR/unknown' in repType:
                    LTR += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'Retroelement' in repType:
                    LTR += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'LINE' in repType:
                    LINE += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'SINE' in repType:
                    SINE += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'MobileElement' in repType:
                    MOBILE_ELEMENT += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'rRNA' in repType:
                    rRNA += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif 'Other' in repType:
                    OTHER += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
                elif repType == "Unknown":
                    UNKNOWN += repeatLen
                    TOTAL_REPEAT_LEN += repeatLen
                elif 'Satellite' in repType:
                    SATELLITE += repeatLen
                    TOTAL_REPEAT_LEN += repeatLen
                else:
                    SIMPLE += repeatLen
                    TOTAL_REPEAT_LEN +=repeatLen
    return(DNA,LTR,LINE,SINE,SATELLITE,SIMPLE,MOBILE_ELEMENT,rRNA,OTHER,UNKNOWN,TOTAL_REPEAT_LEN)

def readGFF(gffFile):
    DNA = 0
    LTR = 0
    LINE = 0
    SINE = 0
    SATELLITE = 0
    SIMPLE = 0
    MOBILE_ELEMENT = 0
    rRNA = 0
    OTHER = 0
    UNKNOWN = 0
    TOTAL_REPEAT_LEN = 0
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                if not line.isspace():
                    if not line.startswith('score'):
                        if "position in query" not in line:
                            line = line.strip().split()
                            seqID = line[4]
                            start = line[5]
                            stop = line[6]
                            start = int(start)
                            stop = int(stop)
                            repeatLen = stop - start + 1
                            repeatMotif = line[9]
                            repeatClass = line[10]
                            if "DNA" in repeatClass:
                                DNA += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if 'Satellite'in repeatClass:
                                SATELLITE += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "Simple_repeat" in repeatClass:
                                SIMPLE += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "Low_complexity" in repeatClass:
                                SIMPLE += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "LTR" in repeatClass:
                                LTR += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "LINE" in repeatClass:
                                LINE += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "SINE" in repeatClass:
                                SINE += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "Helitron" in repeatClass:
                                DNA += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "rRNA" in repeatClass:
                                rRNA += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if "Other" in repeatClass:
                                OTHER += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if repeatClass == "Unknown":
                                UNKNOWN += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if 'MobileElement' in repeatClass:
                                MOBILE_ELEMENT += repeatLen
                                TOTAL_REPEAT_LEN +=repeatLen
    return(DNA,LTR,LINE,SINE,SATELLITE,SIMPLE,MOBILE_ELEMENT,rRNA,OTHER,UNKNOWN,TOTAL_REPEAT_LEN)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hop repeat GFF file> <hop assembly length> <Arabidopsis repeat GFF file> <Arabidopsis assembly length> <Maize repeat GFF file> <Maize assembly length>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()

hopGFF = sys.argv[1]
hopAssemblyLen = sys.argv[2]
arabidopsisGFF = sys.argv[3]
arabidopsisAssemblyLen = sys.argv[4]
maizeGFF = sys.argv[5]
maizeAssemblyLen = sys.argv[6]

fullFileName = arabidopsisGFF.strip()
reducedFileName = os.path.basename(fullFileName)
baseName,fileExt = os.path.splitext(reducedFileName)

pieRadius = 0.8

hDNA,hLTR,hLINE,hSINE,hSATELLITE,hSIMPLE,hMOBILE_ELEMENT,h_rRNA,hOTHER,hUNKNOWN,hTOTAL_REPEAT_LEN = readHopGFF(hopGFF)
aDNA,aLTR,aLINE,aSINE,aSATELLITE,aSIMPLE,aMOBILE_ELEMENT,a_rRNA,aOTHER,aUNKNOWN,aTOTAL_REPEAT_LEN = readGFF(arabidopsisGFF)
mDNA,mLTR,mLINE,mSINE,mSATELLITE,mSIMPLE,mMOBILE_ELEMENT,m_rRNA,mOTHER,mUNKNOWN,mTOTAL_REPEAT_LEN = readGFF(maizeGFF)

perc_hDNA = round(float(hDNA) / hTOTAL_REPEAT_LEN * 100,2)
perc_hLTR = round(float(hLTR) / hTOTAL_REPEAT_LEN * 100,2)
perc_hLINE = round(float(hLINE) / hTOTAL_REPEAT_LEN * 100,2)
perc_hSINE = round(float(hSINE) / hTOTAL_REPEAT_LEN * 100,2)
perc_hSATELLITE = round(float(hSATELLITE) / hTOTAL_REPEAT_LEN * 100,2)
perc_hSIMPLE = round(float(hSIMPLE) / hTOTAL_REPEAT_LEN * 100,2)
perc_hMOBILE_ELEMENT = round(float(hMOBILE_ELEMENT) / hTOTAL_REPEAT_LEN * 100,2)
perc_h_rRNA = round(float(h_rRNA) / hTOTAL_REPEAT_LEN * 100,2)
perc_hOTHER = round(float(hOTHER) / hTOTAL_REPEAT_LEN * 100,2)
perc_hUNKNOWN = round(float(hUNKNOWN) / hTOTAL_REPEAT_LEN * 100,2)

hopLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element']
hopSizes = [perc_hLTR,perc_hDNA,perc_hLINE,perc_hSINE,perc_hSATELLITE,perc_hSIMPLE,perc_hMOBILE_ELEMENT]

perc_aDNA = round(float(aDNA) / aTOTAL_REPEAT_LEN * 100,2)
perc_aLTR = round(float(aLTR) / aTOTAL_REPEAT_LEN * 100,2)
perc_aLINE = round(float(aLINE) / aTOTAL_REPEAT_LEN * 100,2)
perc_aSINE = round(float(aSINE) / aTOTAL_REPEAT_LEN * 100,2)
perc_aSATELLITE = round(float(aSATELLITE) / aTOTAL_REPEAT_LEN * 100,2)
perc_aSIMPLE = round(float(aSIMPLE) / aTOTAL_REPEAT_LEN * 100,2)
perc_aMOBILE_ELEMENT = round(float(aMOBILE_ELEMENT) / aTOTAL_REPEAT_LEN * 100,2)
perc_a_rRNA = round(float(a_rRNA) / aTOTAL_REPEAT_LEN * 100,2)
perc_aOTHER = round(float(aOTHER) / aTOTAL_REPEAT_LEN * 100,2)
perc_aUNKNOWN =round(float(aUNKNOWN) / aTOTAL_REPEAT_LEN * 100,2)

arabidopsisLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element']
arabidopsisSizes = [perc_aLTR,perc_aDNA,perc_aLINE,perc_aSINE,perc_aSATELLITE,perc_aSIMPLE,perc_aMOBILE_ELEMENT]

perc_mDNA = round(float(mDNA) / mTOTAL_REPEAT_LEN * 100,2)
perc_mLTR = round(float(mLTR) / mTOTAL_REPEAT_LEN * 100,2)
perc_mLINE = round(float(mLINE) / mTOTAL_REPEAT_LEN * 100,2)
perc_mSINE = round(float(mSINE) / mTOTAL_REPEAT_LEN * 100,2)
perc_mSATELLITE = round(float(mSATELLITE) / mTOTAL_REPEAT_LEN * 100,2)
perc_mSIMPLE = round(float(mSIMPLE) / mTOTAL_REPEAT_LEN * 100,2)
perc_mMOBILE_ELEMENT = round(float(mMOBILE_ELEMENT) / mTOTAL_REPEAT_LEN * 100,2)
perc_m_rRNA = round(float(m_rRNA) / mTOTAL_REPEAT_LEN * 100,2)
perc_mOTHER = round(float(mOTHER) / mTOTAL_REPEAT_LEN * 100,2)
perc_mUNKNOWN =round(float(mUNKNOWN) / mTOTAL_REPEAT_LEN * 100,2)

maizeLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element']
maizeSizes = [perc_mLTR,perc_mDNA,perc_mLINE,perc_mSINE,perc_mSATELLITE,perc_mSIMPLE,perc_mMOBILE_ELEMENT]

print("Hop DNA\t%s\t" % (perc_hDNA))
print("Hop LTR\t%s\t" % (perc_hLTR))
print("Hop LINE\t%s\t" % (perc_hLINE))
print("Hop SINE\t%s\t" % (perc_hSINE))
print("Hop Satellite\t%s\t" % (perc_hSATELLITE))
print("Hop Simple\t%s\t" % (perc_hSIMPLE))
print("Hop Mobile Element\t%s\t" % (perc_hMOBILE_ELEMENT))
print("Hop rRNA\t%s\t" % (perc_h_rRNA))
print("Hop Other\t%s\t" % (perc_hOTHER))
print("Hop Unknown\t%s\t" % (perc_hUNKNOWN))

print("Arabidopsis DNA\t%s\t" % (perc_aDNA))
print("Arabidopsis LTR\t%s\t" % (perc_aLTR))
print("Arabidopsis LINE\t%s\t" % (perc_aLINE))
print("Arabidopsis SINE\t%s\t" % (perc_aSINE))
print("Arabidopsis Satellite\t%s\t" % (perc_aSATELLITE))
print("Arabidopsis Simple\t%s\t" % (perc_aSIMPLE))
print("Arabidopsis Mobile Element\t%s\t" % (perc_aMOBILE_ELEMENT))
print("Arabidopsis rRNA\t%s\t" % (perc_a_rRNA))
print("Arabidopsis Other\t%s\t" % (perc_aOTHER))
print("Arabidopsis Unknown\t%s\t" % (perc_aUNKNOWN))

print("Maize mDNA\t%s\t" % (perc_mDNA))
print("Maize mLTR\t%s\t" % (perc_mLTR))
print("Maize mLINE\t%s\t" % (perc_mLINE))
print("Maize mSINE\t%s\t" % (perc_mSINE))
print("Maize mSATELLITE\t%s\t" % (perc_mSATELLITE))
print("Maize mSIMPLE\t%s\t" % (perc_mSIMPLE))
print("Maize mMOBILE_ELEMENT\t%s\t" % (perc_mMOBILE_ELEMENT))
print("Maize m_rRNA\t%s\t" % (perc_m_rRNA))
print("Maize mOTHER\t%s\t" % (perc_mOTHER))
print("Maize mUNKNOWN\t%s\t" % (perc_mUNKNOWN))

# colors = ['#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43']
# ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element']                                                                                                                             

colors = ['#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43']
 
#plt.rcParams["figure.frameon"] = False
plt.rcParams['axes.titlesize'] = 20
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 12


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
for percValue,color in zip(hopSizes,colors):
    if 0 < percValue and percValue < 50:
        #hopExplode += 0.05
        hopExplode += 0.075
        percValue = round(percValue,1)
        hopExplodeList.append(hopExplode)
        filteredHopSizes.append(percValue)
        filteredHopColors.append(color)
    else:
        hopExplode = 0.0
        percValue = round(percValue,1)
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

plt.title('Hop primary contig assembly', fontsize=10)

ax2 = plt.subplot(132)
ax2.set_aspect('equal')

arabidopsisExplode = 0
arabidopsisExplodeList = []
filteredArabidopsisSizes = []
filteredArabidopsisColors = []
for percValue,color in zip(arabidopsisSizes,colors):
    if 0 < percValue and percValue < 25:
        percValue = round(percValue,1)
        arabidopsisExplode += 0.05
        arabidopsisExplodeList.append(arabidopsisExplode)
        filteredArabidopsisSizes.append(percValue)
        filteredArabidopsisColors.append(color)
    elif 25 < percValue:
        percValue = round(percValue,1)
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

plt.title('$\it{Arabidopsis}$ $\it{thaliana}$ TAIR10.1', fontsize=10)

ax3 = plt.subplot(133)
# wedges, texts, autotexts = ax3.pie(maizeSizes, colors=colors, autopct='%1.1f%%', shadow=False, radius=1.25, startangle=90)
# plt.setp(autotexts, size=7)
ax3.set_aspect('equal')

maizeExplode = 0
maizeExplodeList = []
filteredMaizeSizes = []
filteredMaizeColors = []
for percValue,color in zip(maizeSizes,colors):
    if 0 < percValue and percValue < 50:
        percValue = round(percValue,1)
        maizeExplode += 0.05
        maizeExplodeList.append(maizeExplode)
        filteredMaizeSizes.append(percValue)
        filteredMaizeColors.append(color)
    else:
        percValue = round(percValue,1)
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

plt.title('Maize B73 RefGen_v4 assembly', fontsize=10)
plt.title('Total repeat content of assembly', fontsize=16)

#plt.legend([ax1, ax2, ax3],  
#           labels=hopLabels,
#           bbox_to_anchor=(1,1),
#           frameon=False,
#           prop={'size':9}
#       )

#plt.subplots_adjust(left=0.6,right=0.9,wspace=0.2, hspace=0.2)
plt.tight_layout()

#plt.title('Percentage of repeat types', fontsize=12)

plt.savefig("combinedPercentRepeatPieChart.pdf")
plt.savefig("combinedPercentRepeatPieChart.svg")



