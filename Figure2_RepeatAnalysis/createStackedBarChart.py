import sys, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from matplotlib.ticker import AutoMinorLocator

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

hDNA,hLTR,hLINE,hSINE,hSATELLITE,hSIMPLE,hMOBILE_ELEMENT,h_rRNA,hOTHER,hUNKNOWN,hTOTAL_REPEAT_LEN = readHopGFF(hopGFF)
aDNA,aLTR,aLINE,aSINE,aSATELLITE,aSIMPLE,aMOBILE_ELEMENT,a_rRNA,aOTHER,aUNKNOWN,aTOTAL_REPEAT_LEN = readGFF(arabidopsisGFF)
mDNA,mLTR,mLINE,mSINE,mSATELLITE,mSIMPLE,mMOBILE_ELEMENT,m_rRNA,mOTHER,mUNKNOWN,mTOTAL_REPEAT_LEN = readGFF(maizeGFF)

HOP_ASSEMBLY_LEN = int(hopAssemblyLen)
ARABIDOPSIS_ASSEMBLY_LEN = int(arabidopsisAssemblyLen)
MAIZE_ASSEMBLY_LEN = int(maizeAssemblyLen)

perc_hDNA = round(float(hDNA) / HOP_ASSEMBLY_LEN * 100,2)
perc_hLTR = round(float(hLTR) / HOP_ASSEMBLY_LEN * 100,2)
perc_hLINE = round(float(hLINE) / HOP_ASSEMBLY_LEN * 100,2)
perc_hSINE = round(float(hSINE) / HOP_ASSEMBLY_LEN * 100,2)
perc_hSATELLITE = round(float(hSATELLITE) / HOP_ASSEMBLY_LEN * 100,2)
perc_hSIMPLE = round(float(hSIMPLE) / HOP_ASSEMBLY_LEN * 100,2)
perc_hMOBILE_ELEMENT = round(float(hMOBILE_ELEMENT) / HOP_ASSEMBLY_LEN * 100,2)
perc_h_rRNA = round(float(h_rRNA) / HOP_ASSEMBLY_LEN * 100,2)
perc_hOTHER = round(float(hOTHER) / HOP_ASSEMBLY_LEN * 100,2)
perc_hUNKNOWN = round(float(hUNKNOWN) / HOP_ASSEMBLY_LEN * 100,2)
hTotalRepeat = perc_hDNA + perc_hLTR + perc_hLINE + perc_hSINE + perc_hSATELLITE + perc_hSIMPLE + perc_hMOBILE_ELEMENT
perc_hNonRepeat = 100 - hTotalRepeat

hopLabels = ['DNA','LTR','LINE','SINE','Satellite','Simple','Mobile element','Non-repeat']
hopSizes = [perc_hDNA,perc_hLTR,perc_hLINE,perc_hSINE,perc_hSATELLITE,perc_hSIMPLE,perc_hMOBILE_ELEMENT,perc_hNonRepeat]
#hopLabels = ['DNA','LTR','LINE','SINE','Satellite','Simple','Mobile element','Retro element','rRNA','Other','Unknown']
#hopSizes = [perc_hDNA,perc_hLTR,perc_hLINE,perc_hSINE,perc_hSATELLITE,perc_hSIMPLE,perc_hMOBILE_ELEMENT,perc_h_rRNA,perc_hOTHER,perc_hUNKNOWN]

perc_aDNA = round(float(aDNA) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aLTR = round(float(aLTR) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aLINE = round(float(aLINE) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aSINE = round(float(aSINE) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aSATELLITE = round(float(aSATELLITE) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aSIMPLE = round(float(aSIMPLE) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aMOBILE_ELEMENT = round(float(aMOBILE_ELEMENT) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_a_rRNA = round(float(a_rRNA) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aOTHER = round(float(aOTHER) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
perc_aUNKNOWN =round(float(aUNKNOWN) / ARABIDOPSIS_ASSEMBLY_LEN * 100,2)
aTotalRepeat = perc_aDNA + perc_aLTR + perc_aLINE + perc_aSINE + perc_aSATELLITE + perc_aSIMPLE + perc_aMOBILE_ELEMENT
perc_aNonRepeat = 100 - aTotalRepeat


arabidopsisLabels = ['DNA','LTR','LINE','SINE','Satellite','Simple','Mobile element','Non-repeat']
arabidopsisSizes = [perc_aDNA,perc_aLTR,perc_aLINE,perc_aSINE,perc_aSATELLITE,perc_aSIMPLE,perc_aMOBILE_ELEMENT,perc_aNonRepeat]
#arabidopsisLabels = ['DNA','LTR','LINE','SINE','Satellite','Simple','Mobile element','Retro element','rRNA','Other','Unknown']
#arabidopsisSizes = [perc_aDNA,perc_aLTR,perc_aLINE,perc_aSINE,perc_aSATELLITE,perc_aSIMPLE,perc_aMOBILE_ELEMENT,perc_a_rRNA,perc_aOTHER,perc_aUNKNOWN]

perc_mDNA = round(float(mDNA) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mLTR = round(float(mLTR) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mLINE = round(float(mLINE) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mSINE = round(float(mSINE) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mSATELLITE = round(float(mSATELLITE) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mSIMPLE = round(float(mSIMPLE) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mMOBILE_ELEMENT = round(float(mMOBILE_ELEMENT) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_m_rRNA = round(float(m_rRNA) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mOTHER = round(float(mOTHER) / MAIZE_ASSEMBLY_LEN * 100,2)
perc_mUNKNOWN =round(float(mUNKNOWN) / MAIZE_ASSEMBLY_LEN * 100,2)
mTotalRepeat = perc_mDNA + perc_mLTR + perc_mLINE + perc_mSINE + perc_mSATELLITE + perc_mSIMPLE + perc_mMOBILE_ELEMENT
perc_mNonRepeat = 100 - mTotalRepeat

maizeLabels = ['DNA','LTR','LINE','SINE','Satellite','Simple','Mobile element','Non-repeat']
maizeSizes = [perc_mDNA,perc_mLTR,perc_mLINE,perc_mSINE,perc_mSATELLITE,perc_mSIMPLE,perc_mMOBILE_ELEMENT,perc_mNonRepeat]
#maizeLabels = ['DNA','LTR','LINE','SINE','Satellite','Simple','Mobile element','Retro element','rRNA','Other','Unknown']
#maizeSizes = [perc_mDNA,perc_mLTR,perc_mLINE,perc_mSINE,perc_mSATELLITE,perc_mSIMPLE,perc_mMOBILE_ELEMENT,perc_m_rRNA,perc_mOTHER,perc_mUNKNOWN]

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
print("Hop non-repeat\t%s\t" % (perc_hNonRepeat))

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
print("Arabidopsis non-repeat\t%s\t" % (perc_aNonRepeat))

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
print("Maize non-repeat\t%s\t" % (perc_mNonRepeat))

types = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element','Non-repeat']
LTR = [perc_hLTR,perc_aLTR,perc_mLTR]
DNA = [perc_hDNA,perc_aDNA,perc_mDNA]
LINE = [perc_hLINE,perc_aLINE,perc_mLINE]
SINE = [perc_hSINE,perc_aSINE,perc_mSINE]
Satellite = [perc_hSATELLITE,perc_aSATELLITE,perc_mSATELLITE]
Simple = [perc_hSIMPLE,perc_aSIMPLE,perc_mSIMPLE]
Mobile = [perc_hMOBILE_ELEMENT,perc_aMOBILE_ELEMENT,perc_mMOBILE_ELEMENT]
nonRepeat = [perc_hNonRepeat,perc_aNonRepeat,perc_mNonRepeat]
#print nonRepeat

onlyNames = ['Hop primary contig assembly','$\it{Arabidopsis}$ $\it{thaliana}$ TAIR10.1','Maize B73 RefGen_v4 assembly']
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
plt.ylabel('Repeat percentage of assembly',size=16)

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
plt.savefig('stackedBarChart.pdf',bbox_inches='tight')
plt.savefig('stackedBarChart.svg',bbox_inches='tight')


