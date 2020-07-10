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
                            if 'Satellite'in repeatClass:
                                SATELLITE += repeatLen
                            if "Simple_repeat" in repeatClass:
                                SIMPLE += repeatLen
                            if "Low_complexity" in repeatClass:
                                SIMPLE += repeatLen
                            if "LTR" in repeatClass:
                                LTR += repeatLen
                            if "LINE" in repeatClass:
                                LINE += repeatLen
                            if "SINE" in repeatClass:
                                SINE += repeatLen
                            if "Helitron" in repeatClass:
                                DNA += repeatLen
                            if "rRNA" in repeatClass:
                                rRNA += repeatLen
                            if "Other" in repeatClass:
                                OTHER += repeatLen
                            if repeatClass == "Unknown":
                                OTHER += repeatLen
                            if 'MobileElement' in repeatClass:
                                MOBILE_ELEMENT += repeatLen
    return(DNA,LTR,LINE,SINE,SATELLITE,SIMPLE,MOBILE_ELEMENT,rRNA,OTHER)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <arabidopsis repeat GFF file> <arabidopsis assembly length> <maize repeat GFF file> <maize assembly length>\n"
if len(sys.argv) != 5:
    print usage
    sys.exit()

arabidopsisGFF = sys.argv[1]
arabidopsisAssemblyLength = sys.argv[2]
maizeGFF = sys.argv[3]
maizeAssemblyLength = sys.argv[4]

arabidopsisAssemblyLength = int(arabidopsisAssemblyLength)
maizeAssemblyLength = int(maizeAssemblyLength)

fullFileName = arabidopsisGFF.strip()
reducedFileName = os.path.basename(fullFileName)
baseName,fileExt = os.path.splitext(reducedFileName)

aDNA,aLTR,aLINE,aSINE,aSATELLITE,aSIMPLE,aMOBILE_ELEMENT,a_rRNA,aOTHER = readGFF(arabidopsisGFF)
mDNA,mLTR,mLINE,mSINE,mSATELLITE,mSIMPLE,mMOBILE_ELEMENT,m_rRNA,mOTHER = readGFF(maizeGFF)

arabidopsisTotalRepeatCount = aDNA + aLTR + aLINE + aSINE + aSATELLITE + aSIMPLE + aMOBILE_ELEMENT + a_rRNA + aOTHER
maizeTotalRepeatCount = mDNA + mLTR + mLINE + mSINE + mSATELLITE + mSIMPLE + mMOBILE_ELEMENT + m_rRNA + mOTHER

aNonRepeatCount = arabidopsisAssemblyLength - arabidopsisTotalRepeatCount
mNonRepeatCount = maizeAssemblyLength - maizeTotalRepeatCount

perc_aDNA = round(float(aDNA) / arabidopsisAssemblyLength * 100,3)
perc_aLTR = round(float(aLTR) / arabidopsisAssemblyLength * 100,3)
perc_aLINE = round(float(aLINE) / arabidopsisAssemblyLength * 100,3)
perc_aSINE = round(float(aSINE) / arabidopsisAssemblyLength * 100,3)
perc_aSATELLITE = round(float(aSATELLITE) / arabidopsisAssemblyLength * 100,3)
perc_aSIMPLE = round(float(aSIMPLE) / arabidopsisAssemblyLength * 100,3)
perc_aMOBILE_ELEMENT = round(float(aMOBILE_ELEMENT) / arabidopsisAssemblyLength * 100,3)
perc_a_rRNA = round(float(a_rRNA) / arabidopsisAssemblyLength * 100,3)
perc_aOTHER = round(float(aOTHER) / arabidopsisAssemblyLength * 100,3)
arabidopsisTotalRepeatPercentage = perc_aDNA + perc_aLTR + perc_aLINE + perc_aSINE + perc_aSATELLITE + perc_aSIMPLE + perc_aMOBILE_ELEMENT + perc_a_rRNA + perc_aOTHER
aNonRepeatPercentage = 100 - arabidopsisTotalRepeatPercentage


perc_mDNA = round(float(mDNA) / maizeAssemblyLength * 100,3)
perc_mLTR = round(float(mLTR) / maizeAssemblyLength * 100,3)
perc_mLINE = round(float(mLINE) / maizeAssemblyLength * 100,3)
perc_mSINE = round(float(mSINE) / maizeAssemblyLength * 100,3)
perc_mSATELLITE = round(float(mSATELLITE) / maizeAssemblyLength * 100,3)
perc_mSIMPLE = round(float(mSIMPLE) / maizeAssemblyLength * 100,3)
perc_mMOBILE_ELEMENT = round(float(mMOBILE_ELEMENT) / maizeAssemblyLength * 100,3)
perc_m_rRNA = round(float(m_rRNA) / maizeAssemblyLength * 100,3)
perc_mOTHER = round(float(mOTHER) / maizeAssemblyLength * 100,3)
maizeTotalRepeatPercentage = perc_mDNA + perc_mLTR + perc_mLINE + perc_mSINE + perc_mSATELLITE + perc_mSIMPLE + perc_mMOBILE_ELEMENT + perc_m_rRNA + perc_mOTHER
mNonRepeatPercentage = 100 - maizeTotalRepeatPercentage

A = open('arabidopsis.stackedBarChartValues.txt','w')
M = open('maize.stackedBarChartValues.txt','w')

### arabidopsis
A.write("Arabidopsis\t%s\n" % (str(arabidopsisAssemblyLength)))
A.write("type\tcount\tpercentage\n")
A.write("DNA\t%s\t%s\n" % (aDNA,perc_aDNA))
A.write("LTR\t%s\t%s\n" % (aLTR,perc_aLTR))
A.write("LINE\t%s\t%s\n" % (aLINE,perc_aLINE))
A.write("SINE\t%s\t%s\n" % (aSINE,perc_aSINE))
A.write("Satellite\t%s\t%s\n" % (aSATELLITE,perc_aSATELLITE))
A.write("Simple\t%s\t%s\n" % (aSIMPLE,perc_aSIMPLE))
A.write("Mobile Element\t%s\t%s\n" % (aMOBILE_ELEMENT,perc_aMOBILE_ELEMENT))
A.write("rRNA\t%s\t%s\n" % (a_rRNA,perc_a_rRNA))
A.write("Other\t%s\t%s\n" % (aOTHER,perc_aOTHER))
A.write("Total Repeat Content\t%s\t%s\n" % (arabidopsisTotalRepeatCount,arabidopsisTotalRepeatPercentage))
A.write("Non-Repeat Content\t%s\t%s\n" % (aNonRepeatCount,aNonRepeatPercentage))


## maize
M.write("Maize\t%s\n" % (maizeAssemblyLength))
M.write("type\tcount\tpercentage\n")
M.write("DNA\t%s\t%s\n" % (mDNA,perc_mDNA))
M.write("LTR\t%s\t%s\n" % (mLTR,perc_mLTR))
M.write("LINE\t%s\t%s\n" % (mLINE,perc_mLINE))
M.write("SINE\t%s\t%s\n" % (mSINE,perc_mSINE))
M.write("Satellite\t%s\t%s\n" % (mSATELLITE,perc_mSATELLITE))
M.write("Simple\t%s\t%s\n" % (mSIMPLE,perc_mSIMPLE))
M.write("Mobile Element\t%s\t%s\n" % (mMOBILE_ELEMENT,perc_mMOBILE_ELEMENT))
M.write("rRNA\t%s\t%s\n" % (m_rRNA,perc_m_rRNA))
M.write("Other\t%s\t%s\n" % (mOTHER,perc_mOTHER))
M.write("Total Repeat Content\t%s\t%s\n" % (maizeTotalRepeatCount,maizeTotalRepeatPercentage))
M.write("Non-Repeat Content\t%s\t%s\n" % (mNonRepeatCount,mNonRepeatPercentage))

