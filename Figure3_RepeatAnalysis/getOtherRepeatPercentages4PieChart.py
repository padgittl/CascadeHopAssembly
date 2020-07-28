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
                                OTHER += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
                            if 'MobileElement' in repeatClass:
                                MOBILE_ELEMENT += repeatLen
                                TOTAL_REPEAT_LEN += repeatLen
    return(DNA,LTR,LINE,SINE,SATELLITE,SIMPLE,MOBILE_ELEMENT,rRNA,OTHER,TOTAL_REPEAT_LEN)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <Arabidopsis repeat GFF file> <Maize repeat GFF file> \n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

arabidopsisGFF = sys.argv[1]
maizeGFF = sys.argv[2]

fullFileName = arabidopsisGFF.strip()
reducedFileName = os.path.basename(fullFileName)
baseName,fileExt = os.path.splitext(reducedFileName)

aDNA,aLTR,aLINE,aSINE,aSATELLITE,aSIMPLE,aMOBILE_ELEMENT,a_rRNA,aOTHER,aTOTAL_REPEAT_LEN = readGFF(arabidopsisGFF)
mDNA,mLTR,mLINE,mSINE,mSATELLITE,mSIMPLE,mMOBILE_ELEMENT,m_rRNA,mOTHER,mTOTAL_REPEAT_LEN = readGFF(maizeGFF)

#arabidopsisTotalRepeatCount = aDNA + aLTR + aLINE + aSINE + aSATELLITE + aSIMPLE + aMOBILE_ELEMENT + a_rRNA + aOTHER
#maizeTotalRepeatCount = mDNA + mLTR + mLINE + mSINE + mSATELLITE + mSIMPLE + mMOBILE_ELEMENT + m_rRNA + mOTHER

perc_aDNA = round(float(aDNA) / aTOTAL_REPEAT_LEN * 100,3)
perc_aLTR = round(float(aLTR) / aTOTAL_REPEAT_LEN * 100,3)
perc_aLINE = round(float(aLINE) / aTOTAL_REPEAT_LEN * 100,3)
perc_aSINE = round(float(aSINE) / aTOTAL_REPEAT_LEN * 100,3)
perc_aSATELLITE = round(float(aSATELLITE) / aTOTAL_REPEAT_LEN * 100,3)
perc_aSIMPLE = round(float(aSIMPLE) / aTOTAL_REPEAT_LEN * 100,3)
perc_aMOBILE_ELEMENT = round(float(aMOBILE_ELEMENT) / aTOTAL_REPEAT_LEN * 100,3)
perc_a_rRNA = round(float(a_rRNA) / aTOTAL_REPEAT_LEN * 100,3)
perc_aOTHER = round(float(aOTHER) / aTOTAL_REPEAT_LEN * 100,3)
arabidopsisTotalRepeatPercentage = perc_aDNA + perc_aLTR + perc_aLINE + perc_aSINE + perc_aSATELLITE + perc_aSIMPLE + perc_aMOBILE_ELEMENT + perc_a_rRNA + perc_aOTHER

perc_mDNA = round(float(mDNA) / mTOTAL_REPEAT_LEN * 100,3)
perc_mLTR = round(float(mLTR) / mTOTAL_REPEAT_LEN * 100,3)
perc_mLINE = round(float(mLINE) / mTOTAL_REPEAT_LEN * 100,3)
perc_mSINE = round(float(mSINE) / mTOTAL_REPEAT_LEN * 100,3)
perc_mSATELLITE = round(float(mSATELLITE) / mTOTAL_REPEAT_LEN * 100,3)
perc_mSIMPLE = round(float(mSIMPLE) / mTOTAL_REPEAT_LEN * 100,3)
perc_mMOBILE_ELEMENT = round(float(mMOBILE_ELEMENT) / mTOTAL_REPEAT_LEN * 100,3)
perc_m_rRNA = round(float(m_rRNA) / mTOTAL_REPEAT_LEN * 100,3)
perc_mOTHER = round(float(mOTHER) / mTOTAL_REPEAT_LEN * 100,3)
maizeTotalRepeatPercentage = perc_mDNA + perc_mLTR + perc_mLINE + perc_mSINE + perc_mSATELLITE + perc_mSIMPLE + perc_mMOBILE_ELEMENT + perc_m_rRNA + perc_mOTHER

A = open('arabidopsis.pieChartValues.txt','w')
M = open('maize.pieChartValues.txt','w')

### arabidopsis
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
A.write("Total Repeat Content\t%s\t%s\n" % (aTOTAL_REPEAT_LEN,arabidopsisTotalRepeatPercentage))

## maize
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
M.write("Total Repeat Content\t%s\t%s\n" % (mTOTAL_REPEAT_LEN,maizeTotalRepeatPercentage))
