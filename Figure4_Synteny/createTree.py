import sys, os, re
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, faces, AttrFace, PhyloTree, CircleFace, RectFace, faces

def readHitFile(hitFile):
    topHit = {}
    topScore = {}
    with open(hitFile,'r') as HF:
        for line in HF:
            geneID,uniprotID,uniprotDesc,pIdent,eValue,bitScore,queryCov = line.strip().split('\t')
            eValue = float(eValue)
            bitScore = float(bitScore)
            if geneID in topHit:
                if bitScore > topScore[geneID]:
                    topScore[geneID] = bitScore
                    topHit[geneID] = (uniprotID,uniprotDesc,eValue,bitScore)
            else:
                topHit[geneID] = (uniprotID,uniprotDesc,eValue,bitScore)
                topScore[geneID] = bitScore
    return(topHit)

def search(tree,topHit1,topHit2):
    colors = {"hopID":"#a50026", "canID":"#313695"}
    
    hopCount = 0
    canCount = 0

    for leaf in tree.traverse():
        leaf.img_style['size'] = 0
        if "F" in leaf.name:
            if leaf.is_leaf():
                geneID = 'hopID'
                color = colors.get(geneID, None)
                if color:
                    style1 = NodeStyle()
                    style1["fgcolor"] = "#a50026"
                    style1["size"] = 0
                    style1["vt_line_color"] = "#a50026"
                    style1["hz_line_color"] = "#a50026"
                    style1["vt_line_width"] = 2
                    style1["hz_line_width"] = 2
                    style1["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
                    style1["hz_line_type"] = 0
                    leaf.set_style(style1)
                    if leaf.name in topHit1:
                        hopCount += 1
                        uniprotID1,uniprotDesc1,eValue1,bitScore1 = topHit1[leaf.name]
                        if 'Berberine' in uniprotDesc1:
                            newLeafName = 'H. lupulus, BBE-like, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style1["bgcolor"] = "#e0f3f8"
                            leaf.set_style(style1)
                        elif 'Inactive tetrahydrocannabinolic acid synthase' in uniprotDesc1:
                            newLeafName = 'H. lupulus, Inactive THCAS, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style1["bgcolor"] = "#4393c3" 
                            leaf.set_style(style1)
                        elif 'Tetrahydrocannabinolic acid synthase ' in uniprotDesc1:
                            newLeafName = 'H. lupulus, THCAS, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style1["bgcolor"] = "#4393c3"
                            leaf.set_style(style1)
                        elif 'Cannabidiolic acid synthase-like' in uniprotDesc1:
                            newLeafName = 'H. lupulus, CBDAS-like, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style1["bgcolor"] = "#92c5de" 
                            leaf.set_style(style1)
                        elif 'Cannabidiolic acid synthase ' in uniprotDesc1:
                            newLeafName = 'H. lupulus, CBDAS, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style1["bgcolor"] = "#92c5de"
                            leaf.set_style(style1)
                        else:
                            print(leaf.name,uniprotID1,uniprotDesc1)
        else:
            if leaf.is_leaf():
                geneID = 'canID'
                color = colors.get(geneID, None)
                if color:
                    style2 = NodeStyle()
                    style2["fgcolor"] = "#313695"
                    style2["size"] = 0
                    style2["vt_line_color"] = "#313695"
                    style2["hz_line_color"] = "#313695"
                    style2["vt_line_width"] = 2
                    style2["hz_line_width"] = 2
                    style2["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
                    style2["hz_line_type"] = 0
                    leaf.set_style(style2)
                    if leaf.name in topHit2:
                        canCount += 1
                        uniprotID2,uniprotDesc2,eValue2,bitScore2 = topHit2[leaf.name]
                        if 'Berberine' in uniprotDesc2:
                            newLeafName = 'C. sativa, BBE-like, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style2["bgcolor"] = "#e0f3f8"
                            leaf.set_style(style2)
                        elif 'Inactive tetrahydrocannabinolic acid synthase' in uniprotDesc2:
                            newLeafName = 'C. sativa, Inactive THCAS, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style2["bgcolor"] = "#4393c3"
                            leaf.set_style(style2)
                        elif 'Tetrahydrocannabinolic acid synthase ' in uniprotDesc2:
                            newLeafName = 'C. sativa, THCAS, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style2["bgcolor"] = "#4393c3"
                            leaf.set_style(style2)
                        elif 'Cannabidiolic acid synthase-like' in uniprotDesc2:
                            newLeafName = 'C. sativa, CBDAS-like, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style2["bgcolor"] = "#92c5de"
                            leaf.set_style(style2)
                        elif 'Cannabidiolic acid synthase ' in uniprotDesc2:
                            newLeafName = 'C. sativa, CBDAS, ' + leaf.name
                            name_face = TextFace(newLeafName, fgcolor=color, fsize=12)
                            leaf.add_face(name_face, column=0, position='branch-right')
                            style2["bgcolor"] = "#92c5de"
                            leaf.set_style(style2)
                        else:
                            print(leaf.name,uniprotID1,uniprotDesc1)

    ts = TreeStyle()
    #ts.rotation = 90
    ts.mode = "c"
    ts.arc_start = -180 
    ts.arc_span = 180
    # ts.optimal_scale_level = 'mid'
    
    ts.branch_vertical_margin = 10
    # ts.scale = 180
    ts.scale = 225
    ts.show_leaf_name = False
    ts.show_scale = False
    tree.render("hop_vs_can_tree_v5.svg",tree_style=ts,w=600)
    tree.render("hop_vs_can_tree_v5.pdf",tree_style=ts,w=600)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <tree file> <hit file 1> <hit file 2> \n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

treeFile = sys.argv[1]
hitFile1 = sys.argv[2]
hitFile2 = sys.argv[3]

t = Tree(treeFile)

topHitDict1 = readHitFile(hitFile1)
topHitDict2 = readHitFile(hitFile2)
search(t,topHitDict1,topHitDict2)
