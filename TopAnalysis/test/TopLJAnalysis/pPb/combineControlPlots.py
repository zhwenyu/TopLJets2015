import ROOT

from prepareWorkspace import showControlPlots

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(False)
    
fIn=ROOT.TFile.Open('plots/controlplots.root')
cats=[      
    ('2b','1l4j2b','#d73027'),
    ('1b','1l4j1b','#fc8d59'),
    ('0b','1l4j1q','#ffffff'),
    ]

for pName in ['m_bwjj','m_wjj','dr_bwjj','dr_jj']:

    allPlots=[]
    total=0
    for c,tag,color in cats:
        fName=pName+'comb'
        print '%s_%s'%(pName,tag)
        allPlots.append( {fName:fIn.Get('%s_%s'%(pName,tag))} )
        allPlots[-1][fName].SetDirectory(0)
        allPlots[-1][fName].SetFillColor(ROOT.TColor.GetColor(color))
        allPlots[-1][fName].SetFillStyle(1001)
        allPlots[-1][fName].SetTitle(c)
        allPlots[-1][fName].SetBinErrorOption(0)
        total+=allPlots[-1][fName].Integral()
    
    print total
    for p in allPlots:
        for fName in p:       
            if p[fName].GetXaxis().GetNbins()%2==0 : p[fName].Rebin()
            p[fName].Scale(1./total)
            p[fName].GetYaxis().SetTitle('1/N dN/dM')
    
    showControlPlots(allPlots)
