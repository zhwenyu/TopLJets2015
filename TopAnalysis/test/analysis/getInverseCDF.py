import ROOT
import numpy as np

methodList=[('TMVA_HighMJJ.root', 'vbf/TrainTree', 'BDT_VBF0HighMJJ'),
            #('TMVA_LowMJJ.root',  'vbf/TrainTree', 'BDT_VBF0HighMJJ')
            ]

cdfInv=[]
for f,t,m in methodList:

    #build list of mva values for background
    mva=[]
    fIn=ROOT.TFile.Open(f)
    tree=fIn.Get(t)
    for i in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(i)
        if tree.classID==0: continue
        mva.append( getattr(tree,m) )
    fIn.Close()

    #compute quantiles
    p=np.percentile(mva, range(1,99) )
    
    cdfInv.append(ROOT.TGraph())
    cdfInv[-1].SetName(m+'_cdfinv')
    cdfInv[-1].SetMarkerStyle(20)
    cdfInv[-1].SetPoint(0,-1,0)
    for i in xrange(0,len(p)):
        cdfInv[-1].SetPoint(i+1,p[i],(i+1.)/100.)
    cdfInv[-1].SetPoint(len(p)+1,1,1)

#save results
fOut=ROOT.TFile('inverse_cdfs.root','RECREATE')
for gr in cdfInv: gr.Write()
fOut.Close()
