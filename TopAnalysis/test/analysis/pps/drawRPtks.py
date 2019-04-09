import ROOT
import os
from runExclusiveAnalysis import getTracksPerRomanPot

def cmsHeader():
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.97,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(32)
    txt.DrawLatex(0.96,0.97,'#scale[0.9]{37.5 fb^{-1} (13 TeV)}')


def doControlPlots(period):
    baseDir='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/Chunks/'
    tree=ROOT.TChain('tree')
    for f in os.listdir(baseDir):
        if not 'Data13TeV_2017%s_DoubleMuon'%period in f : continue
        tree.AddFile(os.path.join(baseDir,f))
    
    histos={'csi023':ROOT.TH1F('csi023','RP023;#xi;PDF',50,0,0.3),
            'csi123':ROOT.TH1F('csi123','RP123;#xi;PDF',50,0,0.3),
            'csi103':ROOT.TH1F('csi103','RP103;#xi;PDF',50,0,0.3),
            'csi003':ROOT.TH1F('csi003','RP003;#xi;PDF',50,0,0.3),
            'csi2d':ROOT.TH2F('csi2d',';#xi;#xi;PDF',50,0,0.3,50,0,0.3),
            'csineg_2d':ROOT.TH2F('csineg_2d',';#xi;#xi;PDF',50,0,0.3,50,0,0.3),
            'csipos_2d':ROOT.TH2F('csipos_2d',';#xi;#xi;PDF',50,0,0.3,50,0,0.3),
            'n023':ROOT.TH1F('n023','RP023;Proton multiplicity;PDF',5,0,5),
            'n123':ROOT.TH1F('n123','RP123;Proton multiplicity;PDF',5,0,5),        
            'n003':ROOT.TH1F('n003','RP003;Proton multiplicity;PDF',5,0,5),
            'n103':ROOT.TH1F('n103','RP103;Proton multiplicity;PDF',5,0,5),        
            'n2d':ROOT.TH2F('n2d',';Proton multiplicity; Proton multiplicity; PDF',5,0,5,5,0,5),
            'nneg_2d':ROOT.TH2F('nneg_2d',';Proton multiplicity; Proton multiplicity; PDF',5,0,5,5,0,5),
            'npos_2d':ROOT.TH2F('npos_2d',';Proton multiplicity; Proton multiplicity; PDF',5,0,5,5,0,5)
            }
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        if not tree.isZ : continue
        if tree.bosonpt>10 : continue
        rp023,rp123=getTracksPerRomanPot(tree)
        rp003,rp103=getTracksPerRomanPot(tree,False,False)
        histos['n2d'].Fill(len(rp023),len(rp123))
        histos['npos_2d'].Fill(len(rp003),len(rp023))
        histos['nneg_2d'].Fill(len(rp103),len(rp123))
        histos['n023'].Fill(len(rp023))
        histos['n123'].Fill(len(rp123))
        histos['n003'].Fill(len(rp003))
        histos['n103'].Fill(len(rp103))
        for x in rp003: histos['csi003'].Fill(x)
        for x in rp023: histos['csi023'].Fill(x)
        for x in rp103: histos['csi103'].Fill(x)
        for x in rp123: histos['csi123'].Fill(x)
        if len(rp023)==1 and len(rp123)==1:
            histos['csi2d'].Fill(rp023[0],rp123[0])
            csi_rp003=rp003[0] if len(rp003)>0  else 0.
            csi_rp103=rp103[0] if len(rp103)>0  else 0.
            histos['csipos_2d'].Fill(csi_rp003,rp023[0])
            histos['csineg_2d'].Fill(csi_rp103,rp123[0])
            
    for key in histos:
        histos[key].SetDirectory(0)
        histos[key].Sumw2()
        histos[key].Scale(1./histos[key].Integral())

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.15)
    for d in ['csi','n']:
        c.SetRightMargin(0.03)

        drawOpt='e1'
        ctr=0
        for rp in ['003','023','103','123']:
            histos['%s%s'%(d,rp)].Draw(drawOpt)
            histos['%s%s'%(d,rp)].SetLineColor(1)
            histos['%s%s'%(d,rp)].SetMarkerColor(1+ctr%2)
            histos['%s%s'%(d,rp)].SetMarkerStyle(20+ctr/2)
            drawOpt='e1same'
            ctr+=1
        leg=c.BuildLegend(0.75,0.94,0.95,0.8)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        cmsHeader()
        c.Modified()
        c.Update()
        for ext in ['pdf','png']:
            c.SaveAs('%s_%s.%s'%(d,period,ext))
    
        c.SetRightMargin(0.15)
        for pfix in ['','neg_','pos_']:
            histos['%s%s2d'%(d,pfix)].Draw('colztext' if d =='n' else 'colz')
            cmsHeader()
            c.Modified()
            c.Update()
            for ext in ['pdf','png']:
                c.SaveAs('%s%s2d_%s.%s'%(d,pfix,period,ext))

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPaintTextFormat("4.2f");
for period in 'BCDEF':
    doControlPlots(period)


    


