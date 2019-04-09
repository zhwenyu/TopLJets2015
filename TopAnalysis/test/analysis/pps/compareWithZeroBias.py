import ROOT
import os

BASEDIR='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/Chunks/'

def getTree(fList,baseDir=BASEDIR):
    tree=ROOT.TChain('tree')
    pudisc=ROOT.TChain('pudiscr')
    for f in fList: 
        tree.AddFile(baseDir+f)
        pudisc.AddFile(baseDir+'/pudiscr/'+f)
    tree.AddFriend(pudisc)
    return tree


def getHistos(t,name,title,lw,lc,ms,weights=None):
    
    histos={
        'xangle':ROOT.TH1F('xangle'+name,';LHC beam crossing angle;PDF',4,120,160),
        'nvtx_120':ROOT.TH1F('nvtx_120'+name,';Vertex multiplicity;PDF',100,0,100),
        'rho_120':ROOT.TH1F('rho_120'+name,';#rho;PDF',25,0,50),
        'met_120':ROOT.TH1F('met_120'+name,';Missing transverse energy [GeV];PDF',25,0,250),
        'nj_120':ROOT.TH1F('nj_120'+name,';Jet multiplicity;PDF',10,0,10),
        'ht_120':ROOT.TH1F('ht_120'+name,';H_{T} [GeV];PDF',25,0,250),
        'nrptk_120':ROOT.TH1F('nrptk_120'+name,';Proton multiplicity;PDF',10,0,10),
        'hfmult_120':ROOT.TH1F('hfmult_120'+name,';PF multiplicity (HF);PDF',50,0,1000),
        'hfht_120':ROOT.TH1F('hfht_120'+name,';PF sum p_{T}(HF) [GeV];PDF',50,0,1000),
        'hfpz_120':ROOT.TH1F('hfpz_120'+name,';PF sum |p_{z}|(HF) [TeV];PDF',50,0,40),
        'rfc_120':ROOT.TH1F('rfc_120'+name,';Random forest classifier;PDF',25,0,1),
        }
    for key in histos.keys():
        for xangle in ['130','140','150']:
            newKey=key.replace('120',xangle)
            histos[newKey]=histos[key].Clone(newKey)


    selZB=False
    if name=='zb' or name=='wzb': selZB=True

    nsel=0
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        
        nvtx=t.nvtx
        xangle=round(t.beamXangle)

        if not xangle in [120,130,140,150]: continue
        if selZB :
            if not t.hasZBTrigger : continue
        else:
            if t.evcat!=13*13 : continue
            if not t.isZ : continue
            if t.bosonpt>10: continue


        nsel+=1
        if nsel>10000: break

        w=t.evwgt
        if weights and xangle in weights:
            w *= weights[xangle].GetBinContent( weights[xangle].GetXaxis().FindBin(t.rho))

        histos['xangle'].Fill(xangle,w)        
        histos['nvtx_%d'%xangle].Fill(nvtx,w)
        histos['rho_%d'%xangle].Fill(t.rho,w)
        histos['met_%d'%xangle].Fill(t.met_pt,w)
        histos['nj_%d'%xangle].Fill(t.nj+t.nb,w)
        histos['ht_%d'%xangle].Fill(t.htj+t.htb,w)
        histos['nrptk_%d'%xangle].Fill(t.nRPtk,w)
        histos['rfc_%d'%xangle].Fill(getattr(t,'rfc_%d'%xangle),w)
        histos['hfmult_%d'%xangle].Fill(t.PFMultSumHF,w)
        histos['hfht_%d'%xangle].Fill(t.PFHtSumHF,w)
        histos['hfpz_%d'%xangle].Fill(t.PFPzSumHF/1000.,w)

    for key in histos:
        histos[key].SetTitle(title)
        histos[key].Sumw2()
        histos[key].SetDirectory(0)
        histos[key].SetLineColor(lc)
        histos[key].SetMarkerColor(lc)
        histos[key].SetMarkerStyle(ms)
        histos[key].SetLineWidth(lw)
        tot=histos[key].Integral()
        if tot!=0:
            histos[key].Scale(1./tot)

    return histos


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
from TopLJets2015.TopAnalysis.Plot import *

for era in 'C':

    zbfiles=[f for f in os.listdir(BASEDIR+'/pudiscr/') if '2017%s_ZeroBias'%era in f][0:1]
    zb=getTree(zbfiles)

    mmfiles=[f for f in os.listdir(BASEDIR+'/pudiscr/') if '2017%s_DoubleMuon'%era in f][0:1]
    mm=getTree(mmfiles)
    
    print '\t ZeroBias histos'
    zbh  = getHistos(zb, 'zb', 'ZeroBias',                    lw=2,lc=ROOT.kGray,ms=1)

    print '\t Zmm histos'
    mmh  = getHistos(mm, 'mm', 'Z#rightarrow#mu#mu',          lw=1,lc=1,ms=20)

    weights={}
    for xangle in [120,130,140,150]:
        weights[xangle]=mmh['rho_%d'%xangle].Clone('ratio_%d'%xangle)  
        weights[xangle].Divide(zbh['rho_%d'%xangle])
        weights[xangle].SetDirectory(0)

    print '\t weighted ZeroBias histos'
    wzbh= getHistos(zb, 'wzb', 'ZeroBias (weighted)',   lw=1,lc=ROOT.kRed,ms=1,weights=weights)


    for key in zbh:  

        p=Plot('%s_2017%s'%(key,era),com='13 TeV')
        p.doPoissonErrorBars=False
        p.add(zbh[key],  title=zbh[key].GetTitle(),  color=ROOT.kGray, isData=False, spImpose=False, isSyst=False)
        p.add(wzbh[key], title=wzbh[key].GetTitle(), color=ROOT.kRed,  isData=False, spImpose=False, isSyst=False)
        p.add(mmh[key],  title=mmh[key].GetTitle(),  color=1,          isData=True,  spImpose=False, isSyst=False)
        p.show(outDir='./',lumi=41833,noStack=True)

