import ROOT

from TopLJets2015.TopAnalysis.Plot import *
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

def getROC(hsig,hbkg,fwd=True):
    gr=ROOT.TGraph()

    print '*'*50
    print hsig.GetTitle()
    nbins=hsig.GetNbinsX()
    totSig=hsig.Integral(0,nbins+1)
    totBkg=hbkg.Integral(0,nbins+1)
    for xbin in range(nbins):
        if fwd:
            esig=hsig.Integral(0,xbin+1)/totSig
            ebkg=hbkg.Integral(0,xbin+1)/totBkg
        else:
            esig=hsig.Integral(xbin+1,nbins+1)/totSig
            ebkg=hbkg.Integral(xbin+1,nbins+1)/totBkg
        print xbin,hsig.GetXaxis().GetBinCenter(xbin+1),esig,ebkg
        gr.SetPoint(xbin,esig,ebkg)

    return gr


plots={}
profiles={}
rocs=[]

dyT=ROOT.TChain('tree')
for i in range(5):
    dyT.AddFile('/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/Chunks/MC13TeV_2017_DY50toInf_fxfx_%d.root'%i)

for ch,chid,trigReq in [('mm',13*13,'(hasMMTrigger || hasMTrigger)'),
                        ('ee',11*11,'hasEETrigger')]:

    sig='/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/Chunks/MC13TeV_Z%s_m_X_950_xangle_120_2017_postTS2_fullsim_0.root'%ch
    sigF=ROOT.TFile.Open(sig)
    sigT=sigF.Get('tree')

    finalSel='bosonpt>50 && isZ && abs(l1id*l2id)==%d && %s'%(chid,trigReq)
    for t,v,c in [("p_{T}(ll) [GeV]",                 "bosonpt >> {0}(20,0,250)",   'isZ'),
                  ('m(ll) [GeV]',                     "mboson >> {0}(20,81,101)",    finalSel),
                  ('Jet multiplicity',                "nj >> {0}(4,0,4)",            finalSel),              
                  ('Missing transverse energy [GeV]', "met_pt >> {0}(20,0,200)",     finalSel),
                  ('#Delta#phi[E_{T}^{miss},p_{T}(ll)] [rad]', "fabs(TVector2::Phi_mpi_pi(met_phi-bosonphi)) >> {0}(20,0,3.15)",     finalSel),
                  ('Charged particle multiplicity',   "nchPV >> {0}(20,0,100)",       finalSel),                  
                  ('p_{T}^{*}(vtx)',                      'fabs(sumPVChPt-bosonpt) >> {0}(40,0,50)',  finalSel),
                  ('[p_{T}(vtx)-p_{T}(ll)]/[p_{T}(1)+p_{T}(2)]',                  '(sumPVChPt-bosonpt)/(l1pt+l2pt) >> {0}(40,-1,1)',  finalSel),
                  ('p_{T}(vtx)/p_{T}(ll)-1',                  'sumPVChPt/bosonpt-1 >> {0}(40,-1,1)',  finalSel),
                  ('p_{T}(vtx)-p_{T}(ll)',            '(sumPVChPt-bosonpt) >> {0}(20,-200,50)',  finalSel),
                  ('N(muons)',                        'nrawmu-2 >> {0}(10,0,10)',      finalSel),
                  ]:

        pname=v.split(' >>')[0]
        for tkn in ['(',')','-','+','/']: pname=pname.replace(tkn,'')
        pname=pname.replace('TVector2::Phi_mpi_pi','dphi_')
        if not pname in plots: 
            plots[pname]=[]            
            profiles[pname]=[]

        #1D
        sigT.Draw(v.format('sig'),'evwgt*(%s)'%c)
        sigH=ROOT.gDirectory.Get('sig')
        sigH.Scale(1./sigH.Integral())
        sigH.GetXaxis().SetTitle(t)
        sigH.GetYaxis().SetTitle('PDF')
        
        dyT.Draw(v.format('dy'),'evwgt*(%s)'%c)
        dyH=ROOT.gDirectory.Get('dy')
        dyH.Scale(1./dyH.Integral())
        dyH.GetXaxis().SetTitle(t)
        dyH.GetYaxis().SetTitle('PDF')

        plots[pname].append( sigH.Clone('sig%s_%s'%(pname,ch)) )
        plots[pname][-1].SetDirectory(0)
        plots[pname].append( dyH.Clone('dy%s_%s'%(pname,ch)) )
        plots[pname][-1].SetDirectory(0)
        sigH.Delete()
        dyH.Delete()

        #2D
        v2d=v.replace('}(','}(5,0,60,')
        v2d=v2d.replace('>>',':nvtx>>')
        
        sigT.Draw(v2d.format('sig'),'evwgt*(%s)'%c)
        sigH=ROOT.gDirectory.Get('sig')
        sigH.Scale(1./sigH.Integral())
        sigH.GetXaxis().SetTitle('Vertex multiplicity')
        sigH.GetYaxis().SetTitle(t)
        sigH.GetZaxis().SetTitle('PDF')
        
        dyT.Draw(v2d.format('dy'),'evwgt*(%s)'%c)
        dyH=ROOT.gDirectory.Get('dy')
        dyH.Scale(1./dyH.Integral())
        dyH.GetXaxis().SetTitle('Vertex multiplicity')
        dyH.GetYaxis().SetTitle(t)
        dyH.GetZaxis().SetTitle('PDF')

        if len(profiles[pname])==0:
            profiles[pname].append(sigH.Clone('sig2d_%s'%ch))
            profiles[pname][-1].SetDirectory(0)
            profiles[pname].append(dyH.Clone('dy2d_%s'%ch))
            profiles[pname][-1].SetDirectory(0)
        else:
            profiles[pname][0].Add(sigH)
            profiles[pname][1].Add(dyH)

        sigH.Delete()
        dyH.Delete()


for pname in plots:
    p=Plot(pname,com='13 TeV')
    p.frameMin=0
    p.add(plots[pname][0],  title='ppZ#mu#muX(950)',        color=ROOT.kBlack,     isData=False,  spImpose=True, isSyst=False)
    p.add(plots[pname][2],  title='ppZeeX(950)',            color=ROOT.kGray,     isData=False,  spImpose=True, isSyst=False)
    p.add(plots[pname][1],   title='Z/#gamma^{*}#rightarrow#mu#mu', color=ROOT.kTeal+1, isData=False,  spImpose=True, isSyst=False)
    p.add(plots[pname][3],   title='Z/#gamma^{*}#rightarrowee',     color=ROOT.kTeal+3, isData=False,  spImpose=True, isSyst=False)
    p.show(outDir='./', lumi=37500,)

    hsig=plots[pname][0].Clone('hsig')
    hsig.SetTitle(pname)
    #hsig.Add(plots[pname][2])
    hbkg=plots[pname][1].Clone('hbkg')
    #hbkg.Add(plots[pname][3])
    fwd=False
    if pname in ['nj','nchPV','nrawmu2','fabssumPVChPtbosonpt'] : fwd=True
    marker=20+len(rocs)
    rocs.append( getROC(hsig,hbkg,fwd) )
    rocs[-1].SetName(pname)
    roctitle=hsig.GetXaxis().GetTitle()
    roctitle=roctitle.replace('[GeV]','')
    roctitle=roctitle.replace('[rad]','')
    rocs[-1].SetTitle(roctitle)
    rocs[-1].SetMarkerStyle(marker)
    rocs[-1].SetFillStyle(0)
    hsig.Delete()
    hbkg.Delete()


    c=ROOT.TCanvas('c','c',600,600)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)

    p1=profiles[pname][0].ProfileX()
    p1.SetMarkerStyle(20)
    p1.SetTitle('ppZX(950)')
    p1.GetYaxis().SetTitle('<%s>'%plots[pname][0].GetXaxis().GetTitle())
    p1.SetLineColor(1)
    p1.SetMarkerColor(1)
    p1.GetXaxis().SetLabelSize(0.035)
    p1.GetYaxis().SetLabelSize(0.035)
    p1.GetXaxis().SetTitleSize(0.04)
    p1.GetYaxis().SetTitleSize(0.04)
    p1.GetYaxis().SetTitleOffset(0.8)
    p1.Draw()

    p2=profiles[pname][1].ProfileX()
    p2.SetMarkerStyle(24)
    p2.SetTitle('Z/#gamma^{*}')
    p2.SetLineColor(ROOT.kTeal+3)
    p2.SetMarkerColor(ROOT.kTeal+3)
    p2.Draw('same')

    
    minY=min([plots[pname][i].GetMean()-plots[pname][i].GetRMS() for i in range(4)])
    #if minY<0 : minY *=1.5
    #else      : minY *=0.5
    maxY=max([plots[pname][i].GetMean()+plots[pname][i].GetRMS() for i in range(4)])
    #if maxY<0 : maxY *=1.5
    #else      : maxY *=0.5
    p1.GetYaxis().SetRangeUser(minY,maxY)

    leg=c.BuildLegend(0.6,0.15,0.9,0.25)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.SetTextFont(42)

    for ext in ['png','pdf']:
        c.SaveAs(pname+'_vsvtx.%s'%ext)
        
    
c=ROOT.TCanvas('c','c',600,600)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetBottomMargin(0.1)
mg=ROOT.TMultiGraph()
for r in rocs:
    mg.Add(r)
mg.Draw('ap')
mg.GetXaxis().SetRangeUser(9e-1,1)
mg.GetYaxis().SetRangeUser(1e-1,1)
mg.GetXaxis().SetTitle('Signal efficiency')
mg.GetYaxis().SetTitle('Background efficiency')
mg.GetXaxis().SetMoreLogLabels()
mg.GetYaxis().SetMoreLogLabels()

leg=c.BuildLegend(0.15,0.5,0.5,0.12)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.035)
leg.SetTextFont(42)

l=ROOT.TLine()
l.SetLineStyle(2)
l.SetLineWidth(2)
l.SetLineColor(ROOT.kGray)
l.DrawLine(0.99,1e-1,0.99,1)

#c.SetLogx()
#c.SetLogy()

for ext in ['png','pdf']:
    c.SaveAs('exc_rocs.%s'%ext)
