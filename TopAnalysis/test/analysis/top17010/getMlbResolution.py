import ROOT
import sys
import numpy as np
from array import array

summary=[]

#loop over events to get the summary of the lb masses
data=ROOT.TChain('analysis/data')
for f in sys.argv[1].split(','):
    data.AddFile(f)

for i in xrange(0,data.GetEntries()):
    data.GetEntry(i)

    #select leptons matched at gen level
    leptons=[]
    for il in xrange(0,data.nl):
        ig=data.l_g[il]
        if ig<0 : continue
        gp4=ROOT.TLorentzVector(0,0,0,0)
        gp4.SetPtEtaPhiM(data.g_pt[ig],data.g_eta[ig],data.g_phi[ig],data.g_m[ig])
        if gp4.Pt()<20 or abs(gp4.Eta())>2.4: continue
        p4=ROOT.TLorentzVector(0,0,0,0)
        p4.SetPtEtaPhiM(data.l_pt[il],data.l_eta[il],data.l_phi[il],data.l_mass[il])
        leptons.append( (data.g_id[ig],gp4,p4) )
    if len(leptons)<2 : continue
    
    #select b-jets matched at gen level
    bjets=[]
    for ij in xrange(0,data.nj):
        ig=data.j_g[ij]
        if ig<0 : continue
        if abs(data.j_hadflav[ij])!=5: continue
        gp4=ROOT.TLorentzVector(0,0,0,0)
        gp4.SetPtEtaPhiM(data.g_pt[ig],data.g_eta[ig],data.g_phi[ig],data.g_m[ig])
        if gp4.Pt()<30 or abs(gp4.Eta())>2.4: continue
        p4=ROOT.TLorentzVector(0,0,0,0)
        p4.SetPtEtaPhiM(data.j_pt[ij],data.j_eta[ij],data.j_phi[ij],data.j_mass[ij])
        bjets.append( (data.j_flav[ij],gp4, p4) )
    if len(bjets)<2 : continue

    #mimic here the selection by ptlb
    lbPairs=[] 
    for il in xrange(0,2):
        for ij in xrange(0,2):           
            lbPairs.append( (leptons[il][1]+bjets[ij][1],
                             leptons[il][2]+bjets[ij][2],
                             leptons[il][0]*bjets[ij][0],
                             il,
                             ij) 
                            )

    lbPairs.sort(key=lambda x: x[0].Pt(), reverse=True)
    leadPair=lbPairs[0]
    summary.append( (leadPair[0].Pt(),leadPair[0].M(),leadPair[1].M()/leadPair[0].M()-1.,leadPair[2]) )
    for ip in xrange(1,4):
        subLeadPair=lbPairs[ip]
        if subLeadPair[3]==leadPair[3]: continue
        if subLeadPair[4]==leadPair[4]: continue
        break
    summary.append( (subLeadPair[0].Pt(),subLeadPair[0].M(),subLeadPair[1].M()/subLeadPair[0].M()-1.,subLeadPair[2]) )

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.04)
c.SetLogy(True)

#compute quantiles and fill distributions
nq=20
summary=np.array(summary)
xq=np.percentile(summary, [i*100./nq for i in xrange(0,nq+1)], axis=0)
xq[0][0]=0.
xq[0][1]=0.
xq[-1][0]=400.
xq[-1][1]=200.

fOut=ROOT.TFile.Open('mlbresol.root','RECREATE')
for i,obs,title in [(0,'ptlb','p_{T}(lb) [GeV]'),
                    (1,'mlb', 'm(lb) [GeV]')]:


    resolH = ROOT.TH1F('resolvs'+obs, '#Deltam(lb)/m(lb);%s;'%title, nq, array('d',xq[:,i]))
    fcorrH = ROOT.TH1F('fcorrvs'+obs, 'correct assignments;%s;'%title, nq, array('d',xq[:,i]))
    obsH   = ROOT.TH1F(obs,           'PDF;%s;'%title, nq, array('d',xq[:,i]))

    totalPairs=float(len(summary))
    for iq in xrange(0,nq):

        minx=xq[iq][i] if iq>0 else 0.
        minxFilt=(summary[:,i]>=minx)
        filtSummary=summary[minxFilt]

        maxx=xq[iq+1][i]
        maxxFilt=(filtSummary[:,i]<maxx)
        filtSummary=filtSummary[maxxFilt]
        totalPass=float(len(filtSummary))

        rq=np.percentile(filtSummary,[16,50,84], axis=0)
        bias=rq[1][2]
        effWid=0.5*((rq[2]-rq[0])[2])
        fcor=float( (filtSummary[:,3]<0).sum() )/totalPass
        
        resolH.SetBinContent(iq+1,effWid)
        fcorrH.SetBinContent(iq+1,fcor)
        obsH.SetBinContent(iq+1,(totalPass/totalPairs)*1./(maxx-minx))

    #plot distribution
    c.cd()
    c.Clear()
    obsH.SetFillStyle(3004)
    obsH.SetFillColor(ROOT.kGray)
    obsH.SetLineColor(1)
    obsH.SetMarkerColor(ROOT.kGray)
    obsH.GetYaxis().SetTitle('Percentage')
    obsH.Draw('hist')
    obsH.SetTitle('PDF (a.u.)')
    obsH.Scale( 40./obsH.GetMaximum() )
    obsH.GetYaxis().SetRangeUser(1,100.)
    resolH.SetLineColor(ROOT.kMagenta+1)
    resolH.SetLineWidth(2)
    resolH.SetLineStyle(9)
    resolH.SetFillStyle(0)
    resolH.Scale(100.)
    resolH.Draw('histsame')    
    fcorrH.SetLineColor(ROOT.kAzure+3)
    fcorrH.SetLineWidth(2)
    fcorrH.SetFillStyle(0)
    fcorrH.Scale(100.)
    fcorrH.Draw('histsame')    
    leg=c.BuildLegend(0.15,0.35,0.4,0.15)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.05)
    leg.SetFillStyle(0)
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,'13 TeV')  
    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']:
        c.SaveAs('mlbresol_%s.%s'%(obs,ext))

    fOut.cd()
    obsH.Write()
    resolH.Write()
    fcorrH.Write()

fOut.Close()
