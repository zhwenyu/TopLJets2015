#!/usr/bin/env python

import sys
import ROOT

"""
associate different metrics to the dijet systems in an event
"""
def dijetsWithMetric(jets):

    #sort by decreasing pT
    sorted(jets, key=lambda p4 : p4.Pt(), reverse=True)
    
    #build all possible dijet pairs with different metrics for ranking
    dijets=[]
    for i in xrange(0,len(jets)):
        for j in xrange(i+1,len(jets)):
            alpha = [ 60./(jets[i].Pt()+jets[j].Pt()), 60./(jets[i]+jets[j]).Pt() ]
            for p in [-2, -1, 0, 1 ,2 ]:
                for metric in ['dr','deta']:
                    kti2p=ROOT.TMath.Power(jets[i].Pt(),2*p)
                    ktj2p=ROOT.TMath.Power(jets[j].Pt(),2*p)
                    deltaij=jets[i].DeltaR(jets[j]) if metric=='dr' else ROOT.TMath.Abs(jets[i].Eta()-jets[j].Eta())
                    dij=ROOT.TMath.Min(kti2p,ktj2p)*deltaij
                    diB=kti2p
                    alpha.append(dij/diB)
            dijets.append( (jets[i]+jets[j], alpha) )

    return dijets

"""
summarize resolutions
"""
def summarizeRankingPerformance(h,window,relativeresol):

    nbinsx=h.GetNbinsX()
    resolH=ROOT.TH1F('resol',';Ranking mode;Resolution',nbinsx,0,nbinsx)
    if relativeresol : resolH.GetYaxis().SetTitle('Relative resolution')
    resolH.SetMarkerStyle(20)
    offPeakFracH=ROOT.TH1F('offpeakfrac',';Ranking mode;N(off-peak)/N(on-peak)',nbinsx,0,nbinsx)
    offPeakFracH.SetMarkerStyle(20)
    

    for xbin in xrange(0,nbinsx):
        py=h.ProjectionY('py',xbin+1,xbin+1)
        resolH.GetXaxis().SetBinLabel(xbin+1,h.GetXaxis().GetBinLabel(xbin+1))
        offPeakFracH.GetXaxis().SetBinLabel(xbin+1,h.GetXaxis().GetBinLabel(xbin+1))

        #determine maximum
        ycenBin=py.GetMaximumBin()
        ycen=py.GetXaxis().GetBinCenter(ycenBin)

        #count off-peak to on-peak fraction
        offPeak,onPeak=0.,0.
        offPeakUnc,onPeakUnc=0.,0.
        for ybin in xrange(0,py.GetNbinsX()):
            yval=py.GetXaxis().GetBinCenter(ybin+1)
            cts,ctsUnc=py.GetBinContent(ybin+1),py.GetBinError(ybin+1)
            if ROOT.TMath.Abs(yval-ycen)>window : 
                offPeak += cts
                offPeakUnc += ctsUnc**2
            else : 
                onPeak += cts
                onPeakUnc += ctsUnc
        frac=offPeak/onPeak
        fracUnc=ROOT.TMath.Sqrt(onPeakUnc*(offPeak**2)+offPeakUnc*(onPeak**2))/onPeak
        offPeakFracH.SetBinContent(xbin+1,frac)
        offPeakFracH.SetBinError(xbin+1,fracUnc)
            
        #resolution
        resol=py.GetRMS()
        if relativeresol : resol=resol/ycen
        resolUnc=resol*py.GetRMSError()/py.GetRMS()
        resolH.SetBinContent(xbin+1,resol)
        resolH.SetBinError(xbin+1,resolUnc)

        py.Delete()

    #
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',800,500)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetTopMargin(0)
    c.SetBottomMargin(0)
    p1 = ROOT.TPad('p1','p1',0.0,0.5,1.0,0.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.008)
    p1.SetBottomMargin(0.12)
    p1.Draw()
    p1.cd()
    resolH.Draw('e1')
    resolH.GetXaxis().SetLabelSize(0.06)
    resolH.GetXaxis().SetTitleSize(0.06)
    resolH.GetYaxis().SetLabelSize(0.06)
    resolH.GetYaxis().SetTitleSize(0.06)
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.5,1.0,1.0)
    p2.SetBottomMargin(0.005)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.05)
    p2.Draw()
    p2.cd()
    offPeakFracH.Draw('e1')
    offPeakFracH.GetXaxis().SetLabelSize(0.08)
    offPeakFracH.GetXaxis().SetTitleSize(0.08)
    offPeakFracH.GetYaxis().SetLabelSize(0.06)
    offPeakFracH.GetYaxis().SetTitleSize(0.06)
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.08)
    label.DrawLatex(0.6,0.85,'#bf{CMS} #it{simulation preliminary}')
    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s_resol.%s'%(h.GetName(),ext))


"""
"""
def main():

    histos={}
    histos['ranked_mjj']      = ROOT.TH2F('ranked_mjj',      ';Ranking mode;Dijet invariant mass [GeV]', 12,0,12,50,0,200)
    histos['ranked_mw']       = ROOT.TH2F('ranked_mw',       ';Ranking mode;Closest W mass [GeV]',       12,0,12,50,0,200)
    histos['ranked_mjjresol'] = ROOT.TH2F('ranked_mjjresol', ';Ranking mode;M(jj)-M(W) [GeV]',           12,0,12,50,-50,50)
    histos['ranked_drresol']  = ROOT.TH2F('ranked_drresol',  ';Ranking mode;#DeltaR(jj,W)',              12,0,12,50,0,1)
    for xbin,label in [(1,'#Sigmap_{T}'),   (2,'p_{T}(jj)'),
                       (3,'#DeltaR,p=-2'), (4,'#Delta#eta,p=-2'),
                       (5,'#DeltaR,p=-1'), (6,'#Delta#eta,p=1'),
                       (7,'#DeltaR,p=0'),  (8,'#Delta#eta,p=0'),
                       (9,'#DeltaR,p=+1'), (10,'#Delta#eta,p=+1'),
                       (11,'#DeltaR,p=+2'),(12,'#Delta#eta,p=+2')]:
        histos['ranked_mjj'].GetXaxis().SetBinLabel(xbin,label)
        histos['ranked_mw'].GetXaxis().SetBinLabel(xbin,label)
        histos['ranked_mjjresol'].GetXaxis().SetBinLabel(xbin,label)
        histos['ranked_drresol'].GetXaxis().SetBinLabel(xbin,label)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetDirectory(0)
        
    #loop over the events
    data=ROOT.TChain('data')
    for url in sys.argv[1].split(','):
        data.AddFile(url)
    nentries=data.GetEntries()
    for i in xrange(0,nentries):
        data.GetEntry(i)

        #get reco and gen jets which are not b-tagged
        jets,genJets=[],[]
        matchedGenJetsIdx=[]
        for j in xrange(0,data.nj):
            if data.j_btag[j]: continue        
            jets.append(ROOT.TLorentzVector())
            jets[-1].SetPtEtaPhiM( data.j_pt[j], data.j_eta[j], data.j_phi[j], data.j_m[j])

            try:
                for k in xrange(0,data.ngj):
                    if k in matchedGenJetsIdx : continue
                    gjp4=ROOT.TLorentzVector()
                    gjp4.SetPtEtaPhiM( data.gj_pt[k], data.gj_eta[k], data.gj_phi[k], data.gj_m[k])
                    if gjp4.DeltaR(jets[-1])>0.4: continue
                    matchedGenJetsIdx.append(k)
                    genJets.append( ROOT.TLorentzVector(gjp4) )
                    break
            except:
                pass

        #get Ws
        Wbosons=[]
        try:
            for k in xrange(0,data.ngp):
                Wbosons.append(ROOT.TLorentzVector())
                Wbosons[-1].SetPtEtaPhiM( data.gp_pt[k], data.gp_eta[k], data.gp_phi[k], data.gp_m[k])
        except:
            pass

        #two jets are reco level
        if len(jets)<2 : continue
                    
        #fill mass and direction resolutions
        weight=data.w
        if len(genJets)>=2 and len(Wbosons)==2:
            genDijets=dijetsWithMetric(genJets)
            for mode in xrange(0,12):
                jj,rank=sorted(genDijets, key=lambda jjranks : jjranks[1][mode])[0]

                dRjjw1,dRjjw2=jj.DeltaR(Wbosons[0]),jj.DeltaR(Wbosons[1])
                widx=0 if dRjjw1<dRjjw2 else 1

                histos['ranked_mjj'].Fill(mode,jj.M(),weight)
                histos['ranked_mw'].Fill(mode,Wbosons[widx].M(),weight)
                histos['ranked_mjjresol'].Fill(mode,jj.M()-Wbosons[widx].M(),weight)
                histos['ranked_drresol'].Fill(mode,jj.DeltaR(Wbosons[widx]),weight)

    #build summary plots
    summarizeRankingPerformance(histos['ranked_mjjresol'],15,True)
    summarizeRankingPerformance(histos['ranked_drresol'],0.1,False)

    #save to file
    fOut=ROOT.TFile('workspace.root','RECREATE')
    for h in histos:
        histos[h].SetDirectory(fOut)
        histos[h].Write()
    fOut.Close()
            

if __name__ == "__main__":
    main()
