#!/usr/bin/env python2.7

import ROOT
import optparse
import os,sys
import random

from kinTools import *

weightVar = ROOT.RooRealVar("w","w",1)

EVENTCATEGORIES=[
    '1l4j2b','1l4j1b','1l4j1q','1f4j1q','1f4j1b','1f4j2b',
    '1l3j1b','1l3j1q','1f3j1q','1f3j1b'
    ]
    
DATA={'pPb iso'    :('data/Data13TeV_pPb.root','#4575b4'),
      'Pbp iso'    :('data/Data13TeV_Pbp.root','#fc8d59'),
      'pPb non-iso':('data/Data13TeV_pPb_noniso.root','#91bfdb'),
      'Pbp non-iso':('data/Data13TeV_Pbp_noniso.root','#fee090')}

    
"""
starts the workspace
"""
def startWorkspace(opt):
    
    #start the workspace
    ws  = ROOT.RooWorkspace('w')

    #list of variables of interest
    args=ROOT.RooArgSet()
    varList=['y_l[-4,4]',        'y_bhad[-4,4]',   'y_blep[-4,4]', #particle level
             'm_wjj[0,500]',     'dr_jj[0,6]',                     #w hadronic level
             'pt_wl[0,500]',     'y_wl[-4,4]',                     #w lep level
             'm_bwjj[0,500]',    'y_bwjj[-4,4]',                   #hadronic top level
             'm_bwl[0,500]',     'y_bwl[-4,4]',                    #leptonic top level
             'dy_bwjjbwl[0,5]',  'sy_bwjjbwl[0,5]'                 #top vs anti-top
             ]
    for var in varList :
        args.add( ws.factory(var) )
    if opt.mc:
        args.add( weightVar )
    #categories for the data
    sample=ROOT.RooCategory('sample','sample')
    for name in EVENTCATEGORIES :
        sample.defineType(name)
    getattr(ws,'import')(sample)
    args.add(ws.cat('sample'))

    #define the dataset
    data=ROOT.RooDataSet('data','data',args)
    return ws,data,args

"""
books control distributions and loops over the events to fill the dataset for the fit
"""
def fillDataAndPlots(opt,data,ws,args,tree,baseCat='1l'):

    p4=ROOT.TLorentzVector(0,0,0,0)
    
    plots={}
    for cat in EVENTCATEGORIES:
        plots['ntrksue_%s'%cat] = ROOT.TH1F('ntrksue_%s'%cat,'; UE track multiplicity; Events',30,0,30)
        
        plots['pt_l_%s'%cat]  = ROOT.TH1F('pt_l_%s'%cat,';Lepton transverse momentum [GeV];Events',20,0,100)
        plots['y_l_%s'%cat]   = ROOT.TH1F('y_l_%s'%cat,';Lepton rapidity;Events',20,-2.1,2.1)

        plots['jsf_%s'%cat]     = ROOT.TH1F('jsf_%s'%cat,    '; JSF=M_{W}/M(jj)', 20,0,2)            
        plots['dr_jj_%s'%cat]   = ROOT.TH1F('dr_jj_%s'%cat,  '; #DeltaR(j,j); Events',15,0,6.3)
        plots['y_wjj_%s'%cat]   = ROOT.TH1F('y_wjj_%s'%cat,  '; W_{jj} rapidity;Events',20,-2.1,2.1)
        plots['pt_wjj_%s'%cat]  = ROOT.TH1F('pt_wjj_%s'%cat, '; W_{jj} transverse momentum [GeV]; Events',50,0,500)
        plots['m_wjj_%s'%cat]   = ROOT.TH1F('m_wjj_%s'%cat,  '; W_{jj} mass [GeV]; Events / 10 GeV',30,0,300)

        plots['dr_bwjj_%s'%cat]         = ROOT.TH1F('dr_bwjj_%s'%cat,        '; #DeltaR(b,jj); Events',15,0,6.3)
        plots['y_bwjj_%s'%cat]          = ROOT.TH1F('y_bwjj_%s'%cat,         '; W_{jj}+b rapidity;Events',20,-2.1,2.1)
        plots['m_bwjj_%s'%cat]          = ROOT.TH1F('m_bwjj_%s'%cat,         '; W_{jj}+b mass [GeV];Events / 10 GeV',35,50,400)
        plots['m_bwjj_vs_m_wjj_%s'%cat] = ROOT.TH2F('m_bwjj_vs_m_wjj_%s'%cat, '; W_{jj}+b mass [GeV]; W_{jj} mass [GeV]; Events',35,50,400,30,0,300)
        plots['m_bwjj_vs_m_bwl_%s'%cat] = ROOT.TH2F('m_bwjj_vs_m_bwl_%s'%cat,  '; W_{jj}+b [GeV]; W_{l}+b mass [GeV]; Events',35,50,400,30,0,300)
   
        plots['m_bwl_%s'%cat]         = ROOT.TH1F('m_bwl_%s'%cat,'; W_{l}+b mass [GeV]; Events',35,50,400)
        plots['y_bwl_%s'%cat]         = ROOT.TH1F('y_bwl_%s'%cat,'; W_{l}+b rapidity;Events',20,-2.1,2.1)

    #loop over all the events
    for i in xrange(0,tree.GetEntriesFast()):

        tree.GetEntry(i)

        #categorize jets according to the b-tag discriminator
        lightJets,bJets=[],[]
        for j in xrange(0,tree.nj):
            p4.SetPtEtaPhiM(tree.j_pt[j],tree.j_eta[j],tree.j_phi[j],tree.j_m[j])
            csvL=True if (((tree.j_btag[j]>>0) & 0x1)==1) else False
            csvM=True if (((tree.j_btag[j]>>1) & 0x1)==1) else False
            if not csvL : lightJets.append( ROOT.TLorentzVector(p4) )
            elif csvM   : bJets.append( ROOT.TLorentzVector(p4) )

        #check #jets and #b-jets
        bJetType='b'
        if len(bJets)==0:
            bJetType='q'
            #randomly assign the b-jets from among the lightJets
            random.shuffle(lightJets)
            if len(lightJets)<3 : continue
            bJets.append( lightJets.pop() )
        else :
            if len(lightJets)<2 : continue

        #assign the category
        cat=baseCat
        if len(lightJets)+len(bJets)==3 : cat += '3j'
        else : cat += '4j'
        if len(bJets)==1    : cat += '1'+bJetType
        else                : cat += '2'+bJetType

        #build the dijet system and rank by increasing dR
        dijets=[]
        for j in xrange(0,len(lightJets)):
            for k in xrange(j+1,len(lightJets)):
                drjj=lightJets[j].DeltaR(lightJets[k])
                dijets.append( (lightJets[j]+lightJets[k], drjj , (j,k) ) )
        dijets.sort(key=lambda x: x[1],reverse=False)
        wjj,drwjj,wj1,wj2=dijets[0][0],dijets[0][1],lightJets[ dijets[0][2][0] ],lightJets[ dijets[0][2][1] ]
        jsf=80.4/wjj.M()
        scaledwjj=wjj
        #scaledwjj=scaleP4(vec=wjj,scale=jsf)
        #wj1=scaleP4(vec=j1,scale=jsf)
        #wj2=scaleP4(vec=j2,scale=jsf)
        
        #build the trijet system and rank by increasing dR
        trijets=[]
        for j in xrange(0,len(bJets)):
            if j>1 : break
            drjjb=scaledwjj.DeltaR(bJets[j])
            trijets.append( (scaledwjj+bJets[j],drjjb,j) )
        trijets.sort(key=lambda x: x[1],reverse=False)
        bwjj,drbwjj,bhad=trijets[0][0],trijets[0][1],bJets[ trijets[0][2] ]

        #solve the leptonic side
        p4.SetPtEtaPhiM(tree.l_pt,tree.l_eta,tree.l_phi,0.)
        lepton=ROOT.TLorentzVector(p4)
        p4.SetPtEtaPhiM(tree.met_pt,0.,tree.met_phi,0.)
        MET=ROOT.TLorentzVector(p4)
        p4nu=getNeutrinoP4(lepton,MET)
        lepbidx = 1 if trijets[0][2]==0 else 0
        blep=bJets[lepbidx] if len(bJets)>1 else None
        bwl=p4nu+lepton+blep if blep else None

        if tree.event==271889813:
            print wjj.M(),bwjj.M(),bwl.M()
        
        #ue tracks
        nueTrks=tree.ntracks-tree.ntracks_hp

        #fill histograms
        plots['ntrksue_%s'%cat].Fill(nueTrks,tree.w)
        plots['y_l_%s'%cat].Fill(lepton.Rapidity(),tree.w)
        plots['pt_l_%s'%cat].Fill(lepton.Pt(),tree.w)
        plots['dr_jj_%s'%cat].Fill(drjj,tree.w)
        plots['y_wjj_%s'%cat].Fill(wjj.Rapidity(),tree.w)
        plots['pt_wjj_%s'%cat].Fill(wjj.Pt(),tree.w)
        plots['m_wjj_%s'%cat].Fill(wjj.M(),tree.w)
        plots['jsf_%s'%cat].Fill(jsf,tree.w)
        plots['dr_bwjj_%s'%cat].Fill(drbwjj,tree.w)
        plots['y_bwjj_%s'%cat].Fill(bwjj.Rapidity(),tree.w)
        plots['m_bwjj_%s'%cat].Fill(bwjj.M(),tree.w)
        plots['m_bwjj_vs_m_wjj_%s'%cat].Fill(bwjj.M(),wjj.M(),tree.w)
        if bwl:
            plots['m_bwjj_vs_m_bwl_%s'%cat].Fill(bwjj.M(),bwl.M(),tree.w)
            plots['m_bwl_%s'%cat].Fill(bwl.M(),tree.w)
            plots['y_bwl_%s'%cat].Fill(bwl.Rapidity(),tree.w)

        #add entry to dataset
        ws.cat('sample').setLabel(cat)
        args.setCatLabel(cat)
        argValList=[('y_l',lepton.Rapidity()),
                    ('y_bhad',bhad.Rapidity()),
                    ('y_blep',blep.Rapidity() if blep else -999.),
                    ('m_wjj',wjj.M()),
                    ('dr_jj',drwjj),
                    ('pt_wl',wjj.Pt()),
                    ('y_wl',wjj.Rapidity()),
                    ('m_bwjj',bwjj.M()),
                    ('y_bwjj',bwjj.Rapidity()),
                    ('m_bwl',bwl.M() if bwl else -999),
                    ('y_bwl',bwl.Rapidity() if bwl else -999),
                    ('dy_bwjjbwl',abs(bwl.Rapidity()-bwjj.Rapidity()) if bwl else -999),
                    ('sy_bwjjbwl',abs(bwl.Rapidity()+bwjj.Rapidity()) if bwl else -999),
                   ]
        for argName,argVal in argValList:
            args.find(argName).setVal(argVal)
        if opt.mc:
            weightVar.setVal(tree.w)
        data.add(args)

                           
    return plots


"""
display the results of the fit
"""
def showControlPlots(allPlots,fOut=None):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)

    for p in allPlots[0]:

        c.Clear()
        totalH,stackH=None,None
        for coll in allPlots:
            if not p in coll: continue
            if totalH is None:
                totalH=coll[p].Clone('%s_total'%p)
                totalH.SetDirectory(0)
                totalH.SetFillStyle(0)
                totalH.SetFillColor(0)
                totalH.SetLineWidth(2)
                totalH.SetLineColor(1)
                totalH.SetTitle('total')
                totalH.SetMarkerStyle(20)
                totalH.Reset('ICE')
            totalH.Add( coll[p] )
            if stackH is None:
                stackH=ROOT.THStack('%s_stack'%p,'')
            stackH.Add( coll[p], 'hist' )
        if totalH.Integral()==0 : continue

        is2D=True if totalH.InheritsFrom('TH2') else False

        leg = ROOT.TLegend(0.14,0.88,0.5,0.75)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)        
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.AddEntry(totalH,'total','p')
        
        if is2D:
            totalH.Draw('coltext')
        else:
            frame=totalH.Clone('frame')
            frame.Reset('ICE')
            frame.Draw()
            frame.GetYaxis().SetTitleOffset(1.4)
            frame.GetYaxis().SetRangeUser(0,totalH.GetMaximum()*1.5)
            stackH.Draw('histsame')
            totalH.Draw('same')

            for i in xrange(0,stackH.GetStack().GetEntriesFast()):
                h=stackH.GetStack().At(i)
                leg.AddEntry(h,h.GetTitle(),'f')

        leg.Draw()

        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
        label.DrawLatex(0.58,0.96,'#scale[0.8]{180 nb^{-1} (pPb at #sqrt{s}=8.16 TeV)}')
        tag='#mu+jets'
        if '1f' in p : tag ='non-iso '+tag
        if '3j' in p : tag='=3j,'
        if '4j' in p : tag='#geq4j,'
        if '1b' in p or '1q' in p : tag+='=1b'
        if '2b' in p or '2q' in p : tag+='#geq2b'
        if '1q' in p or '2q' in p : tag+='(side-band)'
        label.DrawLatex(0.68,0.9,'#scale[0.8]{#it{%s}}'%tag)
        
        c.cd()
        c.Modified()
        c.Update()

        c.SaveAs('plots/%s.png'%p)

        if fOut:
            fOut.cd()
            totalH.SetDirectory(fOut)
            totalH.Write(p)
        


"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0, type=int,   help='Verbose mode [%default]')
    parser.add_option('-t','--mc',   dest='mc',   default=False,   help='Run on MC [%default]')
    (opt, args) = parser.parse_args()

    #init the workspace
    ws,data,args=startWorkspace(opt)

    #prepare control plots and fill dataset
    allPlots=[]
    for coll in DATA:
        url,color=DATA[coll]
        fIn=ROOT.TFile(url)
        tree=fIn.Get('data')
        baseCat='1f' if 'non-iso' in coll else '1l'
        allPlots.append( fillDataAndPlots(opt,data,ws,args,tree,baseCat) )

        color=ROOT.TColor.GetColor(color)
        for p in allPlots[-1]:
            allPlots[-1][p].SetBinErrorOption(ROOT.TH1.kPoisson)
            allPlots[-1][p].SetFillStyle(1001)
            allPlots[-1][p].SetFillColor(color)
            allPlots[-1][p].SetLineColor(color)
            allPlots[-1][p].SetTitle(coll)
            #allPlots[-1][p].Sumw2()
            allPlots[-1][p].SetDirectory(0)
        
        fIn.Close()

    #import data and save workspace
    if opt.mc:
        datasetWeight = ROOT.RooDataSet(data.GetName() + "_withWeights",                                       
                                        "Weighted "+data.GetTitle(),
                                        args,
                                        ROOT.RooFit.Import(data),
                                        ROOT.RooFit.WeightVar(weightVar)                                 
                                        )                               
        datasetWeight.Print('v')

        getattr(ws,'import')(datasetWeight)
    else:
        getattr(ws,'import')(data)

    ws.writeToFile('workspace.root',True)
    
    #show results
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(False)
    fOut=ROOT.TFile.Open('plots/controlplots.root','RECREATE')
    showControlPlots(allPlots,fOut)
    fOut.Close()


if __name__ == "__main__":
    main()

