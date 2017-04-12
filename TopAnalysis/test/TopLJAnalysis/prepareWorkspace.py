#!/usr/bin/env python2.7

import ROOT
import optparse
import os,sys
import random
from collections import OrderedDict
from kinTools import *
from datasets import *

EVENTCATEGORIES=[
    '1l4j2b','1l4j1b','1l4j1q','1f4j1q','1f4j1b','1f4j2b',
    '1l3j1b','1l3j1q','1f3j1q','1f3j1b'
    ]


"""
"""
def getMCTruthCat(tree,thad):

    mcTruthCat='ww'

    p4=ROOT.TLorentzVector(0,0,0,0)
    allJets=[]
    for j in xrange(0,tree.ngj):
        p4.SetPtEtaPhiM(tree.gj_pt[j],tree.gj_eta[j],tree.gj_phi[j],tree.gj_m[j])
        allJets.append( ROOT.TLorentzVector(p4) )

    #require a gen-dijet peaking at the W mass
    dijets=[]
    for i in xrange(0,len(allJets)):
        for j in xrange(i+1,len(allJets)):
            dijets.append( (allJets[i]+allJets[j],
                           ((allJets[i]+allJets[j]).M()-80.4)**2 ,
                           (i,j) ) )
    dijets.sort(key=lambda x: x[1])
    if len(dijets)==0 : return mcTruthCat
    wjj,wk,ij=dijets[0]
    if ROOT.TMath.Sqrt(wk)>10: return mcTruthCat

    #require a gen-trijet peaking at the top mass
    trijets=[]
    for i in xrange(0,len(allJets)):
        if i in ij : continue
        trijets.append( ( (allJets[i]+wjj),
                        ((allJets[i]+wjj).M()-172.5)**2,
                        (i,ij[0],ij[1]) ) )
    if len(trijets)==0 : return mcTruthCat
    bwjj,tk,ijk=trijets[0]
    if ROOT.TMath.Sqrt(tk)>10: return mcTruthCat

    #match geometrically to reco level
    print ijk,allJets[ijk[1]].DeltaR(thad[1]),allJets[ijk[2]].DeltaR(thad[2]),allJets[ijk[0]].DeltaR(thad[0])
    if allJets[ijk[1]].DeltaR(thad[1])<0.5:
        if allJets[ijk[2]].DeltaR(thad[2])<0.5 :
            mcTruthCat='cw'
            if allJets[ijk[0]].DeltaR(thad[0])<0.5:
                mcTruthCat='cc'


    print bwjj.M(),wjj.M(),"->",(thad[0]+thad[1]+thad[2]).M(),(thad[1]+thad[2]).M(),mcTruthCat
    raw_input()
    return mcTruthCat





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
             'pt_wjj[0,500]',    'y_wjj[-4,4]',                    #w had level
             'm_bwjj[0,500]',    'y_bwjj[-4,4]',                   #hadronic top level
             'm_bwl[0,500]',     'y_bwl[-4,4]',                    #leptonic top level
             'dy_bwjjbwl[0,5]',  'sy_bwjjbwl[0,5]'                 #top vs anti-top
             ]
    for var in varList :
        args.add( ws.factory(var) )

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
def fillDataAndPlots(data,ws,args,tree,leptonSel,mcTruth=False,orderByMass=True):

    p4=ROOT.TLorentzVector(0,0,0,0)

    plots={}
    baseCat='1l'
    for lid in leptonSel:
        if lid>1000 :
            baseCat='1f'
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
        plots['m_bwjj_%s'%cat]          = ROOT.TH1F('m_bwjj_%s'%cat,         '; W_{jj}+b mass [GeV];Events / 10 GeV',40,50,450)
        plots['m_bwjj_vs_m_wjj_%s'%cat] = ROOT.TH2F('m_bwjj_vs_m_wjj_%s'%cat, '; W_{jj}+b mass [GeV]; W_{jj} mass [GeV]; Events',40,50,450,30,0,300)
        plots['m_bwjj_vs_m_bwl_%s'%cat] = ROOT.TH2F('m_bwjj_vs_m_bwl_%s'%cat,  '; W_{jj}+b [GeV]; W_{l}+b mass [GeV]; Events',40,50,450,30,0,300)

        plots['m_bwl_%s'%cat]         = ROOT.TH1F('m_bwl_%s'%cat,'; W_{l}+b mass [GeV]; Events',40,50,450)
        plots['y_bwl_%s'%cat]         = ROOT.TH1F('y_bwl_%s'%cat,'; W_{l}+b rapidity;Events',20,-2.1,2.1)

    #add MC truth if required
    if mcTruth:
        for p in plots.keys():
            for truth in ['cc','cw','ww']:
                newp='%s_%s'%(p,truth)
                plots[newp]=plots[p].Clone(newp)

    #loop over all the events
    totalEntries=tree.GetEntries()
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)
        if i%100==0 : sys.stdout.write('\r [ %d/%d ] done' %(i,totalEntries))

        if tree.l_pt<30 or abs(tree.l_eta)>2.1: continue
        #if tree.l_pt<30 or abs(tree.l_eta)>1.45: continue
        if len(leptonSel)!=0 and not tree.l_id in leptonSel : continue

        #categorize jets according to the b-tag discriminator
        lightJets,loosebJets,bJets=[],[],[]
        for j in xrange(0,tree.nj):
            p4.SetPtEtaPhiM(tree.j_pt[j],tree.j_eta[j],tree.j_phi[j],tree.j_m[j])
            if p4.Pt()<25 or abs(p4.Eta())>2.5 : continue
            csvL=True if (((tree.j_btag[j]>>0) & 0x1)==1) else False
            csvM=True if (((tree.j_btag[j]>>1) & 0x1)==1) else False
            if not csvL            : lightJets.append( ROOT.TLorentzVector(p4) )
            elif csvL and not csvM : loosebJets.append( ROOT.TLorentzVector(p4) )
            else                   : bJets.append( ROOT.TLorentzVector(p4) )
        lightJets.sort(key=lambda x: x.Pt())
        loosebJets.sort(key=lambda x: x.Pt())
        bJets.sort(key=lambda x: x.Pt())

        #build the dijet system and rank by increasing dR
        dijets=buildDijets(lightJets,orderByMass)
        if len(dijets)==0 : continue
        wjj,drjj,wj1,wj2=dijets[0][0],dijets[0][1],lightJets[ dijets[0][2][0] ],lightJets[ dijets[0][2][1] ]
        if wjj.M()<1 :
            print i,tree.nj
            print dijets
        jsf=80.4/wjj.M()


        #build the trijet system
        bJetType='b'
        trijets=buildTrijets(scaleP4(wj1,jsf),scaleP4(wj2,jsf),bJets,orderByMass)
        if len(trijets)==0:

            bJetType='q'
            if baseCat!='1f' and len(loosebJets)!=0:
                bJets.append( loosebJets.pop(random.randint(0,len(loosebJets)-1) ) )
            else:
                #require at least 3 jets for one to be popped as a b-jet candidate
                if len(lightJets)<3: continue
                bJets.append( lightJets.pop( random.randint(0,len(lightJets)-1) ) )

                #re-build the light jets
                dijets=buildDijets(lightJets,orderByMass)
                wjj,drjj,wj1,wj2=dijets[0][0],dijets[0][1],lightJets[ dijets[0][2][0] ],lightJets[ dijets[0][2][1] ]

            #rebuild the trijets
            #trijets=buildTrijets(wj1,wj2,bJets,orderByMass)
            trijets=buildTrijets(scaleP4(wj1,jsf),scaleP4(wj2,jsf),bJets,orderByMass)

        #if no trijet candidate found, reject event
        if len(trijets)==0 : continue

        #get the objects
        bwjj,drbwjj,bhad=trijets[0][0],trijets[0][1],bJets[ trijets[0][2] ]

        #assign the category
        cat=baseCat
        if len(lightJets)+len(bJets)==3 : cat += '3j'
        else                            : cat += '4j'
        if len(bJets)==1                : cat += '1'+bJetType
        else                            : cat += '2'+bJetType

        #solve the leptonic side
        p4.SetPtEtaPhiM(tree.l_pt,tree.l_eta,tree.l_phi,0.)
        lepton=ROOT.TLorentzVector(p4)
        p4.SetPtEtaPhiM(tree.met_pt,0.,tree.met_phi,0.)
        MET=ROOT.TLorentzVector(p4)
        p4nu=getNeutrinoP4(lepton,MET)
        wl=lepton+p4nu
        lepbidx = 1 if trijets[0][2]==0 else 0
        blep=bJets[lepbidx] if len(bJets)>1 else None
        bwl=wl+blep if blep else None

        #ue tracks
        nueTrks=tree.ntracks-tree.ntracks_hp

        #fill histograms
        weight=tree.w
        catsToFill=[cat]
        if mcTruth:
                mcTruthCat=getMCTruthCat(tree=tree,thad=(bhad,wj1,wj2))
                catsToFill.append( cat+'_'+mcTruthCat )
        for c in catsToFill:
            plots['ntrksue_%s'%c].Fill(nueTrks,weight)
            plots['y_l_%s'%c].Fill(lepton.Rapidity(),weight)
            plots['pt_l_%s'%c].Fill(lepton.Pt(),weight)
            plots['dr_jj_%s'%c].Fill(drjj,weight)
            plots['y_wjj_%s'%c].Fill(wjj.Rapidity(),weight)
            plots['pt_wjj_%s'%c].Fill(wjj.Pt(),weight)
            plots['m_wjj_%s'%c].Fill(wjj.M(),weight)
            plots['jsf_%s'%c].Fill(jsf,weight)
            plots['dr_bwjj_%s'%c].Fill(drbwjj,weight)
            plots['y_bwjj_%s'%c].Fill(bwjj.Rapidity(),weight)
            plots['m_bwjj_%s'%c].Fill(bwjj.M(),weight)
            plots['m_bwjj_vs_m_wjj_%s'%c].Fill(bwjj.M(),wjj.M(),weight)
            if bwl:
                plots['m_bwjj_vs_m_bwl_%s'%c].Fill(bwjj.M(),bwl.M(),weight)
                plots['m_bwl_%s'%c].Fill(bwl.M(),weight)
                plots['y_bwl_%s'%c].Fill(bwl.Rapidity(),weight)

        #add entry to dataset
        ws.cat('sample').setLabel(cat)
        args.setCatLabel(cat)
        argValList=[('y_l',lepton.Rapidity()),
                    ('y_bhad',bhad.Rapidity()),
                    ('y_blep',blep.Rapidity() if blep else -999.),
                    ('m_wjj',wjj.M()),
                    ('dr_jj',drjj),
                    ('pt_wjj',wjj.Pt()),
                    ('y_wjj',wjj.Rapidity()),
                    ('pt_wl',wl.Pt()),
                    ('y_wl',wl.Rapidity()),
                    ('m_bwjj',bwjj.M()),
                    ('y_bwjj',bwjj.Rapidity()),
                    ('m_bwl',bwl.M() if bwl else -999),
                    ('y_bwl',bwl.Rapidity() if bwl else -999),
                    ('dy_bwjjbwl',abs(bwl.Rapidity()-bwjj.Rapidity()) if bwl else -999),
                    ('sy_bwjjbwl',abs(bwl.Rapidity()+bwjj.Rapidity()) if bwl else -999),
                   ]
        for argName,argVal in argValList:
            args.find(argName).setVal(argVal)
        data.add(args)

    return plots


"""
display the results of the fit
"""
def showControlPlots(allPlots,fOut=None,outDir='plots/'):

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
            if coll[p].Integral()==0 : continue
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
        if '_pp' in outDir:
            label.DrawLatex(0.58,0.96,'#scale[0.8]{881 pb^{-1} (pp at #sqrt{s}=8.0 TeV)}')
        else:
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
        c.SaveAs('%s/%s.png'%(outDir,p))

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
    parser.add_option('-d', '--data',         dest='data',        default='Data8TeV_pp',   type='string',          help='dataset to use [%default]')
    parser.add_option('-v', '--verbose',      dest='verbose',     default=0,               type=int,               help='Verbose mode [%default]')
    parser.add_option(      '--orderByMass',  dest='orderByMass', default=False,           action='store_true',    help='Verbose mode [%default]')

    (opt, args) = parser.parse_args()

    #init the workspace
    ws,data,args=startWorkspace(opt)

    #prepare control plots and fill dataset
    allPlots=[]
    dataset=defineDataset(data=opt.data)
    for coll in dataset:
        urlList,color,leptonSel,mcTruth=dataset[coll]
        tree=ROOT.TChain('data')
        for url in urlList: tree.AddFile(url)
        allPlots.append( fillDataAndPlots(data,ws,args,tree,leptonSel,mcTruth,opt.orderByMass) )

        color=ROOT.TColor.GetColor(color)
        for p in allPlots[-1]:
            if 'Data' in opt.data:
                allPlots[-1][p].SetBinErrorOption(ROOT.TH1.kPoisson)
            allPlots[-1][p].SetFillStyle(1001)
            allPlots[-1][p].SetFillColor(color)
            allPlots[-1][p].SetLineColor(color)
            allPlots[-1][p].SetTitle(coll)
            allPlots[-1][p].Sumw2()
            allPlots[-1][p].SetDirectory(0)


    #import data and save workspace
    getattr(ws,'import')(data)
    ws.writeToFile('workspace_%s.root'%opt.data,True)

    #show results
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(False)
    outDir='plots/%s'%opt.data
    os.system('mkdir -p %s'%outDir)
    fOut=ROOT.TFile.Open('%s/controlplots.root'%outDir,'RECREATE')
    showControlPlots(allPlots,fOut,outDir=outDir)
    fOut.Close()


if __name__ == "__main__":
    main()
