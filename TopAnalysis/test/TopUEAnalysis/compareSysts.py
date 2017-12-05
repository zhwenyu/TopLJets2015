import ROOT

from os import listdir
from os.path import isfile, join
import re
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from TopLJets2015.TopAnalysis.Plot import Plot


def getMC(ana='chmult/inc'):

    histos=[]
    fIn=ROOT.TFile.Open('store/TOP-17-015/%s/MC13TeV_TTJets.root'%ana)

    systs=[('nominal','_0',1,True,False),
           #('SF_{vtx.} (#mu) ','_22',ROOT.kGray,False,False),
           ('SF_{#eta} (from #mu)', '_25',ROOT.kRed,False,False),
           ('SF_{#eta} (from D*)','_26',ROOT.kCyan,False,False)] 

    p=Plot('chmult_tkeffsyst')
    p.ratiorange=(0.7,1.3)
    p.doPoissonErrorBars=False
    p.cmsLabel='#bf{CMS} #it{simulation}'
    p.doMCOverData = False
    p.ratioTitle='Ratio '
    for s,pf,color,isref,spimpose in systs:
        h=fIn.Get('reco'+pf)
        try :
            h.Scale(1./h.Integral())
        except:
            print 'Skipping ',s
            continue
        
        p.add(h,s,color,isref,spimpose,False)
        histos.append(h.Clone())
        histos[-1].SetDirectory(0)

    p.show(outDir=' ~/www/TOP-17-015',lumi=35900,noStack=True,noRatio=False)

    fIn.Close()


def getData(ana='chmult/inc'):

    #list files
    path='store/TOP-17-015/'+ana
    onlyfiles = [f for f in listdir(path) if 'MuonEG' in f]

    dist='reco_0'
    
    plots={}
    total=None
    for f in onlyfiles:
        era=re.search('(?<=2016)\w+', f).group(0)[0]
        era='BCDEF' if era in ['B','C','D','E','F']  else 'GH'

        #getdistribution from file
        fIn=ROOT.TFile.Open(join(path,f))
        h=fIn.Get(dist)
        if not era in plots:
            plots[era]=h.Clone(era)
            plots[era].SetTitle(era)
            plots[era].SetDirectory(0)
        else:
            plots[era].Add(h)
        if total is None:
            total=h.Clone('total')
            total.SetTitle('2016')
            total.SetDirectory(0)
        else:
            total.Add(h)
        fIn.Close()


        p=Plot('%s_eras'%ana.replace('/','_'))
        colors={'B':ROOT.kRed+1,
                'C':ROOT.kRed+2,
                'D':ROOT.kRed+3,
                'E':ROOT.kRed+4,
                'F':ROOT.kRed+5,
                'BCDEF':ROOT.kRed,
                'G':ROOT.kCyan-1,
                'H':ROOT.kCyan-2,
                'GH':ROOT.kCyan,
                }
        for era in plots:
            plots[era].Scale(1./plots[era].Integral())
            p.add(plots[era],era,colors[era],False,False,False)
        total.Scale(1./total.Integral())
        p.add(total,'2016',1,True,True,False)
        p.ratiorange=(0.7,1.3)
        p.cmsLabel='#bf{CMS} #it{preliminary}'
        p.doPoissonErrorBars=False
        p.ratioTitle='Ratio'
        p.show(outDir=' ~/www/TOP-17-015/',lumi=35900,noStack=True,noRatio=False) 

        return plots


mcHistos=getMC()
dataHistos=getData()    
