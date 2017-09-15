import ROOT
import sys
import pickle

from runDataFit import lumi
from runQCDestimation import getHistos,QCDNORM,PROCNORM
from TopLJets2015.TopAnalysis.Plot import Plot

def showPlotFor(args,var):
    
    histos=getHistos(args,var)
    for key in histos:
        for subKey in histos[key]:
            xtitle=histos[key][subKey].GetXaxis().GetTitle()
            if 'W_{jj}'  in xtitle: histos[key][subKey].GetXaxis().SetTitle("m_{jj'} [GeV]")
            if 't_{had}' in xtitle: histos[key][subKey].GetXaxis().SetTitle("M_{top} [GeV]")
            if 't_{lep}' in xtitle: histos[key][subKey].GetXaxis().SetTitle("M_{l#nub} [GeV]")
            histos[key][subKey].GetYaxis().SetTitle('Events')

    for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:


        histTotal={}
        for ch in ['e','mu']:

            qcdKey=(ch,'1l4j2q')
            key=(ch,cat)

            #scale up processes
            histos['ttbar'][key].Scale(lumi[0]*PROCNORM[('ttbar',ch)])
            histos['wjets'][key].Scale(lumi[0]*PROCNORM[('wjets',ch)])
            histos['dy'][key].Scale(lumi[0]*PROCNORM[('dy',ch)])
            qcdHisto=histos['qcd'][qcdKey].Clone('qcd_%s%s'%(ch,cat))
            qcdHisto.Scale(QCDNORM[key][0]/qcdHisto.Integral(0,qcdHisto.GetNbinsX()+1))
            for xbin in xrange(1,qcdHisto.GetNbinsX()+1):
                totalUnc=qcdHisto.GetBinError(xbin)**2+(qcdHisto.GetBinContent(xbin)*QCDNORM[key][1]/QCDNORM[key][0])**2
                qcdHisto.SetBinError(xbin,ROOT.TMath.Sqrt(totalUnc))

            if not 'ttbar' in histTotal:
                histTotal['ttbar']=histos['ttbar'][key].Clone('tot_ttbar')
                histTotal['dy']=histos['dy'][key].Clone('tot_dy')
                histTotal['wjets']=histos['wjets'][key].Clone('tot_wjets')
                histTotal['qcd']=qcdHisto.Clone('tot_qcd')
                histTotal['data']=histos['data'][key].Clone('tot_dat')
                for x in histTotal: histTotal[x].SetDirectory(0)
            else:
                histTotal['ttbar'].Add( histos['ttbar'][key] )
                histTotal['dy'].Add( histos['dy'][key] )
                histTotal['wjets'].Add( histos['wjets'][key] )
                histTotal['qcd'].Add( qcdHisto )
                histTotal['data'].Add( histos['data'][key] )

            p=Plot('%s%s_%s_control'%(ch,cat,var))
            p.add(histos['ttbar'][key],'t#bar{t}',0,False,False,False)
            p.add(histos['wjets'][key],'W+jets',ROOT.TColor.GetColor('#fee090'),False,False,False)
            p.add(histos['dy'][key],'DY',ROOT.TColor.GetColor('#fc8d59'),False,False,False)
            p.add(qcdHisto,'Multijets',ROOT.TColor.GetColor('#e0f3f8'),False,False,False)
            p.add(histos['data'][key],'Data',1,True,False,False)
        
            p.ratiorange=(0.4,1.7)
            p.com='8.16 TeV'
            p.show(outDir='./',lumi=174*1e-3)
            p.reset()

        #combined plots
        p=Plot('l%s_%s_control'%(cat,var))
        p.add(histTotal['ttbar'],'t#bar{t}',ROOT.TColor.GetColor('#d7191c'),False,False,False)
        p.add(histTotal['wjets'],'W+jets',ROOT.TColor.GetColor('#3ed2e0'),False,False,False)
        p.add(histTotal['dy'],'DY',ROOT.TColor.GetColor('#2baeba'),False,False,False)
        p.add(histTotal['qcd'],'Multijets',ROOT.TColor.GetColor('#2b83ba'),False,False,False)
        p.add(histTotal['data'],'Data',1,True,False,False)
        p.cmsLabel='#scale[1.2]{#bf{CMS}} #scale[0.75]{#it{Supplementary}}'
        p.ratiorange=(0.4,1.7)
        p.ratioTitle='#scale[0.9]{Data/Exp.}'
        p.legSize=0.05
        p.com='#sqrt{s_{NN}} = 8.16 TeV'
        tag='e^{#pm} / #mu^{#pm} + #geq4j'
        if '1f'   in cat: tag ='non-iso '+tag
        if '2q'   in cat: tag+=' (=0b)'
        if '1b1q' in cat: tag+=' (=1b)'
        if '2b'   in cat: tag+=' (#geq2b)'
        tag='#scale[1.1]{#bf{%s}}'%tag
        p.show(outDir='./',lumi=174*1e-3,extraText=tag)
        p.reset()

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #load qcd normalization
    global QCDNORM
    with open('qcdnorm.pck','r') as fIn:
        QCDNORM=pickle.load(fIn)

    args=sys.argv[1:]
    for var in ['met','mtw','minmlb','maxmlb','mjj','drjj','ntracks','ntrackshp','mthad','mtlep']: 
        showPlotFor(args,var)


if __name__ == "__main__":
    main()
