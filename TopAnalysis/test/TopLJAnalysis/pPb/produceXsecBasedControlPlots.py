import ROOT
import sys
import os
import pickle
import optparse

from runDataFit import lumi
from runQCDestimation import getHistos,QCDNORM,PROCNORM
from TopLJets2015.TopAnalysis.Plot import Plot

def showPlotFor(args,var,opt):
    
    histos=getHistos(args,var)
    for key in histos:
        for subKey in histos[key]:
            xtitle=histos[key][subKey].GetXaxis().GetTitle()
            if 'W_{jj}'  in xtitle: histos[key][subKey].GetXaxis().SetTitle("m_{jj'} [GeV]")
            if 't_{had}' in xtitle: histos[key][subKey].GetXaxis().SetTitle("m_{top} [GeV]")
            if 't_{lep}' in xtitle: histos[key][subKey].GetXaxis().SetTitle("m_{l#nub} [GeV]")
            histos[key][subKey].GetYaxis().SetTitle('Events')

    for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:


        histTotal={}
        for ch in ['e','mu']:

            qcdKey=(ch,'1l4j2q')
            key=(ch,cat)

            #scale up processes
            histos['ttbar'][key].Scale(opt.lumi*PROCNORM[('ttbar',ch)])
            histos['wjets'][key].Scale(opt.lumi*PROCNORM[('wjets',ch)])
            histos['dy'][key].Scale(opt.lumi*PROCNORM[('dy',ch)])
            qcdHisto=histos['qcd'][qcdKey].Clone('qcd_%s%s'%(ch,cat))
            qcdHisto.Scale((opt.lumi/lumi[0])*QCDNORM[key][0]/qcdHisto.Integral(0,qcdHisto.GetNbinsX()+1))
            for xbin in xrange(1,qcdHisto.GetNbinsX()+1):
                totalUnc=qcdHisto.GetBinError(xbin)**2+(qcdHisto.GetBinContent(xbin)*QCDNORM[key][1]/QCDNORM[key][0])**2
                qcdHisto.SetBinError(xbin,ROOT.TMath.Sqrt(totalUnc))

            p=Plot('%s%s_%s_control'%(ch,cat,var))
            p.add(histos['ttbar'][key],'t#bar{t}',0,False,False,False)
            p.add(histos['wjets'][key],'W+jets',ROOT.TColor.GetColor('#fee090'),False,False,False)
            p.add(histos['dy'][key],'DY',ROOT.TColor.GetColor('#fc8d59'),False,False,False)
            p.add(qcdHisto,'Multijets',ROOT.TColor.GetColor('#e0f3f8'),False,False,False)

            if opt.pseudoData:
                print 'Generating pseudo-data from total'
                histos['data'][key].Reset('ICE')
                for m in p.mc:
                    nevts=p.mc[m].Integral()
                    for i in xrange(0,ROOT.gRandom.Poisson(nevts)):
                        histos['data'][key].Fill( p.mc[m].GetRandom() )                        
                p.add(histos['data'][key],'Pseudo data',1,True,False,False)
                p.cmsLabel='#scale[0.9]{#bf{CMS} #it{simulation preliminary}}'
            else:
                p.add(histos['data'][key],'Data',1,True,False,False)

            #create the totals
            for pkey in ['ttbar','dy','wjets','qcd','data']:
                if not pkey in histTotal:
                    if pkey=='qcd':
                        histTotal[pkey]=qcdHisto.Clone('tot_qcd')
                    else:
                        histTotal[pkey]=histos[pkey][key].Clone('tot_%s'%pkey)
                    histTotal[pkey].SetDirectory(0)
                    histTotal[pkey].Reset('ICE')     
                if pkey=='qcd':               
                    histTotal[pkey].Add( qcdHisto )
                else:
                    histTotal[pkey].Add( histos[pkey][key] )

            p.ratioFrameDrawOpt='pX0'
            p.doMCOverData=False
            p.ratiorange=(0.4,1.7)
            p.com='8.16 TeV'
            p.show(outDir='./',lumi=opt.lumi*1e-3)
            p.appendTo('xsec_plotter.root')
            p.reset()
            
        #combined plots
        p=Plot('l%s_%s_control'%(cat,var))
        p.add(histTotal['ttbar'],'t#bar{t}',ROOT.TColor.GetColor('#d7191c'),False,False,False)
        p.add(histTotal['wjets'],'W+jets',ROOT.TColor.GetColor('#3ed2e0'),False,False,False)
        p.add(histTotal['dy'],'DY',ROOT.TColor.GetColor('#2baeba'),False,False,False)
        p.add(histTotal['qcd'],'Multijets',ROOT.TColor.GetColor('#2b83ba'),False,False,False)
        p.add(histTotal['data'],'Data',1,True,False,False)
        if opt.pseudoData:
            p.cmsLabel='#scale[0.9]{#bf{CMS} #it{simulation preliminary}}'            
        else:
            p.cmsLabel='#scale[1.2]{#bf{CMS}} #scale[0.75]{#it{Supplementary}}'
        p.ratiorange=(0.4,1.7)
        p.ratioTitle='#scale[0.9]{Data/Exp.}'
        p.legSize=0.05
        p.ratioFrameDrawOpt='pX0'
        p.doMCOverData=False
        p.com='#sqrt{s_{NN}} = 8.16 TeV'
        tag='e^{#pm} / #mu^{#pm} + #geq4j'
        if '1f'   in cat: tag ='non-iso '+tag
        if '2q'   in cat: tag+=' (=0b)'
        if '1b1q' in cat: tag+=' (=1b)'
        if '2b'   in cat: tag+=' (#geq2b)'
        tag='#scale[1.1]{#bf{%s}}'%tag
        p.show(outDir='./',lumi=opt.lumi*1e-3,extraText=tag)
        p.appendTo('xsec_plotter.root')
        p.reset()

"""
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-l', '--lumi',       dest='lumi',       help='luminosity [nb] [%default]',   default=lumi[0], type=float)
    parser.add_option(      '--pseudoData', dest='pseudoData', help='pseudo-data [%default]',       default=False,   action='store_true')
    parser.add_option('-b', '--batch',      dest='batch',      help='run in batch mode [%default]', default=False,   action='store_true')
    (opt, args) = parser.parse_args()

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(opt.batch)
    os.system('rm xsec_plotter.root')

    #load qcd normalization
    global QCDNORM
    with open('qcdnorm.pck','r') as fIn: QCDNORM=pickle.load(fIn)

    #do plots
    for var in ['met','mtw','minmlb','maxmlb','mjj','pt_l','y_l','drjj','ntracks','ntrackshp','mthad','mtlep']: 
        print args
        showPlotFor(args,var,opt)


if __name__ == "__main__":
    main()
