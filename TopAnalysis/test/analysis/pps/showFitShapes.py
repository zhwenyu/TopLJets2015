import sys
import ROOT
import re
import os
from TopLJets2015.TopAnalysis.Plot import *

def showShapes(resultsDir,name,title,mass,boson,lumi,plotData=True,showPseudoData=True,showAllBkgs=True,subCatTitles=None):

    """ show the shapes corresponding to a given analysis """
    
    shapeFiles=[os.path.join(resultsDir,f) for f in os.listdir(resultsDir) if 'shapes_' in f and '.root' in f]
    for f in shapeFiles:

        fIn=ROOT.TFile.Open(f)

        #check consistency of the file with the boson to plot
        if 'shapes_22' in f and boson=='g':
            v='g'
            channel='#gamma'
        elif 'shapes_169' in f and boson=='z': 
            v='zmm'
            channel='Z#rightarrow#mu#mu'
        elif 'shapes_121' in f and boson=='z': 
            v='zee'
            channel='Z#rightarrowee'
        else:
            continue
        
        ncats=len(subCatTitles) if subCatTitles else 1
        for icat in range(ncats):
            bkgH       = fIn.Get('bkg_%s'%v)
            fidsigH    = fIn.Get('fidsig_%s_m%d'%(v,mass))
            outfidsigH = fIn.Get('outfidsig_%s_m%d'%(v,mass))
            dataH      = fIn.Get('data_obs_%s'%(v))

            if showPseudoData:
                dataH.Reset('ICE')
                nexp=bkgH.Integral()
                if fidsigH: nexp+=bkgH.Integral()
                if outfidsigH: nexp+=outfidsigH.Integral()
                for iev in range( ROOT.gRandom.Poisson( nexp ) ):
                    dataH.Fill( bkgH.GetRandom() )
            try:
                
                p=Plot('%s_%s_inc'%(name,v))
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                if fidsigH:
                    p.add(fidsigH,            title='fiducial #scale[0.8]{(%d)}'%mass, color=ROOT.TColor.GetColor('#fdc086'), isData=False, spImpose=False, isSyst=False)
    
                if showAllBkgs:
                    if outfidsigH:
                        p.add(outfidsigH,         title='non-fiducial',  color=ROOT.TColor.GetColor('#a6cee3'), isData=False, spImpose=False, isSyst=False)
                    p.add(bkgH,               title='background',    color=ROOT.TColor.GetColor('#1f78b4'), isData=False, spImpose=False, isSyst=False)
                else:
                    allBkg=bkgH.Clone('allbkg')
                    if outfidsigH : 
                        allBkg.Add(outfidsigH)
                    p.add(allBkg,               title='background',    color=ROOT.TColor.GetColor('#1f78b4'), isData=False, spImpose=False, isSyst=False)

                if plotData:
                    dtitle='pseudo-data' if showPseudoData else 'Data'
                    p.add(dataH, title=dtitle,   color=1, isData=True, spImpose=False, isSyst=False)

                if fidsigH:
                    p.add(fidsigH.Clone(),    title=title+'#scale[0.8]{(%d)}'%mass, color=ROOT.TColor.GetColor('#fdc086'), isData=False, spImpose=True,  isSyst=False)

                if showAllBkgs and outfidsigH:
                    p.add(outfidsigH.Clone(), title='non-fiducial',  color=ROOT.TColor.GetColor('#a6cee3'), isData=False, spImpose=True,  isSyst=False)
                p.ratiorange=[0.68,1.32]

                extraText=[subCatTitles[icat] if subCatTitles else '']
                p.show('./',lumi*1000,extraText='\\'.join(extraText))

                colors=[ROOT.kGreen+1,ROOT.kAzure+3,ROOT.kRed+2,ROOT.kOrange+2]

                #background systs
                p=Plot('%s_%s_inc_bkgunc'%(name,v))
                p.noErrorsOnRatio=True
                p.doPoissonErrorBars=False
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                p.add(bkgH, title='background', color=1, isData=True,spImpose=False, isSyst=False)
                ic=0
                for syst,title in [('Up',            'e#mu mix'),
                                   ('SingleDiffUp',  'single arm mix')]:
                    h=fIn.Get('bkg_%s_bkgShape%s'%(v,syst))
                    p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                    ic+=1
                p.ratiorange=[0.76,1.24]                
                p.show('./',lumi,noStack=True,extraText=extraText)

                #signal systs
                if fidsigH:
                    p=Plot('%s_%s_sigunc'%(name,v))
                    #fidsigH.Scale(1./5.)
                    p.doPoissonErrorBars=False
                    p.noErrorsOnRatio=True
                    p.xtit='Missing mass [GeV]'
                    p.ytit='Events'
                    p.add(fidsigH, title='signal', color=1, isData=True,spImpose=False, isSyst=False)
                    ic=0
                    for syst,title in [('ShapeUp',   'e#mu mix.'),
                                       ('CalibUp',   'Time dependence'),
                                       ('PPSEffUp',  'Strip. efficiency'),
                                       ('PzModelUp', 'p_{z}(pp)')]:
                        h=fIn.Get('fidsig_%s_m%d_sig%s'%(v,mass,syst))
                        print ic,syst,title,h
                        p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                        ic+=1
                    p.ratiorange=[0.5,1.46]
                    p.show('./',lumi,noStack=True,extraText=extraText)

                #out fiducial signal systs
                if outfidsigH:
                    p=Plot('%s_%s_outfidsigunc'%(name,v))
                    #outfidsigH.Scale(1./5.)
                    p.doPoissonErrorBars=False
                    p.noErrorsOnRatio=True
                    p.xtit='Missing mass [GeV]'
                    p.ytit='Events'
                    p.add(outfidsigH, title='out-fid. signal', color=1, isData=True,spImpose=False, isSyst=False)
                    ic=0
                    for syst,title in [('ShapeUp',   'e#mu mix.'),
                                       ('CalibUp',   'Time dependence'),
                                       ('PPSEffUp',  'Strip. efficiency'),
                                       ('PzModelUp', 'p_{z}(pp)')]:
                        h=fIn.Get('outfidsig_%s_m%d_sig%s'%(v,mass,syst))
                        p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                        ic+=1

                    p.ratiorange=[0.5,1.46]                    
                    p.show('./',lumi,noStack=True,extraText=extraText)

            except Exception as e:
                print e
                pass

        fIn.Close()


def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    lumi=37.5
    bosonName=sys.argv[3]
    if bosonName=='g' : 
        bosonName='#gamma'
        lumi=2.76

    from prepareOptimScanCards import OPTIMLIST
    optimPt=int(re.search('optim_(\d+)',sys.argv[1]).group(1))-1
    cuts=OPTIMLIST[optimPt][2].split(',')

    showShapes(resultsDir=sys.argv[1],
               name='shapes',
               title='pp%sX'%bosonName,
               mass=int(sys.argv[2]),
               boson=sys.argv[3],
               lumi=lumi,
               plotData=False, #True,
               showPseudoData=False,
               showAllBkgs=True,  #False)
               subCatTitles=None)

if __name__ == "__main__":
    main()
