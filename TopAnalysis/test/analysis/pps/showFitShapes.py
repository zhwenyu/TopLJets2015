import sys
import ROOT
import re
import os
from TopLJets2015.TopAnalysis.Plot import *

def showShapes(resultsDir,name,plotTitle,mass,boson,lumi,plotData=True,showPseudoData=True,showAllBkgs=True,plotpfix=''):

    """ show the shapes corresponding to a given analysis """

    colors=[ROOT.kGreen+1,ROOT.kAzure+3,ROOT.kOrange+2,ROOT.kGray,ROOT.kRed+2]
    
    shapeFiles=[]
    for f in os.listdir(resultsDir):
        if not '.root' in f : continue
        if len(f.split('_'))!=2: continue
        shapeFiles.append(os.path.join(resultsDir,f))
    print shapeFiles

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
            p=Plot('%s_%s_inc%s'%(name,v,plotpfix))
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
                p.add(fidsigH.Clone(),    title=plotTitle+'#scale[0.8]{(%d)}'%mass, color=ROOT.TColor.GetColor('#fdc086'), isData=False, spImpose=True,  isSyst=False)

            if showAllBkgs and outfidsigH:
                p.add(outfidsigH.Clone(), title='non-fiducial',  color=ROOT.TColor.GetColor('#a6cee3'), isData=False, spImpose=True,  isSyst=False)
            
            p.ratiorange=[0.68,1.32]
            p.show('./',lumi*1000,extraText=plotTitle)

            

            #background systs
            p=Plot('%s_%s_inc_bkgunc%s'%(name,v,plotpfix))
            p.noErrorsOnRatio=True
            p.doPoissonErrorBars=False
            p.xtit='Missing mass [GeV]'
            p.ytit='Events'
            p.add(bkgH, title='background', color=1, isData=True,spImpose=False, isSyst=False)
            ic=0
            for syst,title in [('bkgShapeEMUp',            'e#mu mix'),
                               ('bkgShapeSingleDiffUp',  'single arm mix')]:
                h=fIn.Get('bkg_%s_%s'%(v,syst))
                p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                ic+=1
            p.ratiorange=[0.76,1.24]                
            p.show('./',lumi,noStack=True,extraText=plotTitle)

            #signal systs
            if fidsigH:
                p=Plot('%s_%s_sigunc%s'%(name,v,plotpfix))
                #fidsigH.Scale(1./5.)
                p.doPoissonErrorBars=False
                p.noErrorsOnRatio=True
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                p.add(fidsigH, title='signal', color=1, isData=True,spImpose=False, isSyst=False)
                ic=0
                for syst,title in [('ShapeEMUp',   'e#mu mix.'),
                                   ('CalibUp',   'Time dependence'),                     
                                   ('PPSEffUp',  'PPS efficiency'),
                                   ('PzModelUp', 'p_{z}(pp)')]:
                    h=fIn.Get('fidsig_%s_m%d_sig%s'%(v,mass,syst))
                    print ic,syst,title,h
                    p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                    ic+=1
                p.ratiorange=[0.5,1.46]
                p.show('./',lumi,noStack=True,extraText=plotTitle)

            #out fiducial signal systs
            if outfidsigH:
                p=Plot('%s_%s_outfidsigunc%s'%(name,v,plotpfix))
                #outfidsigH.Scale(1./5.)
                p.doPoissonErrorBars=False
                p.noErrorsOnRatio=True
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                p.add(outfidsigH, title='out-fid. signal', color=1, isData=True,spImpose=False, isSyst=False)
                ic=0
                for syst,title in [('ShapeEMUp',   'e#mu mix.'),
                                   ('CalibUp',   'Time dependence'),
                                   ('PPSEffUp',  'PPS efficiency'),
                                   ('PzModelUp', 'p_{z}(pp)')]:
                    h=fIn.Get('outfidsig_%s_m%d_sig%s'%(v,mass,syst))
                    p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                    ic+=1

                p.ratiorange=[0.5,1.46]                    
                p.show('./',lumi,noStack=True,extraText=plotTitle)

        except Exception as e:
            print e
            pass

        fIn.Close()


def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    from prepareOptimScanCards import OPTIMLIST
    optimPt=int(re.search('optim_(\d+)',sys.argv[1]).group(1))-1
    cuts=OPTIMLIST[optimPt][2].split(',')

    lumi=37.5
    bosonName=sys.argv[3]
    bosonName=bosonName.replace('mm','#mu#mu')
    bosonName=bosonName.replace('z','Z')
    if bosonName=='g' : 
        bosonName='#gamma'
        lumi=2.76
    plotpfix=sys.argv[4] if len(sys.argv)>4 else ''

    title='pp%sX'%bosonName
    if plotpfix=='mm':
        title+= ' (multi-multi)'
    if plotpfix=='ms':
        title+= ' (multi-single)'
    if plotpfix=='sm':
        title+= ' (single-multi)'
    if plotpfix=='ss':
        title+= ' (single-single)'

    showShapes(resultsDir=sys.argv[1],
               name='shapes',
               plotTitle=title,
               mass=int(sys.argv[2]),
               boson=sys.argv[3],
               lumi=lumi,
               plotData=False,
               showPseudoData=False,
               showAllBkgs=True, 
               plotpfix=plotpfix)

if __name__ == "__main__":
    main()
