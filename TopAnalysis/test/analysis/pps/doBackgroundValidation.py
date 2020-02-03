import ROOT
import os
import sys
import argparse
import pickle
import copy
from generateBinnedWorkspace import defineProcessTemplates,smoothMissingMass,VALIDLHCXANGLES
from TopLJets2015.TopAnalysis.Plot import *
import numpy as np

def rebinUnequalBinSize(h,newBins):
    nBin = len(newBins)
    newname = h.GetName()+'_new'
    hnew = h.Rebin(nBin-1, newname, np.array(newBins))
    return hnew


def fillShapes(inputDir,selCuts,proc='MuonEG'):

    """fills the shapes from a set of files"""

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(inputDir,x) for x in os.listdir(inputDir) if proc in x]:
        data.AddFile(f)
        
    histos={}
    sf=None
    for dist,hist,title in [('mmiss',                    '(50,0,2000)',             'Missing mass [GeV]'),
                            ('(bosony>=0?mmiss:-mmiss)', '(100,-2000,2000)',        'Missing mass x sgn[y(e#mu)] [GeV]'),
                            ('csi1',                     '(25,0.025,0.18)',         '#xi(+)'),
                            ('csi2',                     '(25,0.025,0.18)',         '#xi(-)'),
                            ('mpp',                      '(25,500,2500)',           'Di-proton mass [GeV]'),
                            ('ypp',                      '(50,-1,1)',               'Di-proton rapidity'),
                            ('mmiss:mpp',                '(25,500,2500,30,0,2500)', 'Missing mass [GeV];Di-proton mass [GeV]'),
                            ('ypp:mpp',                  '(25,500,2500,50,-1,1)',   'Di-proton rapidity;Di-proton mass [GeV]'),                        
                            ]:

        histos[dist]={}

        bkgHistos=[]
        for tag,pf,mixType in [ ('data',               '',        0),   #observed data
                                ('bkg',                '',        1),   #mixed data
                                ('bkgshape',           'syst',    1),   #mixed data from e-mu
                                ('bkgsinglediffUp',    '',        2),   #mix only pos. arm
                                ('bkgsinglediffDown',  'syst',    2),   #mix only neg. arm
                                ]:
            finaldist=dist
            for c in ['mmiss','ypp','mpp','csi1','csi2']:
                finaldist=finaldist.replace(c,pf+c)
                
            #if dist=='mmiss' or dist=='(bosony>=0?mmiss:-mmiss)':
            #    from array import array
            #    newBins=range(0,2000,40)
            #    #newBins=[0,100]+range(200,1000,40)+[1000,1100,1200,1500,2000]
            #    if dist!='mmiss':
            #        newBins=range(-2000,2000,40)
            #        #nbins=len(newBins)
            #        #newBins=[-newBins[ix] for ix in range(nbins-1,0,-1)] + newBins
            #    h=ROOT.TH1F( 'hmiss',';Missing mass [GeV];Events',len(newBins)-1,array('d',newBins) )
            #    data.Draw('{0} >> hmiss'.format(finaldist),
            #              'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(selCuts,pf,mixType),
            #              'goff')
            #else:

            finalSelCuts=selCuts
            if pf=='syst' : 
                finalSelCuts=finalSelCuts.replace('protonCat','systprotonCat')
            print finalSelCuts
            data.Draw('{0} >> h{1}'.format(finaldist,hist),
                      'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(selCuts,pf,mixType),
                      'goff')
                
            h=data.GetHistogram()
            histos[dist][tag]=h.Clone('{0}_{1}_{2}_obs'.format(dist,tag,proc))
            if ';' in title:
                xtit,ytit=title.split(';')
                histos[dist][tag].GetZaxis().SetTitle('Events')
                histos[dist][tag].GetYaxis().SetTitle(ytit)
                histos[dist][tag].GetXaxis().SetTitle(xtit)
            else:
                histos[dist][tag].GetYaxis().SetTitle('Events')
                histos[dist][tag].GetXaxis().SetTitle(title)
            histos[dist][tag].SetDirectory(0)
            h.Reset('ICE')

            if tag=='data' : 
                continue

            #apply a smoothing procedure for 1D mmiss
            #if 'mmiss' in dist and not h.InheritsFrom('TH2'):
            #    rawH=histos[dist][tag]
            #    histos[dist][tag]=smoothMissingMass(rawH)
            #    histos[dist][tag].SetName(rawH.GetName())
            #    histos[dist][tag].SetDirectory(0)
            #    rawH.Delete()
            
            bkgHistos.append(histos[dist][tag])

        h.Delete()

        #use low missing mass to fix possible normalization differences
        if sf is None:
            sf=histos[dist]['data'].Integral()/histos[dist]['bkg'].Integral()
            #upBin=histos[dist]['data'].GetXaxis().FindBin(500.)
            #sf=histos[dist]['data'].Integral(0,upBin)/histos[dist]['bkg'].Integral(0,upBin)
            print 'Scale factor:',sf

        #scale background estimates
        for tag in histos[dist]:
            if not 'bkg' in tag: continue
            histos[dist][tag].Scale(sf)

        #mirror shapes
        finalTemplates=defineProcessTemplates(bkgHistos)
        for h in finalTemplates:
            histos[dist][h.GetName()]=h

    return histos



def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--sig',
                        dest='sig',
                        default='{boson}_m_X_{mass}_xangle_{xangle}_2017_preTS2_opt_v1_simu_reco.root',
                        help='signal point [%default]')
    parser.add_argument('--massList',
                        dest='massList',
                        default='800,1000,1200',
                        help='signal mass list (CSV) [%default]')
    parser.add_argument('--selCuts',
                        dest='selCuts', 
                        default='bosonpt>40 && l1pt>30 && l2pt>20',
                        help='preselection for Z categories [default: %default]')
    parser.add_argument('--doPerAngle',
                        dest='doPerAngle', 
                        default=False,
                        help='do per crossing angle [default: %default]',
                        action='store_true')
    parser.add_argument('--protonCat',
                        dest='protonCat', 
                        default=None,
                        type=int,
                        help='proton reco category [default: %default]')
    parser.add_argument('--lumi',
                        dest='lumi',
                        default=37500.,
                        type=float,
                        help='integrated luminosity [default: %default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='analysis/bkg',
                        help='Output directory [default: %default]')
    opt=parser.parse_args(args)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    os.system('mkdir -p %s'%opt.output)

    #upate pre-selection
    catTitle='inclusive'
    if opt.protonCat:
        opt.selCuts += ' && protonCat==%d'%opt.protonCat
        if opt.protonCat==1:
            catTitle='multi-multi'
        if opt.protonCat==2:
            catTitle='multi-single'
        if opt.protonCat==3:
            catTitle='single-multi'
        if opt.protonCat==4:
            catTitle='single-single'

    selCuts=[ (opt.selCuts,'','e#mu, %s'%catTitle) ]
    if opt.doPerAngle:
        for xangle in VALIDLHCXANGLES:
            selCuts.append( (opt.selCuts + ' && xangle==%d'%xangle,'_%d'%xangle,'e#mu, %s, %d #murad'%(catTitle,xangle)) )

    for cuts,pfix,catTitle in selCuts:

        if opt.protonCat:
            pfix += 'pp%d'%opt.protonCat

        os.system('rm %s/plotter_%s.root'%(opt.output,pfix))

        data=fillShapes(inputDir=opt.input,selCuts=cuts,proc='MuonEG')

        #sigs={}
        #for m in opt.massList.split():
        #    if not m in sigs: sigs[m]={}
        #    for xangle in VALIDLHCXANGLES:
        #        newSigs=fillShapes(inputDir=opt.input,selCuts=opt.selCuts,tag=opt.sig.format(mass=m,xangle=xangle))
        #        for dist in newSigs:
        #            if not dist in sigs[m]:
        #                sigs[m]=newSigs[dist]['data'].Clone('{0}_{1}'.format(dist,m))
        #            #FIXME scale me according to xsec
        #FIXME plot me

        for dist in data:

            pname=dist+pfix
            for c in [':',',','>','=','(',')','-','<','?']: pname=pname.replace(c,'')
            if opt.protonCat:
                pname += 'pp%d'%opt.protonCat

            p=Plot(pname)
            p.doChi2=False #True
            p.nominalDistForSystsName='background'

            #if dist=='mmiss':
            #    newBins=range(0,1000,40)+[1000,1100,1200,1500,2000]
            #    print 'Rebinning',dist,'to',newBins
            #    for k in data[dist]:
            #        data[dist][k]=rebinUnequalBinSize(data[dist][k],newBins)

            #main distributions
            #doDivideByBinWidth=True if dist=='mmiss' or dist=='(bosony>=0?mmiss:-mmiss)' else False
            #if doDivideByBinWidth:
            #    p.doPoissonErrorBars=False

            p.add(data[dist]['data'], title='data',       color=ROOT.kBlack, isData=True, spImpose=False, isSyst=False)#,doDivideByBinWidth=doDivideByBinWidth)
            p.add(data[dist]['bkg'],  title='background', color=ROOT.TColor.GetColor('#1f78b4'), isData=False, spImpose=False, isSyst=False) #,doDivideByBinWidth=doDivideByBinWidth)

            #background systematics
            for syst in ['{0}_bkgshape_MuonEG_obsUp',        
                         '{0}_bkgshape_MuonEG_obsDown',
                         '{0}_bkgsinglediffUp_MuonEG_obs', 
                         '{0}_bkgsinglediffDown_MuonEG_obs', 
                         #'{0}_bkgsinglediffpos_MuonEG_obsUp',
                         #'{0}_bkgsinglediffpos_MuonEG_obsDown',
                         #'{0}_bkgsinglediffneg_MuonEG_obsDown', 
                         #'{0}_bkgsinglediffneg_MuonEG_obsUp'
                         ]:
                p.add(data[dist][syst.format(dist)],  
                      title=syst, 
                      color=ROOT.TColor.GetColor('#1f78b4'), 
                      isData=False, 
                      spImpose=False, 
                      isSyst=True)
                      #doDivideByBinWidth=doDivideByBinWidth)       

            #p.ratiorange=[0.78,1.22]
            p.ratiorange=[0.58,1.42]
            #p.ratiorange=[0.,2.]
            p.show(opt.output,opt.lumi,extraText=catTitle)
            p.appendTo('%s/plotter_%s.root'%(opt.output,pfix))
            p.reset()

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



