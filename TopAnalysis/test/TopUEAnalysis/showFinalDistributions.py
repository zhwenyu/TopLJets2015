import ROOT
import os
import sys
import optparse
import numpy as np
import pickle
from collections import OrderedDict
from UEAnalysisHandler import VARTITLES,SYSTS
from UEPlot import *
from UETools import getGraphExtremes

COMPARISONSETS=[
    ('PW+PY8 CUETP8M2T4', [ ('nominal',         ['t#bar{t}']), 
                            ('#deltaCUET8P2MT4',['t#bar{t} UEup',     't#bar{t} UEdn']),
                            ('FSR',             ['t#bar{t} fsr up',   't#bar{t} fsr dn']),
                            ('ISR',             ['t#bar{t} isr up',   't#bar{t} isr dn']),
                            ('hdamp',           ['t#bar{t} hdamp up', 't#bar{t} hdamp dn']),
                            ('CR',              ['t#bar{t} QCDbased', 't#bar{t} ERDon', 't#bar{t} gluon move']) ] 
     ),
    ('FSR up', [ ('nominal',         ['t#bar{t} fsr up'])]),
    ('FSR dn', [ ('nominal',         ['t#bar{t} fsr dn'])]),
    ('ISR up', [ ('nominal',         ['t#bar{t} isr up'])]),
    ('ISR dn', [ ('nominal',         ['t#bar{t} isr dn'])]),
    ('QCD based CR', [ ('nominal',         ['t#bar{t} QCDbased'])]),
    ('ERD on CR',    [ ('nominal',         ['t#bar{t} ERDon'])]),
    ('Gluon move CR',[ ('nominal',         ['t#bar{t} gluon move'])]),
    ('aMC@NLO+PY8 CUETP8M2T4', [ ('nominal', ['t#bar{t} aMC@NLO']) ]),
    ('PW+HW++ EE5C'      , [ ('nominal', ['t#bar{t} Herwig++']) ]),
    ]

PLOTTINGSETS=[
    ('Data',              '2',  1001, '#26e278', 1 , True,  None),
    ('PW+PY8 CUETP8M2T4', 'ep0', 0,    '#000000', 20, False, 0.2),
    ('QCD based CR',      'ep0', 0,    '#ef8a62', 24, False, 0.5),
    ('ERD on CR',         'ep0', 0,    '#ef8a62', 25, False, 0.5),
    ('Gluon move CR',     'ep0', 0,    '#b2182b', 27, False, 0.5),
    ('FSR up',            'ep0', 0,    '#67a9cf', 22, False, 0.8),
    ('FSR dn',            'ep0', 0,    '#67a9cf', 23, False, 0.8)
]

def compareUEPlots(uePlots,outDir,cuts,obs):
    """This method dumps the formatted plots to the canvas"""

    #start the canvas
    c=ROOT.TCanvas('c','c',600,600)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetBottomMargin(0.0)

    c.cd()
    p1=ROOT.TPad('p1','p1',0,0.3,1.0,1.0)
    p1.SetTopMargin(0.06)
    p1.SetRightMargin(0.03)
    p1.SetLeftMargin(0.12)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.cd()
    frame=ROOT.TH1F('frame','frame',uePlots['Data'].trueAxis.GetNbins(),uePlots['Data'].trueAxis.GetXbins().GetArray())
    frame.Draw()
    frame.GetYaxis().SetTitle('1/N dN/d%s'%VARTITLES[obs])
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.04)


    #plot and add to the legend
    leg=ROOT.TLegend(0.6,0.92,0.96,0.92-len(PLOTTINGSETS)*0.06) 
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry( uePlots['Data'].plot[0], 'Data','f' )
    maxY,minY=-10,10
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in PLOTTINGSETS:
        try:
            uePlots[p].format(fill,color,marker,keepXUnc,shiftX)
            uePlots[p].plot[0].Draw(drawOpt)
            if p!='Data': leg.AddEntry(uePlots[p].plot[0],p,drawOpt)
            iminY,imaxY=getGraphExtremes(uePlots[p].plot[0])
            maxY=max(imaxY,maxY)
            minY=min(iminY,minY)
        except:
            pass
    frame.GetYaxis().SetRangeUser(minY*0.8,maxY*1.25)
    leg.Draw()

    #standard label
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.87,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.73,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
    icut=0
    for cutKey in cuts:
        y=0.89-0.04*icut
        if cutKey=='region': tex.DrawLatex(0.55,y,'%s region=%s'%(VARTITLES[cuts[cutKey][0]],cuts[cutKey][1]) )
        else               : tex.DrawLatex(0.55,y,'%3.1f#leq%s<%3.1f'%(cuts[cutKey][0],VARTITLES[cutKey],cuts[cutKey][1]))
        icut+=1

    p1.RedrawAxis()

    c.cd()
    p2=ROOT.TPad('p2','p2',0,0.0,1.0,0.3)
    p2.SetTopMargin(0.01)
    p2.SetRightMargin(0.03)
    p2.SetLeftMargin(0.12)
    p2.SetBottomMargin(0.3)
    p2.Draw()
    p2.cd()
    ratioframe=frame.Clone('ratioframe')
    ratioframe.GetYaxis().SetTitle('Ratio to Data')
    ratioframe.GetXaxis().SetTitle(VARTITLES[obs])
    ratioframe.GetYaxis().SetTitleOffset(0.5)
    ratioframe.GetYaxis().SetTitleSize(0.11)
    ratioframe.GetXaxis().SetTitleSize(0.11)
    ratioframe.GetYaxis().SetLabelSize(0.1)
    ratioframe.GetXaxis().SetLabelSize(0.1)

    ratioframe.Draw()
    uePlotRatios=getRatiosWithRespectTo(uePlots,'Data')
    maxY,minY=1.55,0.45
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in PLOTTINGSETS:
        try:
            uePlotRatios[p].Draw(drawOpt)
            iminY,imaxY=getGraphExtremes(uePlotRatios[p])
            maxY=max(imaxY,maxY)
            minY=min(iminY,minY)
        except:
            pass
    ratioframe.GetYaxis().SetRangeUser(max(minY*0.8,0.15),min(maxY*1.25,1.95))
    ratioframe.GetYaxis().SetNdivisions(5)

    p2.RedrawAxis()
    
    # all done
    c.Modified()
    c.Update()
    for ext in ['pdf','png']: c.SaveAs('%s/%s_unfolded.%s'%(outDir,obs,ext))


"""
"""
def showSystsSummary(systsH,outdir,cuts,obs):

    #show systematics
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.03)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    leg=ROOT.TLegend(0.15,0.88,0.65,0.62)
    leg.SetNColumns(3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    nhistos=len(systsH)
    frame=systsH['Data'].Clone('frame')
    frame.Reset('ICE')
    frame.Draw()
    frame.GetYaxis().SetRangeUser(0.,0.5)
    frame.GetYaxis().SetTitle('Relative uncertainty')
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetTitle('Bin number')
    systsH['Data'].Draw('histsame')
    systsH['Data'].SetLineWidth(2)
    leg.AddEntry(systsH['Data'],'Stats','l')
    icontrib=0
    for key in systsH:
        if key=='Data': continue
        systsH[key].Draw('histsame')
        systsH[key].SetLineColor(icontrib%3+1)
        systsH[key].SetLineStyle(icontrib/3+1)
        leg.AddEntry(systsH[key],key,'l')
        icontrib+=1
    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.67,0.96,'#scale[0.7]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/%s_reluncertainty.%s'%(outdir,obs,ext))

    #save summary to file
    #with open('%s/%s_reluncertainty.dat'%(outdir,outname),'w') as fout:
    #    fout.write('Source & <Rel. unc.> \\\\\n')
    #    for s in systSummary:
    #        fout.write('%20s & %3.1f \\\\\n'%(s,systSummary[s]))


def readParticlePlotsFrom(args,obsAxis,cuts,obs):
    """
    reads the results contained in the unfold summary
    """

    uePlots={}
        
    #list of systematics
    systList={
        ('Pileup'         , 'exp') : ['puup','pudn'],
        ('Trigger/Sel.'   , 'exp') : ['effup','effdn'],
        ('p_{T}(top)  '   , 'th')  : ['toppt'],
        ('#mu_{R}/#mu_{F}', 'th')  : ['murup','murdn','mufup','mufdn','qup','qdn'],
        ('b-tag'          , 'exp') : ['btagup','btagdn'],
        ('JES'            , 'exp') : ['jesup','jesdn'],
        ('JER'            , 'exp') : ['jerup','jerdn'],
        ('LES'            , 'exp') : ['eesup','eesdn','mesup','mesdn'],
        ('Trk. eff.'      , 'exp') : ['tkeff','tkeffeta','tkeffdstar'],
        ('UE'             , 'th')  : ['t#bar{t} UEdn',     't#bar{t} UEup'],
        ('FSR'            , 'th')  : ['t#bar{t} fsr dn',   't#bar{t} fsr up'],
        ('ME-PS'          , 'th')  : ['t#bar{t} hdamp dn', 't#bar{t} hdamp up'],
        ('ISR'            , 'th')  : ['t#bar{t} isr dn',   't#bar{t} isr up'],
        ('m_{t}'          , 'th')  : ['t#bar{t} m=169.5',  't#bar{t} m=175.5'],
        ('Background'     , 'exp') : ['bckpfakes'],
        #('PY8-HW++'      , 'th')  : ['t#bar{t} Herwig++'],
        }

    finalSystList={}
    for key in systList:
        newVars=[]
        for v in systList[key]:
            newv=v
            for i in xrange(0,len(SYSTS)):
                if SYSTS[i][0]!=v : continue
                newv='_%d'%i
            newVars.append(newv)
        finalSystList[key]=newVars    

    #read unfolded data
    fIn=ROOT.TFile(args[0])
    uePlots['Data']=UEPlot(obs,VARTITLES[obs],obsAxis)
    uePlots['Data'].addVariation('Data',None, fIn.Get('corrected_data') )
    for key in finalSystList:
        for v in finalSystList[key]:
            uePlots['Data'].addVariation(key[0],key[1],fIn.Get('corrected_data%s'%v))
    fIn.Close()

    #read sets to compare (use only the nominal one)          
    fGen=ROOT.TFile.Open(args[1])
    fSyst=ROOT.TFile.Open(args[2])    
    for varTitle,subVars in COMPARISONSETS:

        varName='gen_%d'%len(uePlots)
        uePlots[varTitle]=UEPlot(varName,varTitle,obsAxis)
        for x,xvars in subVars:

            for mc in xvars:
                h=None
                key='gen/gen_%s'%mc
                for f in [fGen,fSyst]:
                    try:
                        h=f.Get(key).Clone(varName+mc)
                        break
                    except:                       
                        pass
                
                if h is None: continue
                uePlots[varTitle].addVariation(mc,None if x=='nominal' else 'th',h)    
    fGen.Close()
    fSyst.Close()

    #all done here, return result
    return uePlots


def readRecoPlotsFrom(args,opt):
    """
    FIXME
    """
    outdir=opt.output

    fIn=ROOT.TFile.Open(args[0])
    fSyst=ROOT.TFile.Open(args[1])


    #list of systematics
    systList={
        'Pileup'         : ['puup','pudn'],
        'Trigger/Sel.'   : ['effup','effdn'],
        'b-tag'          : ['btagup','btagdn'],
        'JES'            : ['jesup','jesdn'],
        'JER'            : ['jerup','jerdn'],
        'LES'            : ['eesup','eesdn','mesup','mesdn'],
        'Trk. eff.'      : ['tkeff','tkeffeta']
        }

    #build plots
    for obs in OBSERVABLES:

        obsAxis=analysisaxis[(obs,True)]
        
        for s in SLICES:

            outname=obs
            sliceAxis=None if s is None else analysisaxis[(s,True)]
            if sliceAxis: outname += '%s'%s
            #if not sliceAxis is None: continue
            
            key='%s_%s_inc_None_True'%(obs,s)
        
            #read data, signal and total background
            data,signal,bkg=None,None,None
            t=fIn.Get(key)
            for pkey in t.GetListOfKeys():
                h=t.Get(pkey.GetName())
                if not h.InheritsFrom('TH1') : continue
                if 'Data' in h.GetTitle():
                    data=h.Clone('data')
                elif h.GetTitle() in opt.signal:
                    signal=h.Clone('%s_nominal'%h.GetName())
                else:
                    if bkg is None: bkg=h.Clone('bkg')
                    else : bkg.Add(h)

            data.Add(bkg,-1)
            normalizePerSlice(data,obsAxis,sliceAxis)    

            #experimental systematics
            expSystsKey='%s_%s_inc_syst_True'%(obs,s)
            expSystsH=fIn.Get('{0}/{0}_{1}'.format(expSystsKey,opt.signal))
            addSystematics(signal,expSystsH,systList,obsAxis,sliceAxis,outdir,outname,True)

            #read sets to compare            
            signalVars=[]
            for varTitle,subVars in COMPARISONSETS:

                subVarHistColl=[]
                for x,xvars in subVars:
                    
                    histoColl=[]
                    for ixvar in xvars:
                        if varTitle==MAINMC[0] and ixvar==MAINMC[1]:  
                            histoColl.append(signal.Clone(ixvar))
                        else :
                            histoColl.append( fSyst.Get('{0}/{0}_{1}'.format(key,ixvar)).Clone(ixvar) )
                            normalizePerSlice(histoColl[-1],obsAxis,sliceAxis)
                        histoColl[-1].SetDirectory(0)
                    subVarHistColl.append( (x,histoColl) )
                signalVars.append( (varTitle,subVarHistColl) )

            buildPlot(data,signalVars,obsAxis,sliceAxis,opt,outname,True)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--signal',
                      dest='signal',
                      help='signal [%default]',
                      type='string',
                      default='t#bar{t}')
    parser.add_option('--reco',
                      dest='reco',
                      help='reco flag [%default]',
                      default=False,
                      action='store_true')
    parser.add_option('--cfg',
                      dest='analysisCfg',
                      help='analysis configuration file [%default]',
                      type='string',
                      default='%s/src/TopLJets2015/TopAnalysis/UEanalysis/chmult/inc/analysiscfg.pck'%os.environ['CMSSW_BASE'])
    (opt, args) = parser.parse_args()

    #analysis configuration
    analysiscfg,cuts,obs=None,None,None
    with open(opt.analysisCfg,'r') as cachefile: 
        analysiscfg = pickle.load(cachefile)
        cuts        = pickle.load(cachefile)
        obs         = pickle.load(cachefile)
    obsAxis=analysiscfg[('reco' if opt.reco else 'gen','axis')]


    uePlots={}
    if opt.reco:
        readRecoPlotsFrom(args,opt)
    else:
        uePlots=readParticlePlotsFrom(args,obsAxis,cuts,obs)


    #finalize the plots
    for key in uePlots: uePlots[key].finalize()

    #show plots
    outDir=os.path.dirname(args[0])
    compareUEPlots(uePlots=uePlots,
                   outDir=outDir,
                   cuts=cuts,obs=obs)
    showSystsSummary(uePlots['Data'].relUncertaintyH,
                     outdir=outDir,
                     cuts=cuts,
                     obs=obs)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
