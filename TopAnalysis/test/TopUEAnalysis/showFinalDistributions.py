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
    ('PW+PY8', [ ('nominal',         ['t#bar{t}']), 
                 ('#deltaCUET8P2MT4',['t#bar{t} UEup',     't#bar{t} UEdn']),
                 ('FSR',             ['t#bar{t} fsr up',   't#bar{t} fsr dn']),
                 ('ISR',             ['t#bar{t} isr up',   't#bar{t} isr dn']),
                 ('hdamp',           ['t#bar{t} hdamp up', 't#bar{t} hdamp dn']),
                 ('CR',              ['t#bar{t} QCDbased', 't#bar{t} ERDon', 't#bar{t} gluon move']) ] 
     ),
    ('UE up',  [ ('nominal',          ['t#bar{t} UEup'])]),
    ('UE dn',  [ ('nominal',          ['t#bar{t} UEdn'])]),
    ('FSR up', [ ('nominal',         ['t#bar{t} fsr up'])]),
    ('FSR dn', [ ('nominal',         ['t#bar{t} fsr dn'])]),
    ('ISR up', [ ('nominal',         ['t#bar{t} isr up'])]),
    ('ISR dn', [ ('nominal',         ['t#bar{t} isr dn'])]),
    ('QCD based', [ ('nominal',         ['t#bar{t} QCDbased'])]),
    ('ERD on',    [ ('nominal',         ['t#bar{t} ERDon'])]),
    ('Gluon move',[ ('nominal',         ['t#bar{t} gluon move'])]),
    ('aMC@NLO+PY8', [ ('nominal', ['t#bar{t} aMC@NLO']) ]),
    ('PW+HW++'    , [ ('nominal', ['t#bar{t} Herwig++']) ]),
    ]

EXTRASETS = [
    ('Sherpa', 'MC13TeV_TTJets_sherpa.root'),
    ('PW+HW7', 'MC13TeV_TTJets_herwig7.root'),
    ('PW+HW7 dipole', 'MC13TeV_TT_herwig7dipole_meNo.root'),
    ('PW+HW7dipMEon', 'MC13TeV_TT_herwig7dipole_meYes.root'),
    ('no MPI', 'MC13TeV_TTJets_MPIoff.root'),
    ('no CR',  'MC13TeV_TTJets_CRoff.root'),
    ('#alpha_{S}^{FSR}=0.070','MC13TeV_TTJets_pythia8_asfsr0.070_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.080','MC13TeV_TTJets_pythia8_asfsr0.080_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.090','MC13TeV_TTJets_pythia8_asfsr0.090_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.100','MC13TeV_TTJets_pythia8_asfsr0.100_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.105','MC13TeV_TTJets_pythia8_asfsr0.105_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.110','MC13TeV_TTJets_pythia8_asfsr0.110_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.115','MC13TeV_TTJets_pythia8_asfsr0.115_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.120','MC13TeV_TTJets_pythia8_asfsr0.120_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.125','MC13TeV_TTJets_pythia8_asfsr0.125_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.130','MC13TeV_TTJets_pythia8_asfsr0.130_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.135','MC13TeV_TTJets_pythia8_asfsr0.135_meon_crdefault.root'),
    #('#alpha_{S}^{FSR}=0.1365','MC13TeV_TTJets_pythia8_asfsr0.1365_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.140','MC13TeV_TTJets_pythia8_asfsr0.140_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.150','MC13TeV_TTJets_pythia8_asfsr0.150_meon_crdefault.root'),
    ('#alpha_{S}^{FSR}=0.160','MC13TeV_TTJets_pythia8_asfsr0.160_meon_crdefault.root'),                              
    ('#alpha_{S}^{ISR}=0.09', 'MC13TeV_TTJets_pythia8_asisr0.09.root'),
    ('#alpha_{S}^{ISR}=0.10', 'MC13TeV_TTJets_pythia8_asisr0.10.root'),
    ('#alpha_{S}^{ISR}=0.11', 'MC13TeV_TTJets_pythia8_asisr0.11.root'),
    ('#alpha_{S}^{ISR}=0.12', 'MC13TeV_TTJets_pythia8_asisr0.12.root'),
    ('#alpha_{S}^{ISR}=0.13', 'MC13TeV_TTJets_pythia8_asisr0.13.root'),
    ('#alpha_{S}^{ISR}=0.14', 'MC13TeV_TTJets_pythia8_asisr0.14.root'),
    ('#alpha_{S}^{ISR}=0.15', 'MC13TeV_TTJets_pythia8_asisr0.15.root'),
    ('#alpha_{S}^{ISR}=0.16', 'MC13TeV_TTJets_pythia8_asisr0.16.root'),
    #('Rope had.',             'MC13TeV_TTJets_pythia8_asfsr0.120_meon_crdefault_flavrope.root'),
    ('Rope',             'MC13TeV_TTJets_pythia8_asfsr0.1365_meon_crdefault_flavrope.root'),
    ('Rope (no CR)',     'MC13TeV_TTJets_pythia8_asfsr0.1365_meon_croff_flavrope.root'),
]

PLOTTINGSET_1=[
    ('Data',              '2',  1001,  '#a6cee3', 1 , True,  None),
    ('PW+PY8',            'ep0', 0,    '#49494c', 24, False, 0.2),
    ('ISR up',            'ep0', 0,    '#fdc086', 22, False, 0.4),
    ('ISR dn',            'ep0', 0,    '#fdc086', 23, False, 0.4),
    ('FSR up',            'ep0', 0,    '#d95f02', 22, False, 0.6),
    ('FSR dn',            'ep0', 0,    '#d95f02', 23, False, 0.6),
    ('UE up',             'ep0', 0,    '#000000', 22, False, 0.80),
    ('UE dn',             'ep0', 0,    '#000000', 23, False, 0.80),
    #('CP5',               'ep0', 0,    '#91b58d', 24, False, 0.9),
]

PLOTTINGSET_2=[
    ('Data',              '2',   1001, '#a6cee3', 1 , True,  None),
    ('no MPI',            'l',   0,    '#91b58d', 1,  False, None),
    ('no CR',             'l',   0,    '#000000', 1,  False, None),
    ('QCD based',         'ep0', 0,    '#fdc086', 20, False, 0.2),
    ('Gluon move',        'ep0', 0,    '#984ea3', 21, False, 0.4),
    ('ERD on',            'ep0', 0,    '#d95f02', 24, False, 0.6),
    ('Rope',              'ep0', 0,    '#000000', 26, False, 0.8),
    ('Rope (no CR)',      'ep0', 0,    '#000000', 32, False, 0.8),
]

PLOTTINGSET_3=[
    ('Data',           '2',   1001, '#a6cee3', 1 , True,  None),
    ('Sherpa',         'ep0', 0,    '#49494c', 24, False, 0.2),
    ('aMC@NLO+PY8',    'ep0', 0,    '#e41a1c', 21, False, 0.5),
    ('PW+HW++',        'ep0', 0,    '#386cb0', 22, False, 0.8),
    ('PW+HW7',         'ep0', 0,    '#386cb0', 23, False, 0.8),
]

PLOTTINGSET_4=[
    ('Data',          '2',   1001, '#a6cee3', 1 , True,  None),
    ('PW+HW7',        'ep0', 0,    '#386cb0', 23, False, 0.8),
    ('PW+HW7 dipole', 'ep0', 0,    '#49494c', 24, False, 0.5),
    ('PW+HW7dipMEon', 'ep0', 0,    '#984ea3', 21, False, 0.2),
]

def compareUEPlots(uePlots,outDir,cuts,obs,plottingSetList=[PLOTTINGSET_1],cmsLabel='#bf{CMS}'):
    """This method dumps the formatted plots to the canvas"""

    logX=True if 'chflux' in obs or 'chavg' in obs or 'chrecoil' in obs else False

    #start the main canvas
    c=ROOT.TCanvas('c','c',600,600)
    c.SetTopMargin(0.06)
    c.SetRightMargin(0.03)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.14)
    c.SetLogx(logX)
    c.SetGridx()
    
    #frame
    xarr=uePlots['Data'].trueAxis.GetXbins().GetArray()
    if logX and xarr[0]==0: xarr[0]=xarr[1]*0.16
    frame=ROOT.TH1F('frame','frame',uePlots['Data'].trueAxis.GetNbins(),xarr)
    frame.Draw()
    xtitle=VARTITLES[obs]
    if obs in ['chflux','chfluxz','chavgpt','chavgpt','chrecoil'] : xtitle += ' [GeV]'
    frame.GetXaxis().SetTitle(xtitle)
    ytitle='1/N dN/d(%s)'%VARTITLES[obs]
    if obs in ['chflux','chfluxz','chavgpt','chavgpt','chrecoil','chavgpz'] : ytitle += ' [GeV^{-1}]'
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetYaxis().SetTitleOffset(1.45)
    frame.GetXaxis().SetTitleOffset(1.25)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.045)
    if logX: 
        frame.GetXaxis().SetMoreLogLabels()

    #plot and add to the legend
    leg=ROOT.TLegend(0.6,0.92,0.96,0.92-len(plottingSetList)*0.09) 
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry( uePlots['Data'].plot[0], 'Data','f' )

    maxY,minY=-10,10
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSetList[0]:
        try:
            uePlots[p].format(fill,color,marker,keepXUnc,shiftX)
            uePlots[p].plot[0].Draw(drawOpt)
            if p!='Data': 
                legOpt=drawOpt if not drawOpt in ['c','2'] else 'l'
                leg.AddEntry(uePlots[p].plot[0],p,legOpt)
            else:
                uePlots[p].plot[1][0].SetFillStyle(3004)
                uePlots[p].plot[1][0].SetFillColor(1)
                uePlots[p].plot[1][0].Draw(drawOpt)

            iminY,imaxY=getGraphExtremes(uePlots[p].plot[0])
            maxY=max(imaxY,maxY)
            minY=min(iminY,minY)
        except:
            pass
    #frame.GetYaxis().SetRangeUser(minY*0.8,maxY*1.25)
    frame.GetYaxis().SetRangeUser(0,maxY*1.25)
    if minY>0 and maxY/minY<100: frame.GetYaxis().SetRangeUser(0.,maxY*1.25)
    leg.Draw()

    #stat component for data
    txt=ROOT.TPaveText(0.615,0.905,0.675,0.89,'brNDC')
    txt.SetFillStyle(3004)
    txt.SetFillColor(1)
    txt.SetBorderSize(0)
    txt.SetLineColor(0)
    txt.Draw()

    #standard label
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.DrawLatex(0.18,0.87,cmsLabel)
    tex.DrawLatex(0.69,0.96,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')
    icut=0
    for cutKey in cuts:   
        y=0.86-len(plottingSetList[0])*0.06-0.06*icut
        if cutKey=='region': 
            regionName='toward'
            if cuts[cutKey][1]==1: regionName='transverse'
            if cuts[cutKey][1]==2: regionName='away'
            tex.DrawLatex(0.62,y,'#scale[0.8]{%s [%s]}'%(regionName,VARTITLES[cuts[cutKey][0]]) )
        else :
            tex.DrawLatex(0.62,y,'#scale[0.8]{%3.1f#leq%s<%3.1f}'%(cuts[cutKey][0],VARTITLES[cutKey],cuts[cutKey][1]))
        icut+=1
    c.RedrawAxis()
    c.Modified()
    c.Update()
    for ext in ['pdf','png','root']: c.SaveAs('%s/%s_unfolded.%s'%(outDir,obs,ext))

    #compute the ratios to data
    uePlotRatios=getRatiosWithRespectTo(uePlots,'Data')
    uePlotRatiosStat=getRatiosWithRespectTo(uePlots,'Data','stat')

    #start the ratio canvas
    dy_xtit=10.
    dy_pad=180.
    y_total=2*dy_xtit+dy_pad*len(plottingSetList)
    rel_dy_xtit=dy_xtit/y_total
    rel_dy_pad=dy_pad/y_total
    cratio=ROOT.TCanvas('c','c',800,int(y_total))
    cratio.SetTopMargin(0.0)
    cratio.SetRightMargin(0.0)
    cratio.SetLeftMargin(0.0)
    cratio.SetBottomMargin(0.0)
    sp=[]
    rf=[]
    lg=[]
    y=1
    for i in xrange(0,len(plottingSetList)):

        cratio.cd()

        #start new sub-pad
        dy_ipad=rel_dy_pad
        if i==0 or i==len(plottingSetList)-1: dy_ipad+=rel_dy_xtit
        sp.append( ROOT.TPad('p%d'%i,'p%d'%i,0,y,1.0,max(y-dy_ipad,0.)) )
        y=max(y-dy_ipad,0.)

        #margins
        if i==0:
            sp[-1].SetTopMargin(0.23)
        else:
            sp[-1].SetTopMargin(0.11)
        sp[-1].SetRightMargin(0.25)
        sp[-1].SetLeftMargin(0.12)
        if i==len(plottingSetList)-1:
            sp[-1].SetTopMargin(0.02)
            sp[-1].SetBottomMargin(0.3) #2*dy_xtit/dy_ipad)
        elif i==0:
            sp[-1].SetBottomMargin(0.11)
        else:
            sp[-1].SetBottomMargin(0.2)
        sp[-1].SetGridx()
        sp[-1].SetLogx(logX)
        sp[-1].Draw()
        sp[-1].cd()        

        rf.append( frame.Clone('ratioframe') )
        rf[-1].GetYaxis().SetTitle('Theory/Data')
        rf[-1].GetXaxis().SetTitle(xtitle)
        rf[-1].GetYaxis().SetTitleOffset(0.3)
        rf[-1].GetYaxis().SetTitleSize(0.13)
        rf[-1].GetYaxis().SetLabelSize(0.11)
        if i==len(plottingSetList)-1:
            rf[-1].GetXaxis().SetTitleSize(0.13)
            rf[-1].GetXaxis().SetLabelSize(0.11)
            rf[-1].GetXaxis().SetTitleOffset(1.0)
            rf[-1].GetXaxis().SetTitle(xtitle)
        else:
            rf[-1].GetXaxis().SetTitleSize(0.0)
            rf[-1].GetXaxis().SetLabelSize(0.0)
        rf[-1].Draw()

        if i==0:
            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.13)
            tex.SetNDC()
            tex.DrawLatex(0.12,0.83,cmsLabel)
            tex.DrawLatex(0.56,0.83,'#scale[0.9]{35.9 fb^{-1} (13 TeV)}')
        
        if i==0:
            lg.append( ROOT.TLegend(0.75,0.1,0.99,0.82) )
        elif i==len(plottingSetList)-1:
            lg.append( ROOT.TLegend(0.75,0.23,0.99,0.95) )
        else:
            lg.append( ROOT.TLegend(0.75,0.1,0.99,0.95) )
        lg[-1].SetFillStyle(0)
        lg[-1].SetBorderSize(0)
        lg[-1].SetTextFont(42)
        lg[-1].SetTextSize(0.1)
        #lg[-1].SetNColumns(len(plottingSetList[i])-1)
    
        maxY,minY=1.55,0.45
        for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSetList[i]:
            try:
                formatGraph(uePlotRatios[p],fill,color,marker,keepXUnc,shiftX)
                uePlotRatios[p].Draw(drawOpt)
                iminY,imaxY=getGraphExtremes(uePlotRatios[p])
                maxY=max(imaxY,maxY)
                minY=min(iminY,minY)
                if p!='Data': 
                    legOpt=drawOpt if not drawOpt in ['c','2'] else 'l'
                    lg[-1].AddEntry(uePlotRatios[p],p,legOpt)
                else:
                    uePlotRatiosStat[p].SetFillStyle(3004)
                    uePlotRatiosStat[p].SetFillColor(1)
                    uePlotRatiosStat[p].Draw(drawOpt)
            except:
                print 'Skip',p
                pass
        rf[-1].GetYaxis().SetRangeUser(max(minY*0.8,0.15),min(maxY*1.25,1.95))
        rf[-1].GetYaxis().SetNdivisions(5)
        if logX: rf[-1].GetXaxis().SetMoreLogLabels()
        
        sp[-1].RedrawAxis()
        lg[-1].Draw()
    
    # all done
    cratio.Modified()
    cratio.Update()
    for ext in ['pdf','png','root']: cratio.SaveAs('%s/%s_unfolded_ratio.%s'%(outDir,obs,ext))

"""
"""
def showSystsSummary(systsH,outdir,cuts,obs,cmsLabel):

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
    frame.GetYaxis().SetRangeUser(0.,0.3)
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
    tex.DrawLatex(0.15,0.9,cmsLabel)
    tex.DrawLatex(0.69,0.96,'#scale[0.7]{35.9 fb^{-1} (13 TeV)}')
    tex.DrawLatex(0.72,0.9,'#scale[0.7]{%s}'%VARTITLES[obs])
    icut=0
    for cutKey in cuts:
        y=0.7-0.04*icut
        cutText='inclusive'
        if cutKey=='region': cutText='%s region=%s'%(VARTITLES[cuts[cutKey][0]],cuts[cutKey][1])
        else               : cutText='%3.1f#leq%s<%3.1f'%(cuts[cutKey][0],VARTITLES[cutKey],cuts[cutKey][1])
        tex.DrawLatex(0.75,y,'#scale[0.7]{%s}'%cutText)
        icut+=1

    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/%s_reluncertainty.%s'%(outdir,obs,ext))


def readParticlePlotsFrom(baseAnaDir,args,obsAxis,cuts,obs):
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
        ('Trk. eff.'      , 'exp') : ['tkeffeta','tkeffdstar'],
        ('UE'             , 'th')  : ['t#bar{t} UEdn',     't#bar{t} UEup'],
        ('FSR'            , 'th')  : ['t#bar{t} fsr dn',   't#bar{t} fsr up'],
        ('ME-PS'          , 'th')  : ['t#bar{t} hdamp dn', 't#bar{t} hdamp up'],
        ('ISR'            , 'th')  : ['t#bar{t} isr dn',   't#bar{t} isr up'],
        ('m_{t}'          , 'th')  : ['t#bar{t} m=171.5',  't#bar{t} m=173.5'],
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

    hdata=fIn.Get('corrected_data')
    uePlots['Data'].addVariation('Data',None, hdata)
    for key in finalSystList:
        for v in finalSystList[key]:
            uePlots['Data'].addVariation(key[0],key[1],fIn.Get('corrected_data%s'%v))

    #statistical covariance matrix
    EmatTotal=fIn.Get("EmatTotal")
    nx,ny=EmatTotal.GetNbinsX(),EmatTotal.GetNbinsY()
    uePlots['Data'].covMatrices['stat']=ROOT.TMatrixF(nx,ny) 
    norm=hdata.Integral()
    for xbin in xrange(1,nx+1):
        if norm<=0 : continue
        for ybin in xrange(1,ny+1):
            uePlots['Data'].covMatrices['stat'][xbin-1][ybin-1]=EmatTotal.GetBinContent(xbin,ybin)/(norm**2)            

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

    #add extra sets (generator level only)
    for varTitle,varUrl in EXTRASETS:
        url=os.path.join( baseAnaDir, varUrl)
        if not os.path.isfile(url): 
            print 'skip',url
            continue
        fIn=ROOT.TFile.Open(url)
        genH=fIn.Get('gen')
        varName='gen_%d'%len(uePlots)
        uePlots[varTitle]=UEPlot(varName,varTitle,obsAxis)
        uePlots[varTitle].addVariation(varTitle,None,genH)
        fIn.Close()
    

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
    parser.add_option('--cmsLabel',
                      dest='cmsLabel',
                      help='cmsLabel [%default]',
                      default='#bf{CMS} #it{preliminary}',
                      type='string')
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
    baseAnaDir=os.path.dirname(opt.analysisCfg)

    uePlots={}
    if opt.reco:
        readRecoPlotsFrom(args,opt)
    else:
        uePlots=readParticlePlotsFrom(baseAnaDir,args,obsAxis,cuts,obs)

    #finalize the plots
    failed=[]
    for key in uePlots: 
        try:
            uePlots[key].finalize(doCov=True)
        except:
            print 'Failed to finalize plot for', key
            failed.append(key)

    #remove if failed
    for key in failed: uePlots.pop(key)


    #show plots
    outDir=os.path.dirname(args[0])
    compareUEPlots(uePlots=uePlots,
                   outDir=outDir,
                   cuts=cuts,
                   obs=obs,
                   plottingSetList=[PLOTTINGSET_1,PLOTTINGSET_2,PLOTTINGSET_3], #PLOTTINGSET_4],
                   cmsLabel=opt.cmsLabel)
    
    showSystsSummary(uePlots['Data'].relUncertaintyH,
                     outdir=outDir,
                     cuts=cuts,
                     obs=obs,
                     cmsLabel=opt.cmsLabel)

    with open(os.path.join(outDir,'mean_summary.dat'),'w') as cachefile:
        cachefile.write(uePlots['Data'].meanUncTable)

    with open(os.path.join(outDir,'unfold_summary.pck'),'w') as cachefile:
        pickle.dump( uePlots, cachefile, pickle.HIGHEST_PROTOCOL)

    

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
