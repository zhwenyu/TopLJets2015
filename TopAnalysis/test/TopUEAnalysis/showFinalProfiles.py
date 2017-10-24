import ROOT
import os
import sys
import optparse
import numpy as np
import pickle
from UEAnalysisHandler import VARTITLES,SYSTS
from UETools import getGraphExtremes,formatGraph
from showFinalDistributions import PLOTTINGSET_1,PLOTTINGSET_2,PLOTTINGSET_3

SLICELISTS=[
    ('',[ 
        ('inc',             'inclusive'),
        ('nj=0,1',          'N_{j}=0'),
        ('nj=1,2',          'N_{j}=1'),
        ('nj=2,999',        'N_{j}#geq2'),
        ('ptll=0,20',       'p_{T}(ll)<20'),
        ('ptll=20,60',      '20<p_{T}(ll)<60'),
        ('ptll=60,120',     '60<p_{T}(ll)<120'),
        ('ptll=120,9999',   'p_{T}(ll)>120'),
        ('mll=0,60',        'M(ll)<60'),
        ('mll=60,120',      '60<M(ll)<120'),
        ('mll=120,200',     '120<M(ll)<200'),
        ('mll=200,9999',    'M(ll)>200'),
        ] ),
    ('_ptll',[
            ('inc_ptll=tow',           'toward'),
            ('ptll=0,20_ptll=tow',     '[0,20['),
            ('ptll=20,60_ptll=tow',    '[20,60['),
            ('ptll=60,120_ptll=tow',   '[60,120['),
            ('ptll=120,9999_ptll=tow', '#geq120'),
            ('inc_ptll=tra',           'transverse'),
            ('ptll=0,20_ptll=tra',     '[0,20['),
            ('ptll=20,60_ptll=tra',    '[20,60['),
            ('ptll=60,120_ptll=tra',   '[60,120['),
            ('ptll=120,9999_ptll=tra', '#geq120'),
            ('inc_ptll=awa',           'away'),
            ('ptll=0,20_ptll=awa',     '[0,20['),
            ('ptll=20,60_ptll=awa',    '[20,60['),
            ('ptll=60,120_ptll=awa',   '[60,120['),
            ('ptll=120,9999_ptll=awa', '#geq120'),
            ]),
    ('_ptllnj',[
            ('inc_ptll=tow',      'toward'),
            ('nj=0,1_ptll=tow',   '=0'),
            ('nj=1,2_ptll=tow',   '=1'),
            ('nj=2,999_ptll=tow', '#geq2'),
            ('inc_ptll=tra',      'transverse'),
            ('nj=0,1_ptll=tra',   '=0'),
            ('nj=1,2_ptll=tra',   '=1'),
            ('nj=2,999_ptll=tra', '#geq2'),
            ('inc_ptll=awa',      'away'),
            ('nj=0,1_ptll=awa',   '=0'),
            ('nj=1,2_ptll=awa',   '=1'),
            ('nj=2,999_ptll=awa', '#geq2'),
            ])
    ]



def getProfileRatiosWithRespectTo(grColl,refKey):
    """
    Compute the ratio using a given key as reference
    """
    grCollRatios,grCollPulls={},{}
    x,y=ROOT.Double(0),ROOT.Double(0)
    xref,yref=ROOT.Double(0),ROOT.Double(0)    
    for key in grColl:
        grCollRatios[key]=grColl[key].Clone('%s_2_%s_ratio'%(grColl[key].GetName(),grColl[refKey].GetName()))
        grCollPulls[key]=grColl[key].Clone('%s_2_%s_pull'%(grColl[key].GetName(),grColl[refKey].GetName()))
        for np in xrange(0,grColl[key].GetN()):

            grColl[key].GetPoint(np,x,y)
            grColl[refKey].GetPoint(np,xref,yref)
            ratio=-1 if yref==0 else y/yref

            ex=grColl[key].GetErrorX(np)
            ey=grColl[key].GetErrorY(np)
            eyref=grColl[refKey].GetErrorY(np)
                
            ratioUnc=(ey*yref)**2
            #if key!=refKey:ratioUnc+=(eyref*y)**2
            ratioUnc=0 if yref==0 else ROOT.TMath.Sqrt(ratioUnc)/(yref**2)
            
            grCollRatios[key].SetPoint(np,x,ratio)
            grCollRatios[key].SetPointError(np,ex,ratioUnc)
            
            pull=(y-yref)/eyref if eyref>0 else 0
            epull=ey/eyref if eyref>0 else 0.
            if key==refKey : epull=1.0
            grCollPulls[key].SetPoint(np,x,pull)
            grCollPulls[key].SetPointError(np,ex,epull)

    return grCollRatios,grCollPulls


def showProfile(grColl,grCollComp,obs,sliceList,plottingSet,outDir,pfix='',isPull=False):

    """
    build the profile plot
    """
    #start the canvas
    c=ROOT.TCanvas('c','c',800,600)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetBottomMargin(0.0)

    c.cd()
    p1=ROOT.TPad('p1','p1',0,0.4,1.0,1.0)
    p1.SetTopMargin(0.06)
    p1.SetRightMargin(0.03)
    p1.SetLeftMargin(0.12)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.SetGridx()
    p1.cd()
    frame=ROOT.TH1F('frame',';Phase space region;<%s>;'%VARTITLES[obs],len(sliceList),0,len(sliceList))
    for xbin in xrange(1,frame.GetNbinsX()+1): frame.GetXaxis().SetBinLabel(xbin,sliceList[xbin-1][1])
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleSize(0.0)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.Draw()

    #plot and add to the legend
    leg=ROOT.TLegend(0.67,0.92,0.94,0.92-len(plottingSet)*0.06)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.05)
    leg.AddEntry( grColl['Data'], 'Data','f' )
    maxY,minY=-10,10
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSet:
        try:
            formatGraph(grColl[p],fill,color,marker,keepXUnc,shiftX)
            grColl[p].Draw(drawOpt)
            if p!='Data': leg.AddEntry(grColl[p],p,drawOpt)
            iminY,imaxY=getGraphExtremes(grColl[p])
            maxY=max(imaxY,maxY)
            minY=min(iminY,minY)
        except:
            pass
    frame.GetYaxis().SetRangeUser(minY*0.8,maxY*1.3)
    leg.Draw()
    
    #standard label
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.06)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.87,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.79,0.955,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

    p1.RedrawAxis()

    c.cd()
    p2=ROOT.TPad('p2','p2',0,0.0,1.0,0.4)
    p2.SetTopMargin(0.01)
    p2.SetRightMargin(0.03)
    p2.SetLeftMargin(0.12)
    p2.SetBottomMargin(0.3)
    p2.SetGridx()
    p2.Draw()
    p2.cd()
    ratioframe=frame.Clone('ratioframe')
    ratioframe.GetYaxis().SetTitle('Ratio to Data')
    ratioframe.GetYaxis().SetTitleOffset(0.7)
    ratioframe.GetYaxis().SetTitleSize(0.09)
    ratioframe.GetYaxis().SetLabelSize(0.08)
    ratioframe.GetXaxis().SetTitleOffset(1.6)
    ratioframe.GetXaxis().SetTitleSize(0.09)
    ratioframe.GetXaxis().SetLabelSize(0.09)
    ratioframe.GetXaxis().SetLabelOffset(0.02)
    ratioframe.Draw()
    minY,maxY=0.99,1.01
    if isPull:
        minY,maxY=3,-3
        ratioframe.GetYaxis().SetTitle('Pull')
        print 'Drawing as pull'
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSet:
        try:
            formatGraph(grCollComp[p],fill,color,marker,keepXUnc,shiftX)
            grCollComp[p].Draw(drawOpt)
            iminY,imaxY=getGraphExtremes(grCollComp[p])
            maxY=max(imaxY,maxY)
            minY=min(iminY,minY)
        except:
            pass
    if isPull:
        ratioframe.GetYaxis().SetRangeUser(-5,5)
    else:
        ratioframe.GetYaxis().SetRangeUser(max(minY*0.8,0.15),min(maxY*1.25,1.95))
    ratioframe.GetYaxis().SetNdivisions(5)

    p2.RedrawAxis()

    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/ueprofile%s.%s'%(outDir,pfix,ext))

    #only ratio
    cratio=ROOT.TCanvas('cratio','cratio',800,240)
    cratio.SetTopMargin(0.01)
    cratio.SetRightMargin(0.03)
    cratio.SetLeftMargin(0.12)
    cratio.SetBottomMargin(0.3)
    cratio.SetGridx()
    cratio.Draw()
    ratioframe.Draw()
    leg=ROOT.TLegend(0.15,0.92,0.95,0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.06)
    leg.SetNColumns(len(plottingSet)-1)
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSet:
        try:
            grCollComp[p].Draw(drawOpt)
            if p!='Data': leg.AddEntry(grCollComp[p],p,drawOpt)
        except:
            pass
    if isPull:
        ratioframe.GetYaxis().SetRangeUser(-5,5)
    else:
        ratioframe.GetYaxis().SetRangeUser(max(minY*0.8,0.15),min(maxY*1.5,1.95))
    ratioframe.GetYaxis().SetNdivisions(5)
    leg.Draw()

    cratio.RedrawAxis()
    cratio.Modified()
    cratio.Update()
    try:
        for ext in ['pdf','png']: cratio.SaveAs('%s/ueprofile%s_ratio.%s'%(outDir,pfix,ext))
    except:
        pass



def main():
    """
    configure and run
    """

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',
                      dest='input',
                      help='input directory [%default]',
                      default=None)
    parser.add_option('-s',
                      dest='slice',
                      help='slice type [%default]',
                      default=1,
                      type=int)
    (opt, args) = parser.parse_args()

    #set the slices to profile
    pfix,sliceList=SLICELISTS[opt.slice-1]
    print sliceList
    #fill the graphs with the mean values
    grColl={}
    finalSliceList=[]
    for np in xrange(0,len(sliceList)):
        var,varTitle=sliceList[np]
        
        #check if summary is available
        pckSummary=os.path.join(opt.input,var,'unfold/unfold_summary.pck')
        if not os.path.isfile(pckSummary) : 
            #print 'Skipping',var,' - failed to find',pckSummary
            continue

        print 'Keeping',var,' with unfolded histos @',pckSummary
        finalSliceList.append( (var,varTitle) )

        #read plots from summary
        with open(pckSummary,'r') as cachefile: 
            uePlots=pickle.load(cachefile)
            for key in uePlots:

                #init new graph if not yet available
                if not key in grColl:
                    grColl[key]=uePlots[key].plot[0].Clone()
                    grColl[key].Set(0)
                    
                #add mean and mean error
                mean=uePlots[key].mean[0]
                meanUnc=ROOT.TMath.Sqrt( sum(x*x for x in uePlots[key].mean[1]) )
                grColl[key].SetPoint(np,np+0.5,mean)
                grColl[key].SetPointError(np,0.5,meanUnc)

    #bail out if nothing found
    if len(finalSliceList)==0: return -1

    #compute the ratios
    grCollRatio,grCollPull=getProfileRatiosWithRespectTo(grColl,'Data')

    #show results
    obs=opt.input.split('/')[-1]
    #showProfile(grColl=grColl,grCollComp=grCollRatio,obs=obs,sliceList=finalSliceList,plottingSet=PLOTTINGSET_1,outDir=opt.input,pfix=pfix)
    #showProfile(grColl=grColl,grCollComp=grCollRatio,obs=obs,sliceList=finalSliceList,plottingSet=PLOTTINGSET_2,outDir=opt.input,pfix=pfix+'_v2')
    #showProfile(grColl=grColl,grCollComp=grCollRatio,obs=obs,sliceList=finalSliceList,plottingSet=PLOTTINGSET_3,outDir=opt.input,pfix=pfix+'_v3')
    showProfile(grColl=grColl,grCollComp=grCollPull,obs=obs,sliceList=finalSliceList,plottingSet=PLOTTINGSET_1,outDir=opt.input,pfix=pfix,isPull=True)
    showProfile(grColl=grColl,grCollComp=grCollPull,obs=obs,sliceList=finalSliceList,plottingSet=PLOTTINGSET_2,outDir=opt.input,pfix=pfix+'_v2',isPull=True)
    showProfile(grColl=grColl,grCollComp=grCollPull,obs=obs,sliceList=finalSliceList,plottingSet=PLOTTINGSET_3,outDir=opt.input,pfix=pfix+'_v3',isPull=True)

    

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
