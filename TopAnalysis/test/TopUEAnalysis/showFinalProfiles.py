import ROOT
import os
import sys
import optparse
import numpy as np
import pickle
from UEAnalysisHandler import VARTITLES,SYSTS
from UETools import getGraphExtremes,formatGraph
from showFinalDistributions import PLOTTINGSET_1,PLOTTINGSET_2,PLOTTINGSET_3,PLOTTINGSET_4

SLICELISTS=[
    ('',
     [ 
            ('inc',             'Inc.'),
            ('nj=0,1',          '=0'),
            ('nj=1,2',          '=1'),
            ('nj=2,999',        '#geq2'),
            ('ptll=0,20',       ']0,20['),
            ('ptll=20,40',      ']20,40['),
            ('ptll=40,80',      ']40,80['),
            ('ptll=80,120',     ']80,120['),
            ('ptll=120,9999',   '>120'),
            ('mll=0,60',        ']0,60['),
            ('mll=60,120',      ']60,120['),
            ('mll=120,200',     ']120,200['),
            ('mll=200,9999',    '>200'),
            ],
     [ 
            ('Extra jets',1,4),
            ("p_{T}(ll) / #it{GeV}",4,9),
            ("m(ll) / #it{GeV}",9,13)
            ]
     ),
    ('_ptll',[
            ('inc_ptll=tow',           'Inc.'),
            ('ptll=0,20_ptll=tow',     ']0,20['),
            ('ptll=20,40_ptll=tow',    ']20,40['),
            ('ptll=40,80_ptll=tow',    ']40,80['),
            ('ptll=80,120_ptll=tow',   ']80,120['),
            ('ptll=120,9999_ptll=tow', '#geq120'),
            ('inc_ptll=tra',           'Inc.'),
            ('ptll=0,20_ptll=tra',     ']0,20['),
            ('ptll=20,40_ptll=tra',    ']20,40['),
            ('ptll=40,80_ptll=tra',    ']40,80['),
            ('ptll=80,120_ptll=tra',   ']80,120['),
            ('ptll=120,9999_ptll=tra', '#geq120'),
            ('inc_ptll=awa',           'Inc.'),
            ('ptll=0,20_ptll=awa',     ']0,20['),
            ('ptll=20,40_ptll=awa',    ']20,40['),
            ('ptll=40,80_ptll=awa',    ']40,80['),
            ('ptll=80,120_ptll=awa',   ']80,120['),
            ('ptll=120,9999_ptll=awa', '#geq120'),
            ],
     [
            ('Toward',0,6),
            ('Transverse',6,12),
            ('Away',12,18)
            ]
     ),
    ('_ptllnj',
     [
            ('inc_ptll=tow',      'Inc.'),
            ('nj=0,1_ptll=tow',   '=0'),
            ('nj=1,2_ptll=tow',   '=1'),
            ('nj=2,999_ptll=tow', '#geq2'),
            ('inc_ptll=tra',      'Inc.'),
            ('nj=0,1_ptll=tra',   '=0'),
            ('nj=1,2_ptll=tra',   '=1'),
            ('nj=2,999_ptll=tra', '#geq2'),
            ('inc_ptll=awa',      'Inc.'),
            ('nj=0,1_ptll=awa',   '=0'),
            ('nj=1,2_ptll=awa',   '=1'),
            ('nj=2,999_ptll=awa', '#geq2'),
            ],
     [
            ('Toward',0,4),
            ('Transverse',4,8),
            ('Away',8,12)
            ]
     )
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
            ratioUnc=0 if yref**2==0 else ROOT.TMath.Sqrt(ratioUnc)/(yref**2)
            
            grCollRatios[key].SetPoint(np,x,ratio)
            grCollRatios[key].SetPointError(np,ex,ratioUnc)
            
            pull=(y-yref)/eyref if eyref>0 else 0
            epull=ey/eyref if eyref>0 else 0.
            if key==refKey : epull=1.0
            grCollPulls[key].SetPoint(np,x,pull)
            grCollPulls[key].SetPointError(np,ex,epull)

    return grCollRatios,grCollPulls


def showProfile(grColl,grCollComp,grCollStat,
                obs,sliceList,plottingSetList=[PLOTTINGSET_1],outDir='./',cmsLabel='#bf{CMS}',isPull=False,pfix='',categList=[]):

    """
    build the profile plot
    """
    #start the canvas
    c=ROOT.TCanvas('c','c',800,600)
    c.SetTopMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.2)
    #c.SetGridx()
    units=''
    if obs in ['chflux','chfluxz','chavgpt','chavgpt','chrecoil'] : units = ' [GeV]'
    
    frame=ROOT.TH1F('frame',';Category;< %s >%s;'%(VARTITLES[obs],units),len(sliceList),0,len(sliceList))
    for xbin in xrange(1,frame.GetNbinsX()+1): frame.GetXaxis().SetBinLabel(xbin,sliceList[xbin-1][1])
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleSize(0.06)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetLabelOffset(0.01)
    frame.GetXaxis().SetTitleOffset(1.5)
    frame.Draw()

    #plot and add to the legend
    #leg=ROOT.TLegend(0.67,0.86,0.94,0.86-len(plottingSetList[0])*0.05)
    leg=ROOT.TLegend(0.55,0.86,0.94,0.86-len(plottingSetList[0])*0.025)
    leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry( grColl['Data'], 'Data','f' )
    maxY,minY=-10,10
    for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSetList[0]:
        try:
            formatGraph(grColl[p],fill,color,marker,keepXUnc,shiftX)
            grColl[p].Draw(drawOpt)
            if p!='Data': leg.AddEntry(grColl[p],p,drawOpt)
            else:
                formatGraph(grCollStat[p],3104,1,marker,keepXUnc,shiftX)
                grCollStat[p].Draw(drawOpt)
            iminY,imaxY=getGraphExtremes(grColl[p])
            maxY=max(imaxY,maxY)
            minY=min(iminY,minY)
        except:
            pass
    frame.GetYaxis().SetRangeUser(minY*0.8,maxY*1.3)
    leg.Draw()

    #stat component for data                                                                                                                                                                                                     
    txt=ROOT.TPaveText(0.555,0.83,0.59,0.84,'brNDC')
    txt.SetFillStyle(3004)
    txt.SetFillColor(1)
    txt.SetBorderSize(0)
    txt.SetLineColor(0)
    txt.Draw()

    #category labels
    catTxt=ROOT.TLatex()
    catTxt.SetTextFont(52)
    catTxt.SetTextSize(0.035)
    catTxt.SetTextAlign(22)
    catTxt.SetNDC(False)
    ar=ROOT.TArrow()
    ar.SetNDC(False)
    ar.SetLineWidth(2)
    for t,xmin,xmax in categList:
        yarr=1.315*maxY
        ytxt=1.36*maxY
        if obs=='detST': yarr,ytxt=0.0415,0.042 #0.44,0.45
        if obs=='sphericity': yarr,ytxt=0.61,0.62 #0.44,0.45
        if obs=='aplanarity': yarr,ytxt=0.188,0.192 #0.1238,0.1268
        if obs=='C': yarr,ytxt=0.88,0.895 #0.66,0.675
        if obs=='D': yarr,ytxt=0.47,0.48 #0.271,0.278
        catTxt.DrawLatex(0.5*(xmax+xmin),ytxt,t)
        ar.DrawArrow(xmin,yarr,xmax,yarr,0.01,"<>")

    #standard label
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.06)
    tex.SetNDC()
    tex.DrawLatex(0.16,0.95,cmsLabel)
    tex.DrawLatex(0.7,0.95,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')

    c.RedrawAxis()
    c.Modified()
    c.Update()
    for ext in ['png','pdf','root']:
        c.SaveAs('%s/ueprofile%s.%s'%(outDir,pfix,ext))

    #start the ratio canvas
    dy_xtit=20.
    dy_pad=180.
    y_total=4*dy_xtit+dy_pad*len(plottingSetList)
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
    ar=ROOT.TArrow()
    ar.SetNDC(False)
    ar.SetLineWidth(2)
    for i in xrange(0,len(plottingSetList)):
        cratio.cd()
        
        #start new sub-pad
        dy_ipad=rel_dy_pad
        if i==0 : dy_ipad += 3*rel_dy_xtit
        if i==len(plottingSetList)-1: dy_ipad+=rel_dy_xtit
        sp.append( ROOT.TPad('p%d'%i,'p%d'%i,0,y,1.0,max(y-dy_ipad,0.)) )
        y=max(y-dy_ipad,0.)

        #margins
        if i==0:
            sp[-1].SetTopMargin(0.33)
        else:
            sp[-1].SetTopMargin(0.11)
        sp[-1].SetRightMargin(0.25)
        sp[-1].SetLeftMargin(0.12)
        if i==len(plottingSetList)-1:
            sp[-1].SetTopMargin(0.06)
            sp[-1].SetBottomMargin(0.28) #2*dy_xtit/dy_ipad)
        else:
            sp[-1].SetBottomMargin(0.11)
        #sp[-1].SetGridx()
        sp[-1].Draw()
        sp[-1].cd()        

        rf.append( frame.Clone('ratioframe') )
        rf[-1].GetYaxis().SetTitle('Theory/Data')
        rf[-1].GetYaxis().SetTitleOffset(0.3)
        rf[-1].GetYaxis().SetTitleSize(0.12)
        rf[-1].GetYaxis().SetLabelSize(0.11)
        if i==len(plottingSetList)-1:
            rf[-1].GetXaxis().SetTitleSize(0.16)
            rf[-1].GetXaxis().SetLabelSize(0.11)
            rf[-1].GetXaxis().SetLabelOffset(0.01)
            rf[-1].GetXaxis().SetTitleOffset(0.8)
            rf[-1].GetXaxis().SetTitle('Category')
        else:
            rf[-1].GetXaxis().SetTitle('')
            for xbin in xrange(0,rf[-1].GetNbinsX()): rf[-1].GetXaxis().SetBinLabel(xbin+1,'')
        if i==0:
            rf[-1].GetYaxis().SetTitleOffset(0.4)
            rf[-1].GetYaxis().SetTitleSize(0.095)
            rf[-1].GetYaxis().SetLabelSize(0.088)
        

        minY,maxY=0.99,1.01
        if isPull:
            minY,maxY=3,-3
            rf[-1].GetYaxis().SetTitle('#frac{(Theory-Data)}{#sigma_{Data}}')
            print 'Drawing as pull'
        rf[-1].Draw()

        if i==0:
            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.13)
            tex.SetNDC()
            tex.DrawLatex(0.12,0.84,cmsLabel)
            tex.DrawLatex(0.5,0.84,'#scale[0.9]{35.9 fb^{-1} (13 TeV)}')
        
        if i==0:
            lg.append( ROOT.TLegend(0.75,0.1,0.99,0.7) )
        elif i==len(plottingSetList)-1:
            lg.append( ROOT.TLegend(0.75,0.23,0.99,0.95) )
        else:
            lg.append( ROOT.TLegend(0.75,0.1,0.99,0.95) )
        lg[-1].SetFillStyle(0)
        lg[-1].SetBorderSize(0)
        lg[-1].SetTextFont(42)
        lg[-1].SetTextSize(0.085 if i==0 else 0.1)
        #lg[-1].SetNColumns(len(plottingSetList[i])-1)
    
        for p,drawOpt,fill,color,marker,keepXUnc,shiftX in plottingSetList[i]:
            try:
                formatGraph(grCollComp[p],fill,color,marker,keepXUnc,shiftX)
                grCollComp[p].Draw(drawOpt)
                iminY,imaxY=getGraphExtremes(grCollComp[p])
                maxY=max(imaxY,maxY)
                minY=min(iminY,minY)
                if p!='Data':
                    legOpt=drawOpt if not drawOpt in ['c','2'] else 'l'
                    lg[-1].AddEntry(grCollComp[p],p,legOpt)
            except:
                pass

        #re-adapt range according to limits found
        if isPull:
            rf[-1].GetYaxis().SetRangeUser(-5,5)
        else:
            rf[-1].GetYaxis().SetRangeUser(max(minY*0.8,0.15),min(maxY*1.25,1.95))
        rf[-1].GetYaxis().SetNdivisions(5)

        sp[-1].RedrawAxis()
        lg[-1].Draw()
    
        #category labels
        if i==0:
            catTxt=ROOT.TLatex()
            catTxt.SetTextFont(52)
            catTxt.SetTextSize(0.08)
            catTxt.SetTextAlign(22)
            catTxt.SetNDC(False)
            maxYLabel=6.5 if isPull else maxY*1.35 
            maxYArr=5.5 if isPull else maxY*1.125
            for t,xmin,xmax in categList:
                catTxt.DrawLatex(0.5*(xmax+xmin),maxYLabel,t)
                ar.DrawArrow(xmin,maxYArr,xmax,maxYArr,0.01,"<>")

    # all done
    cratio.Modified()
    cratio.Update()
    pfix+='_pull' if isPull else '_ratio'
    for ext in ['pdf','png','root']: cratio.SaveAs('%s/ueprofile%s.%s'%(outDir,pfix,ext))


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
    parser.add_option('--doPull',
                      dest='doPull',
                      help='do pulls instead of ratios [%default]',
                      default=False,
                      action='store_true')
    parser.add_option('-s',
                      dest='slice',
                      help='slice type [%default]',
                      default=1,
                      type=int)
    parser.add_option('--cmsLabel',
                      dest='cmsLabel',
                      help='cmsLabel [%default]',
                      default='#bf{CMS} #it{preliminary}',
                      type='string')
    (opt, args) = parser.parse_args()

    #set the slices to profile
    pfix,sliceList,categList=SLICELISTS[opt.slice-1]

    #fill the graphs with the mean values
    grColl={}
    grCollStat={}
    finalSliceList=[]
    for np in xrange(0,len(sliceList)):
        var,varTitle=sliceList[np]
        
        #check if summary is available
        pckSummary=os.path.join(opt.input,var,'unfold/unfold_summary.pck')
        if not os.path.isfile(pckSummary) : 
            print 'Skipping',var,' - failed to find',pckSummary
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
                    grCollStat[key]=grColl[key].Clone( grColl[key].GetName()+'_stat')

                #add mean and mean error
                mean=uePlots[key].mean[0]
                meanUnc=ROOT.TMath.Sqrt( sum(x*x for x in uePlots[key].mean[1]) )
                grColl[key].SetPoint(np,np+0.5,mean)
                grColl[key].SetPointError(np,0.5,meanUnc)
                
                grCollStat[key].SetPoint(np,np+0.5,mean)
                grCollStat[key].SetPointError(np,0.5,uePlots[key].mean[1][0])

    #bail out if nothing found
    if len(finalSliceList)==0: return -1

    #compute the ratios
    grCollRatio,grCollPull=getProfileRatiosWithRespectTo(grColl,'Data')
    grCollRatioStat,_=getProfileRatiosWithRespectTo(grCollStat,'Data')
    
    #show results
    obs=opt.input.split('/')[-1]    
    if len(obs)==0: obs=opt.input.split('/')[-2]
    showProfile(grColl=grColl,grCollComp=grCollPull,grCollStat=grCollStat,
                obs=obs,sliceList=finalSliceList,
                plottingSetList=[PLOTTINGSET_1,PLOTTINGSET_2,PLOTTINGSET_3], #PLOTTINGSET_4],
                outDir=opt.input,isPull=opt.doPull,cmsLabel=opt.cmsLabel,pfix=pfix,
                categList=categList)


    

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
