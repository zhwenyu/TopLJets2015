import ROOT
import os
import sys
import optparse
import numpy as npy
import pickle
from collections import defaultdict
from UEAnalysisHandler import VARTITLES
from UETools import formatGraph

def performChisquareFitTo(paramScan):
    """
    scans the chi^2 values, determines minimum and appropriate uncertainty range
    """

    #turn into a graph
    gr=ROOT.TGraph()
    gr.SetMarkerStyle(20)
    for pval,chi2 in paramScan:
        n=gr.GetN()
        gr.SetPoint(n,pval,chi2)
    gr.Sort()
    
    #finer scan with TSpline3
    xy=[]
    for xval in npy.arange(gr.GetX()[0], gr.GetX()[gr.GetN()-1], 0.0001):
        yval = gr.Eval(xval, 0, 'S')
        xy.append( [xval,yval] )
    
    #chi^2 minimum
    y0=min( [ pt[1] for pt in xy ] )
    x0= [ pt[0] for pt in xy if pt[1]==y0 ][0]

    #scans for values where delta chi^2=1
    xyl=[ pt for pt in xy if pt[1]>y0+1 and pt[0]<x0 ]
    xmin = xyl[-1][0] if len(xyl)>0 else None
    xyr=[ pt for pt in xy if pt[1]>y0+1 and pt[0]>x0 ]
    xmax = xyr[0][0] if len(xyr)>0 else None

    #scan for values where chi^2->2chi^2
    xyl_2=[ pt for pt in xy if pt[1]>y0*2 and pt[0]<x0 ]
    xmin_2 = xyl_2[-1][0] if len(xyl_2)>0 else None
    xyr_2=[ pt for pt in xy if pt[1]>y0*2 and pt[0]>x0 ]
    xmax_2 = xyr_2[0][0] if len(xyr_2)>0 else None

    #maximize
    if xmin_2 : xmin=min(xmin,xmin_2) if xmin else xmin_2
    if xmax_2 : xmax=max(xmax,xmax_2) if xmax else xmax_2

    #return the result
    return (y0,x0,xmax,xmin,gr)


def buildChisquareReportFrom(pckSummary):
    
    """
    Opens the unfold summary and computes chi^2 and p-values between the data and the models
    """

    chi2report={'mean':{},'dist':{}}
    chi2Scan=defaultdict(list)

    #check if it is valid
    if not os.path.isfile(pckSummary) : return chi2report,chi2Scan

    with open(pckSummary,'r') as cachefile: 
        uePlots=pickle.load(cachefile)

        #data
        data=uePlots['Data']
        np=data.plot[0].GetN()


        for model in uePlots:
            if model=='Data' : continue
            modelPlot=uePlots[model]

            #sum up covariance matrices and invert
            sumCov=ROOT.TMatrixF(data.covMatrices['total'])
            if 'total' in modelPlot.covMatrices: sumCov+=modelPlot.covMatrices['total']
            invTotalCov=ROOT.TMatrixF(sumCov)
            invTotalCov.Invert()

            sumNoModelSystCov=ROOT.TMatrixF(data.covMatrices['total'])
            invNoModelSystCov=ROOT.TMatrixF(sumNoModelSystCov)
            invNoModelSystCov.Invert()

            #vector of differences
            diffVec=ROOT.TVectorF(np)
            x,y=ROOT.Double(0),ROOT.Double(0)
            for i in xrange(0,np):

                data.plot[0].GetPoint(i,x,y)
                ydata_i=float(y)

                modelPlot.plot[0].GetPoint(i,x,y)
                ymodel_i=float(y)

                #the difference has to be multiplied by the bin width 
                ex=data.plot[0].GetErrorX(i)
                diffVec[i]=ex*(ymodel_i-ydata_i)

            #chi^2 for distribution analysis
            chi2,chi2nomodelsyst=0,0
            for i in xrange(0,np):
                for j in xrange(0,np):
                    chi2            += diffVec[i]*invTotalCov[i][j]*diffVec[j]            
                    chi2nomodelsyst += diffVec[i]*invNoModelSystCov[i][j]*diffVec[j]            
            pval=ROOT.TMath.Prob(chi2,np)
            chi2report['dist'][model +'*']=(chi2,np,pval)
            pvalnomodelsyst=ROOT.TMath.Prob(chi2nomodelsyst,np)
            chi2report['dist'][model]=(chi2nomodelsyst,np,pvalnomodelsyst)

            if '/inc/' in pckSummary:
                print pckSummary.split('/')[2],model,(chi2,np,pval),'|',(chi2nomodelsyst,np,pvalnomodelsyst)

            if '#alpha_{S}' in model:
                param,valStr=model.split('=')
                val=float(valStr)
                chi2Scan[param].append( (val,chi2/np) )
 

            #chi^2 from mean analysis
            mean=modelPlot.mean[0]
            meanUnc=ROOT.TMath.Sqrt( sum(x*x for x in modelPlot.mean[1]) )
            meanData=data.mean[0]
            meanDataUnc=ROOT.TMath.Sqrt( sum(x*x for x in data.mean[1]) )
            pull=(mean-meanData)/ROOT.TMath.Sqrt(meanUnc**2+meanDataUnc**2)
            pval=ROOT.TMath.Prob(pull**2,1)
            chi2report['mean'][model +'*']=(pull,1,pval)

            meanStatUnc=modelPlot.mean[1][0]
            pullNoModelSyst=(mean-meanData)/ROOT.TMath.Sqrt(meanStatUnc**2+meanDataUnc**2)
            pvalNoModelSyst=ROOT.TMath.Prob(pullNoModelSyst**2,1)
            chi2report['mean'][model]=(pullNoModelSyst,1,pvalNoModelSyst)

    chi2ScanReport={}
    for param in chi2Scan:
        chi2ScanReport[param]=performChisquareFitTo(chi2Scan[param])

    return chi2report,chi2ScanReport


def main():
    """
    configure and run
    """

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',
                      dest='input',
                      help='input directory [%default]',
                      default=None)
    parser.add_option('--cmsLabel',
                      dest='cmsLabel',
                      help='cms label [%default]',
                      default='#bf{CMS} #it{preliminary}')
    (opt, args) = parser.parse_args()

        
    #check if summary is available
    chi2Results={}
    paramScanResults={}
    varList=[]
    for var in os.listdir(opt.input):

        fulld=os.path.join(opt.input,var)
        if not os.path.isdir(fulld) : continue

        varList.append(var)
        varKey='evshape' if var in ['C','D','sphericity','aplanarity'] else 'flux'
        for psSlice in os.listdir(fulld):

            #skip these ones for the moment
            if 'chmult' in psSlice: continue

            pckSummary=os.path.join(fulld,psSlice,'unfold/unfold_summary.pck')

            chi2report,paramScan=buildChisquareReportFrom(pckSummary)

            for ana in chi2report:
                for model in chi2report[ana]:
                    if not model in chi2Results: chi2Results[model]={}
                    anaKey=(ana,varKey)
                    if not anaKey in chi2Results[model]: chi2Results[model][anaKey]=[]
                    chi2Results[model][anaKey].append(chi2report[ana][model])

            for param in paramScan:
                if not param in paramScanResults: paramScanResults[param]={}
                paramScanResults[param][(var,psSlice)]=paramScan[param]
    

    #display the results
    models2Plot=['PW+PY8*','PW+PY8',
                 'ISR up','ISR dn','FSR up','FSR dn','ERD on','QCD based','Gluon move','UE up','UE dn','no MPI','no CR',
                 'aMC@NLO+PY8','PW+HW++','PW+HW7','Sherpa']
    ana2Plot=[('mean','evshape'),('dist','evshape'),('mean','flux'),('dist','flux')]

    c=ROOT.TCanvas('c','c',550,500)
    c.SetTopMargin(0.06)
    c.SetLeftMargin(0.23)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.15)
    c.SetGridy()


    for anaKey in ana2Plot:
        ana,varKey=anaKey

        summaryH=None
        ztitle='event shape' if varKey=='evshape' else 'flux'
        if ana=='dist':
            summaryH=ROOT.TH2F('pvalSummary', 
                               ';p-value;;Number of '+ztitle+' analysis', 
                               20,1e-3,1.0, 
                               len(models2Plot),0,len(models2Plot))        
        else:
            summaryH=ROOT.TH2F('pullSummary', 
                               ';Pull;;Number of '+ztitle+' analysis', 
                               20,-5,5, 
                               len(models2Plot),0,len(models2Plot))        

        for i in xrange(0,len(models2Plot)):
            ybin=len(models2Plot)-i-1
            model=models2Plot[i]
            if not model in chi2Results: continue
            summaryH.GetYaxis().SetBinLabel(ybin+1,model)
            for chi2,np,pval in chi2Results[model][anaKey]:
                if ana=='dist':
                    summaryH.Fill(min(0.999,max(pval,1e-3)),ybin)
                else:
                    summaryH.Fill(min(4.999,max(chi2,-5.0)),ybin)
                            
        c.Clear()            
        summaryH.Draw('colz')
        summaryH.GetXaxis().SetLabelSize(0.04)
        summaryH.GetXaxis().SetTitleSize(0.05)
        summaryH.GetXaxis().SetTitleOffset(0.88)
        summaryH.GetYaxis().SetLabelSize(0.05)
        if ana=='mean':
            line=ROOT.TLine()
            line.SetLineWidth(3)
            line.SetLineStyle(2)
            line.SetLineColor(1)
            line.DrawLine(0,0,0,len(models2Plot))
            line.DrawLine(1,0,1,len(models2Plot))
            line.DrawLine(-1,0,-1,len(models2Plot))
            
            gr=ROOT.TGraphErrors()
            gr.SetLineWidth(2)
            gr.SetMarkerStyle(20)
            pullLabel=ROOT.TLatex()
            pullLabel.SetTextFont(42)
            pullLabel.SetTextSize(0.03)
            pullLabel.SetTextAlign(12)
            for ybin in xrange(1,summaryH.GetNbinsY()+1):
                proj=summaryH.ProjectionX("py",ybin,ybin)
                mean,rms=proj.GetMean(),proj.GetRMS()
                ycen=summaryH.GetYaxis().GetBinCenter(ybin)
                gr.SetPoint(ybin-1,mean,ycen)
                gr.SetPointError(ybin-1,rms,0)
                pullLabel.DrawLatex(-4.5,ycen,
                                     '#scale[0.8]{%3.2f #pm %3.2f}'%(mean,rms))
                proj.Delete()
            gr.Draw('p')
            
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.05)
        tex.SetNDC()
        tex.DrawLatex(0.1,0.96,opt.cmsLabel)
        tex.DrawLatex(0.65,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
        c.RedrawAxis()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/anasummary_%s_%s.%s'%(opt.input,varKey,ana,ext))
            
        summaryH.Delete()
    
    #show the parameter scans now
    c.Clear()
    c.SetTopMargin(0.06)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetGridy(False)
    for p in paramScanResults:

        pname='asfsr'
        if 'ISR' in p : pname='asisr'

        evshapeFits=ROOT.TGraphAsymmErrors()
        evshapeFits.SetMarkerStyle(20)
        evshapeFits.SetName('evshapes')
        fluxFits=evshapeFits.Clone('flux')
        for v,psSlice in paramScanResults[p]:

            chi20,p0,pmax,pmin,gr=paramScanResults[p][(v,psSlice)]
            if pmax and pmin:
                summarygr=evshapeFits if v in ['C','D','sphericity','aplanarity'] else fluxFits
                npts=summarygr.GetN()
                summarygr.SetPoint(npts,p0,npts)
                summarygr.SetPointError(npts,pmax-p0,p0-pmin,0,0)

            if psSlice!='inc' : continue

            #show chi2 scan plots with different comparisons
            def addToLegend(leg,chi2Res):
                pVal,p0Val,pmaxVal,pminVal,chi2Gr=chi2Res
                legTxt='#splitline{%s}{#scale[0.7]{%3.3f '%(chi2Gr.GetTitle(),p0Val)
                if pminVal and pmaxVal:
                    legTxt+='[%3.3f,%3.3f]}}'%(pminVal,pmaxVal)
                elif pminVal:
                    legTxt+='[%3.3f,n/a]}}'%pminVal
                elif pmaxVal:
                    legTxt+='[n/a,%3.3f]}}'%pmaxVal
                else:
                    legTxt+=' n/a}}'
                leg.AddEntry(chi2Gr,legTxt,'lp')
                
            cSlicesList={'ptll':[('inc_ptll=awa','away'), ('inc_ptll=tow','toward'),('inc_ptll=tra','transverse')],
                         'nj':[('nj=0,1','N_{j}=0'),('nj=1,2','N_{j}=1'),('nj=2,999','N_{j}#geq2')],
                         'inc':[('chavgpz','#bar{p}_{z}'),('aplanarity','Aplanarity'),('sphericity','Sphericity')]

                         }            
            for pfix in cSlicesList:
                if pfix=='ptll' and v in ['C','D','sphericity','aplanarity']: continue
                if pfix=='inc'  and v!='chavgpt': continue
                cSlices=cSlicesList[pfix]

                c.Clear()
                gr.Draw('apc')
                gr.SetTitle('inclusive')
                gr.GetXaxis().SetTitle(p)
                gr.GetYaxis().SetTitle('#chi^{2}/ndf')
                gr.GetYaxis().SetRangeUser(0,max(chi20*4,chi20+4))
                gr.GetXaxis().SetTitleSize(0.05)
                gr.GetYaxis().SetTitleSize(0.05)
                gr.GetXaxis().SetLabelSize(0.04)
                gr.GetYaxis().SetLabelSize(0.04)
                gr.SetLineColor(1)
                gr.SetMarkerColor(1)
                gr.SetLineWidth(2)

                leg=ROOT.TLegend(0.16,0.8,0.45,0.54) if pfix=='inc' else ROOT.TLegend(0.16,0.76,0.45,0.32)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetTextFont(42)
                leg.SetTextSize(0.045)
                if pfix=='inc' : leg.AddEntry(gr,VARTITLES[v],'lp')
                else : addToLegend(leg,paramScanResults[p][(v,psSlice)])
                
                #draw comparisons
                compColors=['#92c5de','#f4a582','#ca0020']
                for ic in xrange(0,len(cSlices)):
                    cSlice,cTitle=cSlices[ic]
                    v2comp=v
                    if pfix=='inc' :
                        v2comp=cSlice
                        cSlice=psSlice
                    if not (v2comp,cSlice) in paramScanResults[p] : continue
                    paramScanResults[p][(v2comp,cSlice)][4].SetTitle(cTitle)
                    ci=ROOT.TColor.GetColor(compColors[ic])
                    paramScanResults[p][(v2comp,cSlice)][4].SetLineColor(ci)
                    paramScanResults[p][(v2comp,cSlice)][4].SetMarkerColor(ci)
                    paramScanResults[p][(v2comp,cSlice)][4].SetLineWidth(2)
                    paramScanResults[p][(v2comp,cSlice)][4].SetMarkerStyle(24+ic)
                    paramScanResults[p][(v2comp,cSlice)][4].Draw('pc')
                    if pfix=='inc' : leg.AddEntry(paramScanResults[p][(v2comp,cSlice)][4],cTitle,'lp')
                    else : addToLegend(leg,paramScanResults[p][(v2comp,cSlice)])
                
                leg.Draw()
            
                #draw aux lines for the range
                if pfix!='inc':
                    line.SetLineColor(ROOT.kRed)
                    if pmax:
                        line.DrawLine(pmax,0,pmax,gr.Eval(pmax))
                    if pmin:
                        line.DrawLine(pmin,0,pmin,gr.Eval(pmin))
               
                tex=ROOT.TLatex()
                tex.SetTextFont(42)
                tex.SetTextSize(0.05)
                tex.SetNDC()
                tex.DrawLatex(0.16,0.88,opt.cmsLabel)
                if pfix!='inc' : tex.DrawLatex(0.16,0.8,VARTITLES[v])
                tex.DrawLatex(0.65,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

                c.RedrawAxis()
                c.Modified()
                c.Update()
                for ext in ['png','pdf']: c.SaveAs('%s/chi2scans_%s_%s_%s.%s'%(opt.input,v,pname,pfix,ext))

        for gr in [evshapeFits,fluxFits]:
            c.Clear()

            frame=ROOT.TH1F('frame','frame',30,0.05,0.25)
            frame.Draw()
            frame.GetYaxis().SetNdivisions(0)
            frame.GetYaxis().SetTitle('Analysis')
            frame.GetXaxis().SetTitle(p)  
            frame.GetXaxis().SetTitleSize(0.05)
            frame.GetYaxis().SetTitleSize(0.05)
            frame.GetXaxis().SetLabelSize(0.04)
            frame.GetYaxis().SetLabelSize(0.04)
            frame.GetYaxis().SetRangeUser(-1,gr.GetN()*1.2)

            gr.Draw('p')

            fitRes=[]
            for i in xrange(0,gr.GetN()) : fitRes.append( gr.GetX()[i] )
            try:
                perc=npy.percentile( npy.array(fitRes), [16,50,84] )
                line.SetLineColor(ROOT.kRed)
                line.DrawLine(perc[0],0,perc[0],gr.GetN())
                line.DrawLine(perc[2],0,perc[2],gr.GetN())
            except:
                print 'Unable to retrieve quantiles for',gr.GetName(),pname

            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.05)
            tex.SetNDC()
            tex.DrawLatex(0.16,0.88,opt.cmsLabel)
            tex.DrawLatex(0.6,0.88,'#scale[0.8]{%s=%3.3f}'%(p,perc[1]))
            tex.DrawLatex(0.6,0.8,'#scale[0.8]{[%3.3f,%3.3f]}'%(perc[0],perc[2]))
            tex.DrawLatex(0.65,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
            c.RedrawAxis()
            c.Modified()
            c.Update()
            for ext in ['png','pdf']: c.SaveAs('%s/chi2summary_%s_%s.%s'%(opt.input,gr.GetName(),pname,ext))

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
