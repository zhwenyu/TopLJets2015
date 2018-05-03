import ROOT
import os
import sys
import optparse
import numpy as npy
import pickle
from collections import defaultdict
from UEAnalysisHandler import VARTITLES
from UETools import formatGraph

def performChisquareFitTo(paramScan,bySpline=False):
    """
    scans the chi^2 values, determines minimum and appropriate uncertainty range
    """

    #turn into a graph
    gr=ROOT.TGraph()
    gr.SetMarkerStyle(20)
    minpval,minchi2=0.120,999999.
    for pval,chi2 in paramScan:
        n=gr.GetN()
        if chi2>100 or chi2<0: continue
        gr.SetPoint(n,pval,chi2)
        if chi2<minchi2:
            minchi2=chi2
            minpval=pval
    gr.Sort()

    #finer scan with TSpline3 or pol4 fit
    xy=[]
    if bySpline:
        for xval in npy.arange(gr.GetX()[0], gr.GetX()[gr.GetN()-1], 0.0001):
            yval = gr.Eval(xval, 0, 'S')
            xy.append( [xval,yval] )
    else:

        bestPol=(2,None,None,None)
        for ipol in xrange(2,9):
            gr.Fit('pol%d'%ipol,'FCQR0') 
            polFnc=gr.GetFunction('pol%d'%ipol)        
            sse,dof=polFnc.GetChisquare(),polFnc.GetNDF()
            if sse==0 or dof==0 : continue

            fdistI=None
            prevBestPol,prevSSE,prevDOF,bestipolfdistI=bestPol
            if prevSSE:
                f=((prevSSE-sse)/(prevDOF-dof)) / (sse/dof)
                fdistI=1-ROOT.TMath.FDistI(f,prevDOF-dof,dof)
                if fdistI<0.1:
                    bestPol=(ipol,sse,dof,fdistI)
            else:
                bestPol=(ipol,sse,dof,fdistI)

            print ipol,sse,dof,fdistI,'(',prevBestPol,')'


        polOrd=bestPol[0]
        print '=>chose',polOrd
        gr.Fit('pol%d'%polOrd,'FCQR0') #'FCQR'
        polFnc=gr.GetFunction('pol%d'%polOrd)

        for xval in npy.arange(gr.GetX()[0], gr.GetX()[gr.GetN()-1], 0.0001):
            yval=polFnc.Eval(xval)
            xy.append( [xval,yval] )

        for i in xrange(0,len(xy)):
            gr.SetPoint(i,xy[i][0],xy[i][1])


    #chi^2 minimum
    y0=min( [ pt[1] for pt in xy ] )
    x0= [ pt[0] for pt in xy if pt[1]==y0 ][0]

    #scans for values where delta chi^2=1 (1sigma)
    xyl=[ pt for pt in xy if pt[1]>y0+1 and pt[0]<x0 ]
    xmin = xyl[-1][0] if len(xyl)>0 else None
    xyr=[ pt for pt in xy if pt[1]>y0+1 and pt[0]>x0 ]
    xmax = xyr[0][0] if len(xyr)>0 else None

    #scan for values where delta chi^2=2 (2sigma)
    xyl_2=[ pt for pt in xy if pt[1]>y0+4 and pt[0]<x0 ]
    xmin_2 = xyl_2[-1][0] if len(xyl_2)>0 else None
    xyr_2=[ pt for pt in xy if pt[1]>y0+4 and pt[0]>x0 ]
    xmax_2 = xyr_2[0][0] if len(xyr_2)>0 else None


    #return the result
    return (y0,x0,xmax,xmin,xmin_2,xmax_2,gr)


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
            pval=ROOT.TMath.Prob(chi2,np-1)
            chi2report['dist'][model +'*']=(chi2,np,pval)
            pvalnomodelsyst=ROOT.TMath.Prob(chi2nomodelsyst,np)
            chi2report['dist'][model]=(chi2nomodelsyst,np,pvalnomodelsyst)

            if '/inc/' in pckSummary:
                print pckSummary.split('/')[2],model,(chi2,np,pval),'|',(chi2nomodelsyst,np,pvalnomodelsyst)

            if '#alpha_{S}' in model:
                param,valStr=model.split('=')
                val=float(valStr)
                if val==0.1365 : continue
                #chi2Scan[param].append( (val,chi2/np) )
                chi2Scan[param].append( (val,chi2) )
                print model,param,valStr,chi2Scan[param]

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
        print pckSummary,param
        chi2ScanReport[param]=performChisquareFitTo(chi2Scan[param])

    return chi2report,chi2ScanReport


def main():
    """
    configure and run
    """

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    #ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
    # Viridis palette reversed + white
    from array import array
    stops = array('d', [0.0, 0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0000])
    red   = array('d', [26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255., 1., 1.])
    green = array('d', [9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255., 1., 1.])
    blue  = array('d', [30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255., 1., 1.])
    ROOT.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 255)

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
        if not var in ['chavgpt','C','C_2'] : continue
        if var in ['maxRap','rapDist']: 
            print 'Skipping',var
            continue
    

        varList.append(var)
        varKey='flux'
        if var in ['C','D','sphericity','aplanarity']: varKey='evshape'
        if var in ['C_2','D_2','sphericity_2','aplanarity_2']: varKey='evshape_2'
        for psSlice in os.listdir(fulld):

            #skip these ones for the moment
            if 'chmult' in psSlice: continue
            if not 'inc' in psSlice: continue
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
                 'ISR up','ISR dn','FSR up','FSR dn','ERD on','QCD based','Gluon move','Rope','Rope (no CR)', 'UE up','UE dn','no MPI','no CR',
                 'aMC@NLO+PY8','PW+HW++','PW+HW7','Sherpa']
    ana2Plot=[('mean','evshape'),('dist','evshape'),
              ('mean','evshape_2'),('dist','evshape_2'),
              ('mean','flux'),('dist','flux')]

    c=ROOT.TCanvas('c','c',550,500)
    c.SetTopMargin(0.06)
    c.SetLeftMargin(0.23)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.15)
    c.SetGridy()


    for anaKey in ana2Plot:
        ana,varKey=anaKey

        summaryH=None
        ztitle='flux'
        if varKey=='evshape' : ztitle='event shape'
        if varKey=='evshape_2': ztitle='event shape (quadratic)'
        if ana=='dist':
            summaryH=ROOT.TH2F('pvalSummary', 
                               #';p-value;;Number of '+ztitle+' analysis', 
                               ';p-value;;Number of distributions*',
                               20,1e-3,1.0, 
                               len(models2Plot),0,len(models2Plot))        
        else:
            summaryH=ROOT.TH2F('pullSummary', 
                               #';Pull;;Number of '+ztitle+' analysis', 
                               ';Pull;;Number of distributions*',
                               20,-5,5, 
                               len(models2Plot),0,len(models2Plot))        

        for i in xrange(0,len(models2Plot)):
            ybin=len(models2Plot)-i-1
            model=models2Plot[i]
            if not model in chi2Results: continue
            summaryH.GetYaxis().SetBinLabel(ybin+1,model)
            if not anaKey in chi2Results[model]:
                print 'Failed to find',anaKey,'for model=',model
                continue
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
        summaryH.GetZaxis().SetTitleOffset(1.3)
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
        tex.DrawLatex(0.67,0.96,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')
        varList=', '.join( [VARTITLES['chmult'],VARTITLES['chflux'],VARTITLES['chavgpt'],VARTITLES['chfluxz'],VARTITLES['chavgpz'],VARTITLES['chrecoil']] )
        if not 'flux' in varKey:
            varList=', '.join([VARTITLES['sphericity'],VARTITLES['aplanarity'],VARTITLES['C'],VARTITLES['D']])
        tex.DrawLatex(0.25,0.91,'#scale[0.6]{* - %s}'%varList)
        c.RedrawAxis()
        c.Modified()
        c.Update()
        for ext in ['png','pdf','root']: c.SaveAs('%s/anasummary_%s_%s.%s'%(opt.input,varKey,ana,ext))
            
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
        if 'ISR' in p      : pname='asisr'
        if 'CMW' in p      : pname += '_cmw2loop'
        if '0.5M_{Z}' in p : pname +=' scale0.5'
        if 'M_{Z}' in p    : pname +=' scale1.0'
        if '2M_{Z}' in p   : pname +=' scale2.0'

        evshapeFits=ROOT.TGraphAsymmErrors()
        evshapeFits.SetMarkerStyle(20)
        evshapeFits.SetName('evshapes')
        evshapeFits2s=evshapeFits.Clone('evshapes2s')
        evshapeFits2s.SetMarkerStyle(1)
        evshapeFits2s.SetLineColor(2)
        evshapeFits2s.SetMarkerColor(2)
        evshape_2Fits=evshapeFits.Clone('evshapes_2')
        evshape_2Fits2s=evshapeFits2s.Clone('evshapes_22s')
        fluxFits=evshapeFits.Clone('flux')
        fluxFits2s=evshapeFits2s.Clone('flux2s')
        for v,psSlice in paramScanResults[p]:

            chi20,p0,pmax,pmin,pmax2s,pmin2s,gr=paramScanResults[p][(v,psSlice)]
            if pmax and pmin:
                summarygr,summarygr2s=fluxFits,fluxFits2s
                if v in ['C','D','sphericity','aplanarity']:
                    summarygr,summarygr2s=evshapeFits,evshapeFits2s
                if v in ['C_2','D_2','sphericity_2','aplanarity_2']:
                    summarygr,summarygr2s=evshape_2Fits,evshape_2Fits2s
                npts=summarygr.GetN()
                summarygr.SetPoint(npts,p0,npts)
                summarygr.SetPointError(npts,pmax-p0,p0-pmin,0,0)
                summarygr2s.SetPoint(npts,p0,npts)
                if pmax2s and pmin2s:
                    summarygr2s.SetPointError(npts,pmax2s-p0,p0-pmin2s,0,0)

            if psSlice!='inc' : continue

            #temp method to add stuff to legends
            def addToLegend(leg,chi2Res,opt='lp'):
                pVal,p0Val,pmaxVal,pminVal,pmaxVal2s,pminVal2s,chi2Gr=chi2Res
                legTxt=chi2Gr.GetTitle()
                #legTxt='#splitline{%s}{#scale[0.7]{%3.3f '%(chi2Gr.GetTitle(),p0Val)
                #if pminVal and pmaxVal:
                #    legTxt+='[%3.3f,%3.3f]}}'%(pminVal,pmaxVal)
                #elif pminVal:
                #    legTxt+='[%3.3f,n/a]}}'%pminVal
                #elif pmaxVal:
                #    legTxt+='[n/a,%3.3f]}}'%pmaxVal
                #else:
                #    legTxt+=' n/a}}'
                leg.AddEntry(chi2Gr,legTxt,opt)


            #do different combinations/comparisons
            cSlicesList={'ptll':[('inc_ptll=awa','away'), ('inc_ptll=tow','toward'),('inc_ptll=tra','transverse')],
                         'nj':[('nj=0,1','N_{j}=0'),('nj=1,2','N_{j}=1'),('nj=2,999','N_{j}#geq2')],
                         'inc':[(v,VARTITLES[v]),('chavgpz','#bar{p}_{z}'),('aplanarity','Aplanarity'),('sphericity','Sphericity')]}
            for pfix in cSlicesList:

                if pfix=='ptll' and v in ['C','D','sphericity','aplanarity']: continue
                
                cSlices=cSlicesList[pfix]

                #sum up chi^2
                chi2scanSum={}
                x,y=ROOT.Double(0),ROOT.Double(0)
                for ic in xrange(0,len(cSlices)):
                    cSlice,cTitle=cSlices[ic]
                    
                    v2comp=cSlice if pfix=='inc' else v
                    if v2comp==v and pfix=='inc' : continue
                
                    if pfix=='inc': cSlice=psSlice    
                    if not (v2comp,cSlice) in paramScanResults[p] : continue

                    igr = paramScanResults[p][(v2comp,cSlice)][6]
                    for ip in xrange(0,gr.GetN()):
                        igr.GetPoint(ip,x,y)
                        xval=float(x)
                        if not xval in chi2scanSum: chi2scanSum[xval]=[0,0]
                        chi2scanSum[xval][0] += 1
                        chi2scanSum[xval][1] += float(y)
                chi2scanSum=[ (x,chi2scanSum[x][1]/chi2scanSum[x][0]) for x in chi2scanSum ]

                try:
                    combResult=performChisquareFitTo(chi2scanSum)
                    chi20_comb,p0_comb,pmax_comb,pmin_comb,pmax2s_comb,pmin2s_comb,gr_comb=combResult                
                except:
                    print 'Could not combine',pfix,v
                    continue
                
                c.Clear()                

#                gr_comb.Draw('apc')
#                gr_comb.SetTitle('inclusive')
#                gr_comb.GetXaxis().SetTitle(p)
#                gr_comb.GetYaxis().SetTitle('#chi^{2} / dof')
#                gr_comb.GetYaxis().SetRangeUser(0,10) #max(chi20*4,chi20+4))
#                gr_comb.GetXaxis().SetTitleSize(0.05)
#                gr_comb.GetYaxis().SetTitleSize(0.05)
#                gr_comb.GetXaxis().SetLabelSize(0.04)
#                gr_comb.GetYaxis().SetLabelSize(0.04)
#                gr_comb.SetLineColor(1)
#                gr_comb.SetMarkerColor(1)
#                gr_comb.SetLineWidth(3)
#                
#                gr.Draw('c')
#                ci=ROOT.TColor.GetColor('#889093')
#                gr.SetLineColor(ci)
#                gr.SetMarkerColor(ci)
#                gr.SetFillColor(0)
#                gr.SetLineWidth(2)
#                gr.SetMarkerStyle(1)
#                gr.SetTitle('inclusive')

                #gr.Draw('apc')
                gr.Draw('ac')
                gr.SetTitle('inclusive')
                gr.GetXaxis().SetTitle(p)
                gr.GetYaxis().SetTitle('#chi^{2}')
                #gr.GetYaxis().SetTitle('#chi^{2} / dof')
                gr.GetYaxis().SetRangeUser(0,50) #max(chi20*4,chi20+4))
                gr.GetXaxis().SetTitleSize(0.05)
                gr.GetYaxis().SetTitleSize(0.05)
                gr.GetXaxis().SetLabelSize(0.04)
                gr.GetYaxis().SetLabelSize(0.04)
                gr.SetLineColor(1)
                gr.SetMarkerColor(1)
                gr.SetLineWidth(3)
                

                resultLog=open('%s/chi2scans_%s_%s_%s.dat'%(opt.input,v,pname,pfix),'w')
                resultLog.write('%s %s %s %s\n'%(p,v,psSlice,' '.join( [str(x) for x in combResult[0:6]] )))
                resultLog.write('%s %s %s %s\n'%(p,v,psSlice,' '.join( [str(x) for x in paramScanResults[p][(v,psSlice)][0:6]] )))
                
                leg=ROOT.TLegend(0.16,0.8,0.45,0.54) if pfix=='inc' else ROOT.TLegend(0.16,0.76,0.45,0.32)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetTextFont(42)
                leg.SetTextSize(0.045)
                leg.AddEntry(gr,'inclusive','l')
                #leg.AddEntry(gr,'combination','lp')
                #leg.AddEntry(gr_comb,'combination','lp')
                #if pfix=='inc' : leg.AddEntry(gr,VARTITLES[v],'l')
                #else : addToLegend(leg,paramScanResults[p][(v,psSlice)])

                #draw comparisons
                compColors=['#889093','#92c5de','#f4a582','#ca0020']
                for ic in xrange(0,len(cSlices)):
                    cSlice,cTitle=cSlices[ic]

                    v2comp=cSlice if pfix=='inc' else v
                    if v2comp==v and pfix=='inc' : continue
                
                    if pfix=='inc': cSlice=psSlice    
                    if not (v2comp,cSlice) in paramScanResults[p] : continue

                    paramScanResults[p][(v2comp,cSlice)][6].SetTitle(cTitle)
                    ci=ROOT.TColor.GetColor(compColors[ic])
                    paramScanResults[p][(v2comp,cSlice)][6].SetLineColor(ci)
                    paramScanResults[p][(v2comp,cSlice)][6].SetMarkerColor(ci)
                    paramScanResults[p][(v2comp,cSlice)][6].SetLineWidth(2)
                    paramScanResults[p][(v2comp,cSlice)][6].SetMarkerStyle(24+ic)
                    paramScanResults[p][(v2comp,cSlice)][6].Draw('c')
                    if pfix=='inc' : leg.AddEntry(paramScanResults[p][(v2comp,cSlice)][6],cTitle,'l')
                    else : addToLegend(leg,paramScanResults[p][(v2comp,cSlice)],'l')
                    resultLog.write('%s %s %s %s\n'%(p,v2comp,cSlice,' '.join( [str(x) for x in paramScanResults[p][(v2comp,cSlice)][0:6]] )))
                
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
                tex.DrawLatex(0.67,0.96,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')

                c.RedrawAxis()
                c.Modified()
                c.Update()
                for ext in ['png','pdf','root']: c.SaveAs('%s/chi2scans_%s_%s_%s.%s'%(opt.input,v,pname,pfix,ext))
                resultLog.close()

        for gr,gr2s in [(evshapeFits,evshapeFits2s),
                        (evshape_2Fits,evshape_2Fits2s),
                        (fluxFits,fluxFits2s)]:
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

            gr2s.Draw('p')
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
            try:
                tex.DrawLatex(0.6,0.88,'#scale[0.8]{%s=%3.3f}'%(p,perc[1]))
                tex.DrawLatex(0.6,0.8,'#scale[0.8]{[%3.3f,%3.3f]}'%(perc[0],perc[2]))
            except:
                pass
            tex.DrawLatex(0.67,0.96,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')
            c.RedrawAxis()
            c.Modified()
            c.Update()
            for ext in ['png','pdf','root']: c.SaveAs('%s/chi2summary_%s_%s.%s'%(opt.input,gr.GetName(),pname,ext))

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
