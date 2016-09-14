#!/usr/bin/env python
import os,sys
import re
import ROOT
from optparse import OptionParser
from tdrStyle import setTDRStyle
import CMS_lumi
from array import array
from subprocess import Popen, PIPE, STDOUT
import pprint

"""
converts a graph to a spline
"""
def convertToSpline(gr,color,width,splineMin=0.3,splineMax=4,npoints=500):

    gr.Sort()
    tsp=ROOT.TMVA.TSpline2("%s_spline"%gr.GetName(),gr.Clone())

    # make graph out of spline
    tspGrX=ROOT.TVector(npoints)
    tspGrY=ROOT.TVector(npoints)
    max_x,max_y=splineMin,tsp.Eval(splineMin)
    for step in xrange(0,npoints) :
        tempX=step*(splineMax-splineMin)/npoints + splineMin
        tspGrX[step]=tempX
        tspGrY[step]=tsp.Eval(tempX)
        if tspGrY[step]<max_y: continue
        max_x,max_y = tempX,tspGrY[step]

    #determine the lower and upper limits at 95% and 99% CL
    cls={
        0.95:{'ul':(max_x,max_y),'ll':(max_x,max_y)},
        0.99:{'ul':(max_x,max_y),'ll':(max_x,max_y)}
        }
    for step in xrange(0,npoints) :
        bound='ll' if tspGrX[step]<max_x else 'ul'
        for cl in cls:
            if ROOT.TMath.Abs(tspGrY[step]-(1-cl))>ROOT.TMath.Abs(cls[cl][bound][1]-(1-cl)) : continue
            cls[cl][bound]=(tspGrX[step],tspGrY[step])


    tspSplGr=ROOT.TGraph(tspGrX,tspGrY)
    tspSplGr.SetLineColor(color)
    tspSplGr.SetLineWidth(width)
    tspSplGr.SetName("%s_splinegr"%gr.GetName())

    return tspSplGr,cls

"""
opens fit results file
"""
def buildGraph( dirList,var,ms,color,name,title,doGraph=True):

    #init graph
    gr=None
    if doGraph:
        gr=ROOT.TGraphAsymmErrors()
        gr.SetMarkerStyle(ms)
        gr.SetFillStyle(0)
        gr.SetLineColor(color)
        gr.SetMarkerColor(color)
        gr.SetName(name)
        gr.SetTitle(title)
    else:
        binWidth = [0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.4,3.6,3.7,3.9,4.1,4.3]
        gr=ROOT.TH2F(name,title+';#Gamma/#Gamma_{SM};x=N[#Gamma(alt)/#Gamma(null)];-2 ln (L/L_{max})',len(binWidth)-1,array('d',binWidth),50,0.,1.)

    #read fit results from file and fill graph
    for altHypo,indir in dirList:

        try:
            fIn=ROOT.TFile.Open(indir+'/higgsCombinex_scan_obs.MultiDimFit.mH172.5.root')
            limit=fIn.Get('limit')
            if doGraph:
                fitRes=fIn.Get('w').getSnapshot("MultiDimFit")
                rooVar=fitRes.find(var)
                val,err=rooVar.getVal(),rooVar.getError()

                llgr=ROOT.TGraphErrors()
                for i in xrange(0,limit.GetEntriesFast()):
                    limit.GetEntry(i)
                    llgr.SetPoint(llgr.GetN(),getattr(limit,var),2*limit.deltaNLL)
                    
                llgr.Fit('pol2','RQ+','same',val,val+1.5*err)
                pol2=llgr.GetFunction('pol2')
                xhi=pol2.GetX(1,val,val+1.5*err)
                llgr.Fit('pol2','RQ+','same',ROOT.TMath.Max(val-1.5*err,0.),val)
                pol2=llgr.GetFunction('pol2')
                xlo=pol2.GetX(1,ROOT.TMath.Max(val-1.5*err,0.),val)
                if xlo<1e-3 : xlo=val
                llgr.Delete()
                np=gr.GetN()
                gr.SetPoint(np,altHypo,val)
                gr.SetPointError(np,0,0,xlo,xhi)
            else:
                for i in xrange(0,limit.GetEntriesFast()):
                    limit.GetEntry(i)
                    gr.Fill(altHypo,getattr(limit,var),2*limit.deltaNLL)
            fIn.Close()
        except:
            print 'skipping',altHypo,'@',indir
    
    return gr

"""
"""
def buildQuantilesAndPlot(dirList,opt):
    
    #prepare lists
    qobs,y,eyu,eyl={},{},{},{}
                
    #fill lists
    for d in dirList:
        isData=True if not 'pseudodata' in d else False
        widths=re.findall(r'\d+[\.]?\d*', d)
        mainHypo,altHypo=float(widths[0]),float(widths[1])
        pseudodataHypo=None
        if not isData : pseudodataHypo=float(widths[2])

        key=(isData,mainHypo,pseudodataHypo)
        if not key in qobs:
            qobs[key]=[]
            y[key]={}
            eyu[key]={}
            eyl[key]={}
            for prepost in ['prefit','postfit']:
                y[key][prepost]={}
                eyu[key][prepost]={}
                eyl[key][prepost]={}
                for hyp in ['null','alt']:
                    y[key][prepost][hyp]=[]
                    eyu[key][prepost][hyp]={}
                    eyl[key][prepost][hyp]={}
                    for ci in [1,2,3]:
                        eyu[key][prepost][hyp][ci]=[]
                        eyl[key][prepost][hyp][ci]=[]

        #parse results in file
        for prepost in ['prefit','postfit','obs']:
            statsFileName="%s/%s/%sstats__%s_%s_%s.txt"%(opt.indir,d,prepost,('%3.1fw'%(altHypo)).replace('.','p'),'inclusive','mlb')           
            for line in open(statsFileName,"r"):
                if prepost!='obs':
                    if "nulquant" in line :
                        tline = map(float,line.split(";")[1:8]);
                        y[key][prepost]['null'].append( (altHypo,tline[3]) )
                        for ci in [1,2,3]:
                            eyl[key][prepost]['null'][ci].append( (altHypo,tline[3]-tline[3-ci]) )
                            eyu[key][prepost]['null'][ci].append( (altHypo,tline[3+ci]-tline[3]) )
                    elif "altquant" in line :
                        tline = map(float,line.split(";")[1:8]);
                        y[key][prepost]['alt'].append( (altHypo,tline[3]) )
                        for ci in [1,2,3]:
                            eyl[key][prepost]['alt'][ci].append( (altHypo,tline[3]-tline[3-ci]) )
                            eyu[key][prepost]['alt'][ci].append( (altHypo,tline[3+ci]-tline[3]) )
                else:
                    if 'qobs' in line :
                        qobs[key].append(  (altHypo,float(line.replace("qobs;",""))) ) 
    
    #draw
    for key in qobs:
        tag='data' if key[0] else 'pseudodata'
        if key[2] is not None : tag += '_%3.1f'%key[2]
        tag+= '_mainHypo_%3.1f'%key[1]
        title='Data' if key[0] else '#Gamma/#Gamma_{SM}='
        if key[2] is not None : title += '%3.1f'%key[2]

        drawQuantiles(y[key]['prefit'],eyu[key]['prefit'],eyl[key]['prefit'],None,tag+'_exp',title,opt)
        drawQuantiles(y[key]['postfit'],eyu[key]['postfit'],eyl[key]['postfit'],qobs[key],tag+'_obs',title,opt)


"""
"""
def buildCLsAndPlot(dirList,opt):
    
    #prepare lists
    statList={'prefit':{}, 'postfit':{}, 'obs':{} }
    for prepost in statList:
        statList[prepost]["Separation"]={}
        statList[prepost]["qobs"]={}
        statList[prepost]["$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$"]={}
        statList[prepost]["$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$"]={}
        statList[prepost]["CL$_s$^{\\rm exp.}$"]={}
        statList[prepost]["CL$_s$^{\\rm obs.}$"]={}
                
    #fill lists
    for d in dirList:
        isData=True if not 'pseudodata' in d else False
        widths=re.findall(r'\d+[\.]?\d*', d)
        mainHypo,altHypo=float(widths[0]),float(widths[1])
        pseudodataHypo=None
        if not isData : pseudodataHypo=float(widths[2])

        key=(isData,mainHypo,pseudodataHypo)
        if not key in statList['prefit']['Separation']:
            for prepost in statList:
                for val in statList[prepost]:
                    statList[prepost][val][key]=[]

        #parse results in file
        for prepost in ['prefit','postfit','obs']:
            statsFileName="%s/%s/%sstats__%s_%s_%s.txt"%(opt.indir,d,prepost,('%3.1fw'%(altHypo)).replace('.','p'),'inclusive','mlb')           
            if prepost=='obs' : print statsFileName,altHypo
            for line in open(statsFileName,"r"):
                if "separation" in line :
                    statList[prepost]["Separation"][key].append( (altHypo,line.split('#')[0]) )
                if "qobs" in line :
                    statList[prepost]["qobs"][key].append( float(line.split(';')[1]) )
                elif "null exceeded density" in line :
                    statList[prepost]["$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$"][key].append( (altHypo,line.split('#')[0]) )
                elif "alt exceeded density" in line :
                    statList[prepost]["$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$"][key].append( (altHypo,line.split('#')[0]) )
                elif "cls expected" in line :
                    if mainHypo==altHypo:
                        statList[prepost]["CL$_s$^{\\rm exp.}$"][key].append( (altHypo,'1') )
                    else:
                        statList[prepost]["CL$_s$^{\\rm exp.}$"][key].append( (altHypo,line.split('#')[0]) )
                elif "cls observed" in line :
                    if mainHypo==altHypo:
                        statList[prepost]["CL$_s$^{\\rm obs.}$"][key].append( (altHypo,'1 \\pm 0') )
                    else:
                        statList[prepost]["CL$_s$^{\\rm obs.}$"][key].append( (altHypo,line.split('#')[0]) )
        
        buildHypoTestDist(d,key,altHypo,statList['obs']['qobs'][key][-1],opt)

    #fill the CLs graphs
    clsGraphs={}
    for prepost, val in [('prefit',  'CL$_s$^{\\rm exp.}$'),
                         ('postfit', 'CL$_s$^{\\rm exp.}$'),
                         ('obs',     'CL$_s$^{\\rm obs.}$')]:
        for key in statList[prepost][val]:
            if not key in clsGraphs:
                clsGraphs[key]={}
                for iprepost in ['prefit','postfit','obs']:
                    clsGraphs[key][iprepost]=ROOT.TGraphErrors()
                    clsGraphs[key][iprepost].SetName(iprepost)
                    clsGraphs[key][iprepost].SetTitle(val)

            allhypoWidth=[]
            for res in statList[prepost][val][key]:
                hypoWidth=1.324*res[0]
                if hypoWidth in allhypoWidth : continue
                allhypoWidth.append(hypoWidth)
                cls,clsunc=float(res[1].split('\\pm')[0]),0
                try:
                    clsunc=float(res[1].split('\\pm')[1])
                except:
                    pass
                np=clsGraphs[key][prepost].GetN()
                clsGraphs[key][prepost].SetPoint(np,hypoWidth,cls)
                clsGraphs[key][prepost].SetPointError(np,0,clsunc)

    #show the CLs graphs
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    c.SetLogy()
    for key in clsGraphs:
        c.Clear()

        maxWidth=4
        if key[0]==False : maxWidth=4*key[2]
        frame=ROOT.TH1F('frame',';#Gamma [GeV];CL_{S}',1,0,maxWidth)
        frame.Draw()
        frame.GetYaxis().SetRangeUser(1e-3,1)
        frame.GetYaxis().SetTitleSize(0.045)
        frame.GetXaxis().SetTitleSize(0.045)
        frame.GetYaxis().SetLabelSize(0.035)
        frame.GetXaxis().SetLabelSize(0.035)

        prefitGrSpl,prefitcls=convertToSpline(clsGraphs[key]['prefit'],ROOT.kTeal-1,2)
        prefitGrSpl.Draw('c')

        postfitGrSpl,postfitcls=convertToSpline(clsGraphs[key]['postfit'],ROOT.kRed+1,2)
        postfitGrSpl.Draw('c')
 
        obsGrSpl,obscls=convertToSpline(clsGraphs[key]['obs'],1,1,0.3,maxWidth)
        obsGrSpl.Draw('c')  
        clsGraphs[key]['obs'].SetLineWidth(1)
        clsGraphs[key]['obs'].SetMarkerStyle(20)
        clsGraphs[key]['obs'].Draw('p')

        leg=ROOT.TLegend(0.58,0.85,0.95,0.7)
        leg.SetTextFont(42)
        leg.SetTextSize(0.028)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetFillColor(0)
        dataTitle='Observed'
        if key[0]==False: dataTitle='Observed (#Gamma/#Gamma_{SM}=%3.1f)'%key[2]
        leg.AddEntry(clsGraphs[key]['obs'], dataTitle, 'pl')
        leg.AddEntry(postfitGrSpl, 'Post-fit model (#mu profile)', 'l')
        leg.AddEntry(prefitGrSpl,  'Pre-fit model (#mu=1)',        'l')
        leg.Draw()
        
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.6,0.9,'#bf{CMS} #it{Preliminary}')
        txt.DrawLatex(0.72,0.98,'#scale[0.7]{12.9 fb^{-1} (13 TeV)}')

        l=ROOT.TLine()
        l.SetLineColor(ROOT.kRed)
        l.SetLineStyle(2)
        l.DrawLine(0.0,0.05,maxWidth,0.05)
        l.DrawLine(0.0,0.01,maxWidth,0.01)

        c.Modified()
        c.Update()

        tag='data' if key[0] else 'pseudodata'
        if key[2] is not None : tag += '_%3.1f'%key[2]
        tag+= '_mainHypo_%3.1f'%key[1]

        #save plot
        for ext in ['png','pdf']:
            c.SaveAs('%s/cls_%s.%s'%(opt.outdir,tag,ext))

        #output to file
        fOut=ROOT.TFile.Open('%s/cls_%s.root'%(opt.outdir,tag),'RECREATE')
        for prepost in clsGraphs[key]:
            clsGraphs[key][prepost].Write()
        postfitGrSpl.Write()
        prefitGrSpl.Write()
        obsGrSpl.Write()
        fOut.Close()

        #output limits to file
        with open('%s/cls_%s_limits.txt'%(opt.outdir,tag),'w') as fOut:
            for cl in prefitcls:
                fOut.write('-'*50+'\n')
                fOut.write('Limits at %d%% CL\n'%(cl*100))
                fOut.write('\t pre-fit  (\\mu=1)    : [%3.2f,%3.2f]\n'%(prefitcls[cl]['ll'][0],  prefitcls[cl]['ul'][0]))
                fOut.write('\t post-fit (\\mu prof.): [%3.2f,%3.2f]\n'%(postfitcls[cl]['ll'][0], postfitcls[cl]['ul'][0]))
                fOut.write('\t obs      (\\mu prof.): [%3.2f,%3.2f]\n'%(obscls[cl]['ll'][0],     obscls[cl]['ul'][0]))
            fOut.write('-'*50+'\n')


"""
"""
def buildHypoTestDist(d,key,altHypo,qobs,opt):

    histos={}
    fs={'prefit':1001,'postfit':3004}
    fc={
        'prefit':{False:ROOT.kBlue-9,True:ROOT.kOrange+1},
        'postfit':{False:ROOT.kBlue-7,True:ROOT.kOrange+7}
        }
    for prepost in ['prefit','postfit']:
        histos[prepost]={}
        histos[prepost][True]=ROOT.TH1F(prepost+'null',';-2 ln [ L(alt)/L(SM) ]; Toys;',200,-150,150)
        histos[prepost][True].SetDirectory(0)
        histos[prepost][True].SetFillStyle(fs[prepost])
        histos[prepost][True].SetFillColor(fc[prepost][True])
        histos[prepost][True].SetLineColor(fc[prepost][True])
        
        histos[prepost][False]=histos[prepost][True].Clone(prepost+'alt')
        histos[prepost][False].SetDirectory(0)
        histos[prepost][False].SetFillStyle(fs[prepost])
        histos[prepost][False].SetFillColor(fc[prepost][False])
        histos[prepost][False].SetLineColor(fc[prepost][False])

        #fill histos
        fIn=ROOT.TFile.Open('%s/%s/x_%s.qvals.root'%(opt.indir,d,prepost))
        q=fIn.Get('q')
        for i in xrange(0,q.GetEntries()):
            q.GetEntry(i)
            hyptype=True if q.type>0 else False
            histos[prepost][hyptype].Fill(-2*q.q)
        fIn.Close()
            
    #build the plots
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)  

    leg=ROOT.TLegend(0.15,0.88,0.65,0.75)
    leg.SetTextFont(42)
    leg.SetTextSize(0.023)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)


    histos['prefit'][True].Draw('hist')
    histos['prefit'][True].GetXaxis().SetTitleSize(0.04)
    histos['prefit'][True].GetYaxis().SetTitleSize(0.04)
    histos['prefit'][True].GetYaxis().SetRangeUser(0,histos['prefit'][True].GetMaximum()*1.4)
    histos['prefit'][True].GetYaxis().SetTitleOffset(1.1)
    leg.AddEntry( histos['prefit'][True],'Pre-fit model (SM)','f')

    histos['prefit'][False].Draw('histsame')
    leg.AddEntry( histos['prefit'][False],'Pre-fit model (alt.)','f')

    histos['postfit'][True].Draw('histsame')
    leg.AddEntry( histos['postfit'][True],'Post-fit model (SM)','f')

    histos['postfit'][False].Draw('histsame')
    leg.AddEntry( histos['postfit'][False],'Post-fit model (alt.)','f')

    leg.SetNColumns(2)
    leg.Draw()

    l=ROOT.TArrow(qobs,100,qobs,0,0.05,'>')
    l.SetLineWidth(2)
    l.SetLineColor(1)
    l.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.05)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.9,'#bf{CMS} #it{Preliminary}')
    txt.DrawLatex(0.7,0.98,'#scale[0.7]{12.9 fb^{-1} (13 TeV)}')
    txt.DrawLatex(0.7,0.9,'#scale[0.7]{Null: %3.1f #Gamma/#Gamma_{SM}}'%(key[1]))
    txt.DrawLatex(0.7,0.85,'#scale[0.7]{Alt: %3.1f #Gamma/#Gamma_{SM}}'%(altHypo))
    dataTitle='Observed'
    if key[0]==False: dataTitle='Observed (#Gamma/#Gamma_{SM}=%3.1f)'%key[2]
    txt.DrawLatex(0.65,0.79,'#rightarrow #scale[0.5]{%s}'%dataTitle)
    
    c.Modified()
    c.Update()
    
    tag='data' if key[0] else 'pseudodata'
    if key[2] is not None : tag += '_%3.1f'%key[2]
    tag+= '_mainHypo_%3.1f'%key[1]
    for ext in ['png','pdf']:
        c.SaveAs('%s/hypotest_%3.1fvs%3.1f_%s.%s'%(opt.outdir,key[1],altHypo,tag,ext))

"""
"""
def drawQuantiles(y,eyu,eyl,qobs,tag,title,opt):

    c=ROOT.TCanvas('c','c',800,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.11)
    c.SetGridy()

    #frame
    binWidth = [0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.4,3.6,3.7,3.9,4.1,4.3]
    frame=ROOT.TH2F('frame',';#Gamma/#Gamma_{SM};-2 ln [ L(alt)/L(SM) ]',len(binWidth)-1,array('d',binWidth),1,-160,160.)
    frame.Draw()
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)

    leg=ROOT.TLegend(0.15,0.85,0.5,0.65)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetNColumns(2)

    allGrs=[]
    colors=[ROOT.kOrange,   ROOT.kBlue-9,
            ROOT.kOrange+1, ROOT.kBlue-7, 
            ROOT.kOrange+7, ROOT.kBlue   ]
    for ci in [3,2,1]:
        for hyp,hyptitle in [('null','SM'),('alt','alt.')]:
            ngr=len(allGrs)
            allGrs.append( ROOT.TGraphAsymmErrors() )
            allGrs[-1].SetFillColor(colors[ngr])
            allGrs[-1].SetFillStyle(1001)
            allGrs[-1].SetLineColor(colors[ngr])
            allGrs[-1].SetName('gr%d'%ngr)
            allGrs[-1].SetTitle('%s, %d#sigma'%(hyptitle,ci))
            for i in xrange(0,len(eyl[hyp][ci])):
                altHypo,cen=y[hyp][i]
                u,l=eyu[hyp][ci][i][1],eyl[hyp][ci][i][1]
                bin=frame.GetXaxis().FindBin(altHypo)
                xwid=frame.GetXaxis().GetBinWidth(bin)
                xcen=frame.GetXaxis().GetBinCenter(bin)-0.25*xwid
                xl,xu=0.2*xwid,0.25*xwid
                if hyp=='alt': 
                    xcen += 0.5*xwid
                    xl,xu=0.25*xwid,0.2*xwid
                np=allGrs[-1].GetN()
                allGrs[-1].SetPoint(np,xcen,cen)
                allGrs[-1].SetPointError(np,xl,xu,l,u)
            allGrs[-1].Sort()
            allGrs[-1].Draw('2')
            leg.AddEntry(allGrs[-1],allGrs[-1].GetTitle(),'f')
    leg.Draw()

    if qobs is not None:
        allGrs.append( ROOT.TGraphAsymmErrors() )
        allGrs[-1].SetFillStyle(0)
        allGrs[-1].SetMarkerStyle(20)
        allGrs[-1].SetMarkerColor(1)
        allGrs[-1].SetLineColor(1)
        allGrs[-1].SetTitle(title)
        allGrs[-1].SetName('data')
        allGrs[-1].SetLineWidth(2)        
        for i in xrange(0,len(qobs)):
            altHypo,cen=qobs[i]
            bin=frame.GetXaxis().FindBin(altHypo)
            xwid=frame.GetXaxis().GetBinWidth(bin)
            xcen=frame.GetXaxis().GetBinCenter(bin)
            np=allGrs[-1].GetN()
            allGrs[-1].SetPoint(np,xcen,cen)
            allGrs[-1].SetPointError(np,0.5*xwid,0.5*xwid,0,0)
        allGrs[-1].Draw('p')
        leg2=ROOT.TLegend(0.5,0.85,0.60,0.8)
        leg2.SetFillStyle(0)
        leg2.SetTextFont(42)
        leg2.SetTextSize(0.04)
        leg2.SetBorderSize(0)
        leg2.AddEntry(allGrs[-1],title,'lp')
        leg2.Draw()


    #header
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.05)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.9,'#bf{CMS} #it{Preliminary}')
    txt.DrawLatex(0.8,0.98,'#scale[0.7]{12.9 fb^{-1} (13 TeV)}')

    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/quantiles_%s.%s'%(opt.outdir,tag,ext))


"""
compare prefit and postfit nuisances
"""
def compareNuisances(inDir,opt):
   
    colors=[1, ROOT.kOrange,  ROOT.kRed+1, ROOT.kMagenta-9, ROOT.kBlue-7]
    ires=0
    frame=None
    gr1s,gr2s=ROOT.TGraph(),ROOT.TGraph()
    postFitNuisGr={}
    dx=0.2
    fitResSummary=[]
    for title,inF in [('obs.',       'higgsCombinex_scan_obs.MultiDimFit.mH172.5'),
                      ('exp. (x=0)', 'higgsCombinex_scan_0_exp.MultiDimFit.mH172.5'),
                      ('exp. (x=1)', 'higgsCombinex_scan_1_exp.MultiDimFit.mH172.5')]:
        ires+=1
        inF=ROOT.TFile.Open('%s/%s.root'%(inDir,inF))
        fit_s=inF.Get('w').getSnapshot("MultiDimFit")

        varNames=[]
        iter = fit_s.createIterator()
        var = iter.Next()
        while var :   
            pname=var.GetName()
            if not 'CMS_th1' in pname and not '_In' in pname :
                if not pname in ['x','r']: 
                    varNames.append( (pname,var.getVal(),var.getError()) )
            var=iter.Next()
        npars=len(varNames)
        
        #init frames if not yet available
        if frame is None:
            frame=ROOT.TH2F('frame',';;N x #sigma_{pre-fit}',npars,0,npars,1,-3,3)
            frame.SetDirectory(0)

            gr1s=ROOT.TGraph()
            gr1s.SetName('gr1s')
            gr1s.SetMarkerStyle(1)
            gr1s.SetMarkerColor(19)
            gr1s.SetLineColor(19) 
            gr1s.SetFillStyle(1001)
            gr1s.SetFillColor(19) 
            gr1s.SetPoint(0,0,-1)
            gr1s.SetPoint(1,npars,-1)
            gr1s.SetPoint(2,npars,1)
            gr1s.SetPoint(3,0,1)
            gr1s.SetPoint(4,-0,-1)
            gr2s=gr1s.Clone('gr2s')
            gr2s.SetMarkerColor(18) 
            gr2s.SetLineColor(18)
            gr2s.SetFillStyle(1001)
            gr2s.SetFillColor(18) 
            gr2s.SetPoint(0,0,-2)
            gr2s.SetPoint(1,npars,-2)
            gr2s.SetPoint(2,npars,2)
            gr2s.SetPoint(3,0,2)
            gr2s.SetPoint(4,0,-2)

        #save post fit parameter values in a graph
        postFitNuisGr[title]=ROOT.TGraphErrors()
        postFitNuisGr[title].SetTitle(title)
        postFitNuisGr[title].SetMarkerStyle(19+ires)
        postFitNuisGr[title].SetMarkerColor(colors[ires-1])
        postFitNuisGr[title].SetLineColor(colors[ires-1])
        postFitNuisGr[title].SetLineWidth(2)
        postFitNuisGr[title].SetFillStyle(0)
        for i in xrange(0,len(varNames)):
            pname,val,unc = varNames[i]
            np=postFitNuisGr[title].GetN()
            postFitNuisGr[title].SetPoint(np,np+0.2+ires*dx,val)
            postFitNuisGr[title].SetPointError(np,0.,unc)
            if ires==1:
                frame.GetXaxis().SetBinLabel(np+1,'#color[%d]{%s}'%((np%2)*10+1,pname))

        inF.Close()

    #show nuisances
    c=ROOT.TCanvas('c','c',1500,500)
    c.SetLeftMargin(0.1)
    c.SetTopMargin(0.1)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.15)
    c.SetGridx(True)
    frame.Draw()
    #frame.GetXaxis().SetLabelOffset(+0.1)
    frame.GetYaxis().SetRangeUser(-3,3)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(0.8)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleSize(0.06)
    gr2s.Draw('f')
    gr1s.Draw('f')
    leg=ROOT.TLegend(0.6,0.92,0.8,0.99)
    leg.SetNColumns(len(postFitNuisGr))
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1001)
    leg.SetFillColor(0)
    for ftitle in postFitNuisGr:
        postFitNuisGr[ftitle].Draw('p')
        leg.AddEntry(postFitNuisGr[ftitle],postFitNuisGr[ftitle].GetTitle(),'p')
    leg.Draw()
    
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.025)
    txt.SetTextColor(ROOT.kGray+3)
    for delta,title in [(1.0,'-1#sigma'),(2,'+2#sigma'),(-1,'-1#sigma'),(-2,'-2#sigma')]:
        txt.DrawLatex(frame.GetXaxis().GetXmax()+0.5,delta-0.2,title)      

    #header
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.05)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{Preliminary}')
    txt.DrawLatex(0.3,0.95,'#scale[0.7]{12.9 fb^{-1} (13 TeV)}')

    c.RedrawAxis()
    c.Modified()
    c.Update()
    
    for ext in ['png','pdf']:
        c.SaveAs('%s/nuisances.%s'%(inDir,ext))


"""
"""
def doFitSummary(fitGraphs,opt):
    #build the plots
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.15)
    c.SetTopMargin(0.01)

    for mainHypo in fitGraphs:
        for v,vtit in [('x','Hypothesis fraction')]:

            dataGr=buildGraph(fitGraphs[mainHypo]['data'],v,20,ROOT.kAzure-3,'data','obs.',True)
            dataGr.SetLineWidth(2)

            for pseudoDataHypo in fitGraphs[mainHypo]['sim']:
                print pseudoDataHypo
                c.Clear()
                #c.SetLogz()
                ll2d=buildGraph(fitGraphs[mainHypo]['sim'][pseudoDataHypo],
                                v,
                                1,
                                ROOT.kAzure,
                                'll2d',
                                'exp. #Gamma/#Gamma_{SM}=%3.1f'%pseudoDataHypo,
                                False)
                ll2d.SetFillStyle(1001)
                ll2d.SetFillColor(ROOT.kAzure)
                ll2d.Draw('colz')
                ll2d.GetXaxis().SetTitleSize(0.04)
                ll2d.GetYaxis().SetTitleSize(0.04)
                ll2d.GetZaxis().SetTitleSize(0.04)

                if dataGr.GetN()>0: dataGr.Draw('p')

                leg=ROOT.TLegend(0.45,0.98,0.83,0.78)
                leg.SetFillStyle(1001)
                leg.SetBorderSize(0)
                leg.SetTextFont(42)
                leg.SetTextSize(0.03)
                leg.SetHeader('#splitline{#bf{CMS} #it{preliminary}}{#scale[0.8]{12.9 fb^{-1} (13 TeV)}}')
                if dataGr.GetN()>0 : leg.AddEntry(dataGr,'obs.','p')
                leg.AddEntry(ll2d,'exp. #Gamma/#Gamma_{SM}=%3.1f'%pseudoDataHypo,'f')
                leg.Draw()

                c.Modified()
                c.Update()
                for ext in ['png','pdf','C']:
                    c.SaveAs('%s/%s_fitresults_%3.1f_%3.1f.%s'%(opt.outdir,v,pseudoDataHypo,mainHypo,ext))
                    
"""
steer the script
"""
def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(54)

    #parse user input
    parser = OptionParser(
        usage="%prog [options] [label=datacard.txt | datacard.txt]",
        epilog="Summarizes the results of the likelihood fit"
        )
    parser.add_option("-i",    type="string", dest="indir"  ,            default="datacards/", help="directory to look for stats files in")
    parser.add_option("-o",    type="string", dest="outdir"  ,           default=None,         help="directory save summary files")
    parser.add_option("--doFitSummary",       dest="doFitSummary",       default=False,        action="store_true", help="do fit summary")
    parser.add_option("--doNuisances",        dest="doNuisances"  ,      default=False,        action="store_true", help="compare pre-fit and post-fit nuisances")
    parser.add_option("--doCLs",              dest="doCLs"  ,            default=False,        action="store_true", help="do the CLs summary")
    parser.add_option("--recreateCLsSummary", dest="recreateCLsSummary", default=False,        action="store_true", help="recreate CLs summary")
    (opt, args) = parser.parse_args()

    if opt.outdir is None: opt.outdir=opt.indir

    #create a map of graphs to make
    fitGraphs={}
    dirs = [d for d in os.listdir(opt.indir) if os.path.isdir(os.path.join(opt.indir,d))]
    for d in dirs:
        if not 'hypotest' in d: continue
        isData=True if not 'pseudodata' in d else False
        widths=re.findall(r'\d+[\.]?\d*', d)
        mainHypo,altHypo=float(widths[0]),float(widths[1])
        pseudodataHypo=None
        if not isData : pseudodataHypo=float(widths[2])

        print 'Analysing',d,(isData,pseudodataHypo,mainHypo,altHypo)

        #store fit summary
        if opt.doFitSummary:
            if not mainHypo in fitGraphs:
                fitGraphs[mainHypo]={'data':[],'sim':{}}
            if isData:
                fitGraphs[mainHypo]['data'].append( (altHypo,os.path.join(opt.indir,d)) )
            else:
                if not pseudodataHypo in fitGraphs[mainHypo]['sim']: fitGraphs[mainHypo]['sim'][pseudodataHypo]=[]
                fitGraphs[mainHypo]['sim'][pseudodataHypo].append( (altHypo,os.path.join(opt.indir,d)) )

        #nuisance summary
        if opt.doNuisances: compareNuisances(os.path.join(opt.indir,d),opt)

        #summarize CLs
        if not opt.recreateCLsSummary : continue
        clsFiles=[f for f in os.listdir(os.path.join(opt.indir,d)) if 'higgsCombinecls' in f]
        for f in clsFiles:
            f='%s/%s/%s'%(opt.indir,d,f)
            prepost='prefit'
            if 'postfit_exp' in f : prepost='postfit'
            if 'postfit_obs' in f : prepost='obs'
            summaryF='%s/%s/x_%s.qvals.root'%(opt.indir,d,prepost)
            runClsSummary=Popen(['root', '-l', '-q', '-b', 
                                 f, 
                                 '%s/hypoTestResultTreeTopWid.cxx(\"%s\",172.5,1,\"x\",1000,\"inclusive\",\"%s\",\"mlb\",true,\"%s/%s/%s\")'%(os.path.dirname(os.path.realpath(__file__)),
                                                                                                                                 summaryF,
                                                                                                                                 ('%3.1fw'%(altHypo)).replace('.','p'),
                                                                                                                                 opt.indir,d,prepost)
                                 ])
            runClsSummary.communicate()

    if opt.doFitSummary: doFitSummary(fitGraphs,opt)
    if opt.doCLs: 
        buildQuantilesAndPlot(dirs,opt)
        buildCLsAndPlot(dirs,opt)
    



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
