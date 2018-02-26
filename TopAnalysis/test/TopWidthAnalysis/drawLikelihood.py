import ROOT
import os
import re
import sys
import numpy as np

def buildLikelihoodScan(baseDir,tag='100pseudodata',testStat='LEP_TEV',title='TEV',marker=20,color=1):
    """opens combine outputs and reads off the likelihood scans"""

    print '[buildLikelihoodScan] with',tag

    #regex to parse hypothesis being tested
    p=re.compile('\d+vs\d+')

    #init graph
    gr=ROOT.TGraph()
    gr.SetName(testStat)
    gr.SetTitle(title)
    gr.SetMarkerStyle(marker)
    gr.SetFillStyle(0)
    gr.SetMarkerColor(color)
    gr.SetLineColor(color)

    #loop over sub-directories
    nllScan=[]
    for d in os.listdir(baseDir):
        d=os.path.join(baseDir,d)
        if not os.path.isdir(d): continue
        if not tag in d : continue

        main,alt=[float(x)*1.31/100. for x in p.findall(d)[0].split('vs')]
        try:
            test={0:None,1:None}
            for x in test:
                fIn=ROOT.TFile.Open(os.path.join(d,'testStat_scan%dn_%s.root'%(x,testStat)))
                limit=fIn.Get('limit')
                limit.GetEntry(0)
                test[x]=limit.limit
                fIn.Close()
            deltaq=test[1]-test[0]
            if abs(deltaq)>1e3 : continue
            nllScan.append( (alt,deltaq) )
        except:
            pass
    
    #find the minimum
    nllScan=sorted(nllScan, key=lambda nllPoint: nllPoint[1])
    nllMin=nllScan[0]
    print '\tMinimum @ ',nllMin,'refining'

    #build the graph
    for x,nll in nllScan:
        npts=gr.GetN()
        gr.SetPoint(npts,x, nll)
        gr.Sort()

    #refine minimum
    xstart=nllMin[0]
    for ix in xrange(0,200):
        xup=xstart+ix*0.01
        inll=gr.Eval(xup)
        if inll<nllMin[1] :
            nllMin=(xup,inll)

        xdn=xstart-ix*0.01
        inll=gr.Eval(xdn)
        if inll<nllMin[1] :
            nllMin=(xdn,inll)

    #find 68%CI
    ci=[(0,10),(0,10)]
    for ix in xrange(0,2000):
        xup=nllMin[0]+ix*0.001
        inll=gr.Eval(xup)
        delta=abs(abs(inll-nllMin[1])-1)
        if delta<ci[1][1] :
            ci[1]=(xup,delta)

        xdn=nllMin[0]-ix*0.001
        inll=gr.Eval(xdn)
        delta=abs(abs(inll-nllMin[1])-1)
        if delta<ci[0][1] :
            ci[0]=(xdn,delta)

    print '\tMinimum found @',nllMin
    print '\tCI [%3.3f,%3.3f]'%(ci[0][0],ci[1][0])

    #put in graph
    x,y=ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,gr.GetN()):
        gr.GetPoint(i,x,y)
        gr.SetPoint(i,x,float(y)-nllMin[1])


    #sort values and return
    return gr,(ci[0][0],nllMin[0],ci[1][0])

def showLikelihoods(allGr,tag):
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.06)
    c.SetRightMargin(0.03)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.1)

    leg=ROOT.TLegend(0.17,0.8,0.5,0.65)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.05)
    for i in xrange(0,len(allGr)):        
        allGr[i].Draw('apc' if i==0 else 'pc')
        allGr[i].GetXaxis().SetTitle('#Gamma_{top} [GeV]')
        allGr[i].GetYaxis().SetTitle('-2 #Delta ln L/L_{max}')
        allGr[i].GetYaxis().SetTitleSize(0.05)
        allGr[i].GetXaxis().SetTitleSize(0.05)
        allGr[i].GetYaxis().SetLabelSize(0.04)
        allGr[i].GetXaxis().SetLabelSize(0.04)
        allGr[i].GetXaxis().SetTitleOffset(0.9)
        allGr[i].GetYaxis().SetTitleOffset(1.4)
        allGr[i].GetYaxis().SetRangeUser(0,20)
        allGr[i].GetXaxis().SetRangeUser(0,5)
        leg.AddEntry(allGr[i],allGr[i].GetTitle(),'lp')
    leg.Draw()

    cl=ROOT.TLine()
    cl.SetLineStyle(2)
    cl.SetLineColor(ROOT.kGray+2)
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.025)
    txt.SetTextColor(ROOT.kGray+3)
    for delta,title in [(1.0,'68% CL'),(3.84,'95% CL')]:
        txt.DrawLatex(0.2,delta+0.25,title)
        cl.DrawLine(0,delta,5,delta)

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.05)
    label.DrawLatex(0.17,0.87,'#scale[1.2]{#bf{CMS} #it{preliminary}}')                                                                     
    label.DrawLatex(0.55,0.955,'35.6 fb^{-1} (#sqrt{s} = 13 TeV)')

    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('likelihood%s.%s'%(tag,ext))



#ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


wScan=[90,100,120,130,140,160,180,200,300]
wScan=[50,100,200,300,400]
baseDir=sys.argv[1]

fitResGr={}
for w in wScan:

    tag='%dpseudodata'%w

    allGr=[]
    for testStat, title, marker, color,shift in [ ('LEP','LEP',24,ROOT.kAzure+3,0.1),
                                                  #('TEV','TEV',25,ROOT.kGray,-0.1),
                                                  ('PL','PL',20,1,0) ]:
        gr,fitRes=buildLikelihoodScan(baseDir,tag=tag,testStat=testStat,title=title,marker=marker,color=color) 
        if gr.GetN()==0 : continue

        allGr.append(gr)

        if not testStat in fitResGr:
            fitResGr[testStat]=ROOT.TGraphAsymmErrors()
            fitResGr[testStat].SetTitle(title)
            fitResGr[testStat].SetMarkerStyle(marker)
            fitResGr[testStat].SetMarkerColor(color)
            fitResGr[testStat].SetLineColor(color)
        npts=fitResGr[testStat].GetN()
        fitResGr[testStat].SetPoint(npts,w*1.31/100.+shift,fitRes[1])
        unc=max(abs(fitRes[1]-fitRes[0]),abs(fitRes[2]-fitRes[1]))
        fitResGr[testStat].SetPointError(npts,0,0,unc,unc)

    showLikelihoods(allGr,tag)


c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.06)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.15)
c.SetBottomMargin(0.1)

leg=ROOT.TLegend(0.17,0.8,0.5,0.65)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.05)
drawOpt='ap'
frame=None
for testStat in fitResGr:
    fitResGr[testStat].Draw(drawOpt)
    drawOpt='p'
    frame=frame if frame else fitResGr[testStat]
    leg.AddEntry(fitResGr[testStat],fitResGr[testStat].GetTitle(),'p')

frame.GetXaxis().SetTitle('#Gamma_{top}(exp) [GeV]')
frame.GetYaxis().SetTitle('#Gamma_{top}(fit) [GeV]')
frame.GetYaxis().SetTitleSize(0.05)
frame.GetXaxis().SetTitleSize(0.05)
frame.GetYaxis().SetLabelSize(0.04)
frame.GetXaxis().SetLabelSize(0.04)
frame.GetXaxis().SetTitleOffset(0.9)
frame.GetYaxis().SetTitleOffset(1.4)
frame.GetYaxis().SetRangeUser(frame.GetXaxis().GetXmin(),frame.GetXaxis().GetXmax()*1.2)

leg.Draw()

label = ROOT.TLatex()
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.05)
label.DrawLatex(0.17,0.87,'#scale[1.2]{#bf{CMS} #it{preliminary}}')                                                                     
label.DrawLatex(0.55,0.955,'35.6 fb^{-1} (#sqrt{s} = 13 TeV)')


cl=ROOT.TLine()
cl.SetLineStyle(2)
cl.SetLineColor(ROOT.kGray+2)
cl.DrawLine(frame.GetXaxis().GetXmin(),frame.GetXaxis().GetXmin(),frame.GetXaxis().GetXmax(),frame.GetXaxis().GetXmax())


c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('linearity.%s'%(ext))

