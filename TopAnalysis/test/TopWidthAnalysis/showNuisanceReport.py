#!/usr/bin/env python
import os,sys
import ROOT
from optparse import OptionParser

COLORS=[1, ROOT.kOrange,  ROOT.kRed+1, ROOT.kMagenta-9, ROOT.kBlue-7]

"""
compare prefit and postfit nuisances
"""
def compareNuisances(args,out):


    varList={}
    for i in xrange(0,len(args)):
        title,url=args[i].split('=')

        inF=ROOT.TFile.Open(url)
        fit_s=inF.Get('w').getSnapshot("MultiDimFit")

        varNames=[]
        iter = fit_s.createIterator()
        var = iter.Next()
        while var :   
            pname=var.GetName()
            if not 'CMS_th1' in pname and not '_In' in pname :
                if not pname in ['x','r'] :
                    if not pname in varList: 
                        varList[pname]={}
                    varList[pname][title]=(var.getVal(),var.getErrorLo(),var.getErrorHi())
            var=iter.Next()
    
    keyCtr=0
    keyLists=[]
    for key in sorted(varList):
        if keyCtr%22==0: keyLists.append( [] )
        keyLists[-1].append(key)
        keyCtr+=1


    #show nuisances
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.3)
    c.SetTopMargin(0.1)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(1.0)
    c.SetGridy(True)
    for ilist in xrange(0,len(keyLists)):

        keyList=keyLists[ilist]
        npars=len(keyList)

        c.Clear()

        #prepare a frame
        frame=ROOT.TH2F('frame',';#hat{#theta};Nuisance parameter',1,-3,3,npars,0,npars)
        frame.SetDirectory(0)
        for ybin in xrange(0,len(keyList)): 
            frame.GetYaxis().SetBinLabel(ybin+1,'#color[%d]{%s}'%((ybin%2)*10+1,keyList[ybin]))
        frame.Draw()
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleOffset(0.8)
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().SetTitleOffset(2.4)

        #1-sigma
        gr1s=ROOT.TGraph()
        gr1s.SetName('gr1s')
        gr1s.SetMarkerStyle(1)
        gr1s.SetMarkerColor(19)
        gr1s.SetLineColor(19) 
        gr1s.SetFillStyle(1001)
        gr1s.SetFillColor(19) 
        gr1s.SetPoint(0,-1,0)
        gr1s.SetPoint(1,-1,npars)
        gr1s.SetPoint(2,1,npars)
        gr1s.SetPoint(3,1,0)
        gr1s.SetPoint(4,-1,0)

        #2-sigma
        gr2s=gr1s.Clone('gr2s')
        gr2s.SetMarkerColor(18) 
        gr2s.SetLineColor(18)
        gr2s.SetFillStyle(1001)
        gr2s.SetFillColor(18) 
        gr2s.SetPoint(0,-2,0)
        gr2s.SetPoint(1,-2,npars)
        gr2s.SetPoint(2,2,npars)
        gr2s.SetPoint(3,2,0)
        gr2s.SetPoint(4,-2,0)

        gr2s.Draw('f')
        gr1s.Draw('f')    
    
        txt=ROOT.TLatex()
        txt.SetTextFont(42)
        txt.SetTextSize(0.03)
        txt.SetTextColor(ROOT.kGray+3)
        for delta,title in [(1.0,'-1#sigma'),(2,'+2#sigma'),(-1,'-1#sigma'),(-2,'-2#sigma')]:
            txt.DrawLatex(delta-0.2,frame.GetYaxis().GetXmax()+0.2,title)      

        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.12,0.955,'#bf{CMS} #it{Preliminary}')
        txt.DrawLatex(0.72,0.955,'#scale[0.7]{35.6 fb^{-1} (13 TeV)}')

        nuisGrs={}
        for key in keyList:
            ires=0
            dy=1.0/float(len(varList[key])+1)
            for title in varList[key]:
                if not title in nuisGrs:
                    idx=len(nuisGrs)
                    nuisGrs[title]=ROOT.TGraphAsymmErrors()
                    nuisGrs[title].SetTitle(title)  
                    nuisGrs[title].SetMarkerStyle(20+idx)
                    nuisGrs[title].SetMarkerColor(COLORS[idx])
                    nuisGrs[title].SetLineColor(COLORS[idx])
                    nuisGrs[title].SetLineWidth(2)
                    nuisGrs[title].SetFillStyle(0)    
                np=nuisGrs[title].GetN()
                val,uncLo,uncHi = varList[key][title]
                nuisGrs[title].SetPoint(np,val,np+dy*(ires+1))
                nuisGrs[title].SetPointError(np,abs(uncLo),abs(uncHi),0.,0.)  
                ires+=1

        leg=ROOT.TLegend(0.7,0.9,0.95,0.9-0.08*len(args))    
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for title in sorted(nuisGrs):
            nuisGrs[title].Draw('p')
            leg.AddEntry(nuisGrs[title],title,'p')
        leg.Draw()

        c.RedrawAxis()
        c.Modified()
        c.Update()
    
        for ext in ['png','pdf']:
            c.SaveAs('%s_%d.%s'%(out,ilist,ext))

    return


#    dy=0.5/float(len(args))
#    frame=None
#    gr1s,gr2s=ROOT.TGraph(),ROOT.TGraph()
#    postFitNuisGr={}
#    fitResSummary=[]
#    for i in xrange(0,len(args)):
#        title,url=args[i].split('=')
#
#        inF=ROOT.TFile.Open(url)
#        fit_s=inF.Get('w').getSnapshot("MultiDimFit")
#
#        varNames=[]
#        iter = fit_s.createIterator()
#        var = iter.Next()
#        while var :   
#            pname=var.GetName()
#            if not 'CMS_th1' in pname and not '_In' in pname :
#                if not pname in ['x','r']: 
#                    varNames.append( (pname,var.getVal(),var.getErrorLo(),var.getErrorHi()) )
#            var=iter.Next()
#        npars=len(varNames)
#        
#        #init frames if not yet available
#        if frame is None:
#
#
#        
#        #save post fit parameter values in a graph
#        postFitNuisGr[title]=ROOT.TGraphAsymmErrors()
#        postFitNuisGr[title].SetTitle(title)
#        postFitNuisGr[title].SetMarkerStyle(20+i)
#        postFitNuisGr[title].SetMarkerColor(COLORS[i])
#        postFitNuisGr[title].SetLineColor(COLORS[i])
#        postFitNuisGr[title].SetLineWidth(2)
#        postFitNuisGr[title].SetFillStyle(0)
#        for k in xrange(0,len(varNames)):
#            pname,val,uncLo,uncHi = varNames[k]
#            np=postFitNuisGr[title].GetN()
#            postFitNuisGr[title].SetPoint(np,val,np+(i+1)*dy)
#            postFitNuisGr[title].SetPointError(np,abs(uncLo),abs(uncHi),0.,0.)
#            if i==1:
#                frame.GetYaxis().SetBinLabel(np+1,'#color[%d]{%s}'%((np%2)*10+1,pname))
#        inF.Close()
#
#    frame.Draw()
#    #frame.GetXaxis().SetLabelOffset(+0.1)
#    frame.GetXaxis().SetRangeUser(-3,3)
#    frame.GetYaxis().SetLabelSize(0.05)
#    frame.GetYaxis().SetTitleOffset(0.8)
#    frame.GetXaxis().SetLabelSize(0.04)
#    frame.GetXaxis().SetTitleSize(0.06)
#    
#    leg=ROOT.TLegend(0.6,0.92,0.8,0.99)
#    leg.SetNColumns(len(postFitNuisGr))
#    leg.SetTextFont(42)
#    leg.SetTextSize(0.035)
#    leg.SetBorderSize(0)
#    leg.SetFillStyle(1001)
#    leg.SetFillColor(0)
#    for ftitle in postFitNuisGr:
#        postFitNuisGr[ftitle].Draw('p')
#        leg.AddEntry(postFitNuisGr[ftitle],postFitNuisGr[ftitle].GetTitle(),'p')
#    leg.Draw()
#    
#
                    
"""
steer the script
"""
def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #parse user input
    parser = OptionParser(
        usage="%prog -o nuisances label1=MultDimFit_1.root label2=MultiDimFit_2.root",
        epilog="Summarizes the post-fit results"
        )
    parser.add_option("-o",    type="string", dest="out"  ,          default='nuisances',  help="name of the output file")
    (opt, args) = parser.parse_args()
    
    compareNuisances(args,opt.out)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
