import ROOT
import os
import re

valTitles={
    'r':'#mu = #sigma_{obs}(t#bar{t})/#sigma_{th}(t#bar{t})',
    'x':'#Gamma_{top}^{alt} fraction'
    }
vals={'r':[],'x':[]}
p=re.compile('\d+vs\d+')

testStr='_100pseudodata'
testStr='_400pseudodata'
testStr='_data'

for d in os.listdir('./'):
    if not os.path.isdir(d): continue
    if not testStr in d : continue


    main,alt=[float(x)*1.31/100. for x in p.findall(d)[0].split('vs')]

    for key in vals:
        try:
            fIn=ROOT.TFile.Open(os.path.join(d,'mlfit%s_fit_obs.root'%key))
            fit_s=fIn.Get('fit_s')
            var=fit_s.floatParsFinal().find(key)
            vals[key].append( (alt,
                               var.getVal(),
                               abs(var.getErrorHi()),
                               abs(var.getErrorLo()),
                               True if main==alt else False) )
            fIn.Close()
        except:
            pass

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.06)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.15)
c.SetBottomMargin(0.1)

for key in vals:
    gr=ROOT.TGraphAsymmErrors()
    gr.SetMarkerStyle(24)
    grmain=ROOT.TGraphAsymmErrors()
    grmain.SetMarkerStyle(20)
    for alt,v,verrHi,verrLo,ismain in vals[key]:
        print alt,v,verrHi,verrLo,ismain
        np=gr.GetN()
        gr.SetPoint(np,alt,v)
        gr.SetPointError(np,0,0,verrLo,verrHi)
        if not ismain: continue
        grmain.SetPoint(0,alt,v)
        grmain.SetPointError(0,0,0,verrLo,verrHi)

    gr.Sort()
    gr.Draw('ap')
    gr.GetXaxis().SetTitle('#Gamma_{top} [GeV]')
    gr.GetYaxis().SetTitle(valTitles[key])
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetLabelSize(0.04)
    gr.GetXaxis().SetLabelSize(0.04)
    gr.GetXaxis().SetTitleOffset(0.9)
    gr.GetYaxis().SetTitleOffset(1.4)
    if key=='x' : gr.GetYaxis().SetRangeUser(0,1.2)
    if key=='r' : gr.GetYaxis().SetRangeUser(0.8,1.2)
    grmain.Sort()
    grmain.Draw('p')

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.05)
    label.DrawLatex(0.17,0.87,'#scale[1.2]{#bf{CMS} #it{preliminary}}')                                                                     
    label.DrawLatex(0.55,0.955,'35.6 fb^{-1} (#sqrt{s} = 13 TeV)')

    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s%s_stability.%s'%(key,testStr,ext))
    
