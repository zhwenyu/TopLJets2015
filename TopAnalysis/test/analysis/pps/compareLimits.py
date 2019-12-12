import ROOT
import sys
from TopLJets2015.TopAnalysis.Plot import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

colors=[ROOT.kBlack,ROOT.kGreen+1,ROOT.kAzure+1,ROOT.kRed+1,ROOT.kGray]
p=Plot('limitcomparison')
p.spimposeOpt='l'
p.xtit='m_{X} [GeV]'
p.ytit='95% CL limits on #sigma_{fid}#timesBR [pb]'
for i in range(len(sys.argv)-1):
    name,url=sys.argv[i+1].split(':')
    fIn=ROOT.TFile.Open(url)
    c=fIn.Get('c')
    gr=c.GetListOfPrimitives().At(2).GetListOfGraphs().At(2)
    gr.SetLineWidth(3)
    gr.SetLineStyle(1+i%2)
    p.add(gr,
          title=name, 
          color=colors[i], 
          isData=False,
          spImpose=True,
          isSyst=False)
    fIn.Close()
p.show('./',lumi=37500,noRatio=True)
p.reset()


p=Plot('limitcomparison_obs')
p.spimposeOpt='l'
p.xtit='m_{X} [GeV]'
p.ytit='95% CL limits on #sigma_{fid}#timesBR [pb]'
for i in range(len(sys.argv)-1):
    name,url=sys.argv[i+1].split(':')
    fIn=ROOT.TFile.Open(url)
    c=fIn.Get('c')
    gr=c.GetListOfPrimitives().At(2).GetListOfGraphs().At(3)
    gr.SetLineWidth(3)
    gr.SetLineStyle(1+i%2)
    p.add(gr,
          title=name, 
          color=colors[i], 
          isData=False,
          spImpose=True,
          isSyst=False)
    fIn.Close()
p.show('./',lumi=37500,noRatio=True)
p.reset()



