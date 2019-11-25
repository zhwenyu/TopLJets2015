import ROOT
from TopLJets2015.TopAnalysis.Plot import *

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
lumi=37500
xangle=int(sys.argv[1])
mass=int(sys.argv[2])
dilFactorBkg=0.069  #extrapolation factor: #(pt>40 && RPin && hpur) / #(pt>40)
dilFactorSig=0.192
baseDir="/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/"


#signal
inF = ROOT.TFile.Open(baseDir+"Chunks/Z_m_X_%d_xangle_%d_2017_preTS2_opt_v1_simu_reco.root"%(mass,xangle))
d=inF.Get("data")
d.Draw("mmiss >> hsigmix(50,0,3000)","wgt*(mixType==2 && cat==169 && mmiss>0 && xangle==%d && bosonpt>40)"%xangle,"goff")
hsigmix=ROOT.gDirectory.Get('hsigmix')
hsigmix.SetDirectory(0)
hsigmix.SetTitle('m_{X}=%d GeV'%mass)
d.Draw("mmiss >> hsig(50,0,3000)","wgt*(mixType==1 && cat==169 && mmiss>0 && xangle==%d && bosonpt>40)"%xangle,"goff")
hsig=ROOT.gDirectory.Get('hsig')
hsig.SetDirectory(0)
hsig.SetTitle('m_{X}=%d GeV'%mass)
inF.Close()

#main bkg
d=ROOT.TChain('data')
for f in [x for x in os.listdir(baseDir+"/Chunks") if 'DY50toInf' in x]:
    url=baseDir+"/Chunks/"+f
    d.AddFile(url)
d.Draw("mmiss >> hdy(50,0,3000)","wgt*(mixType==1 && cat==169 && mmiss>0 && xangle==%d && bosonpt>40)"%xangle,"goff")
hdy=ROOT.gDirectory.Get('hdy').Clone()
hdy.SetDirectory(0)
hdy.Scale(5765.4*lumi)
inF.Close()


#construct different scenarios
for f in [0.01,0.1,0.2,0.5]:

    p=Plot('mix%d_mmiss_dysig_a%d_f%d'%(mass,xangle,f*100),com='13 TeV')
    p.doPoissonErrorBars=False
    p.ratiorange=(0.6,1.37)
    p.spimposeWithErrors=True
    p.frameMax=1.3
    p.ratiotitle='#frac{DY+mix.signal}{DY}'
    fhsig=hsig.Clone('sig')
    fhsig.Scale(hdy.Integral()/hsig.Integral())
    fhsig.Scale(f)

    fhsigmix=hsigmix.Clone('sigmix')
    fhsigmix.Scale(hdy.Integral()/hsigmix.Integral())
    fhsigmix.Scale(f*dilFactorSig)

    totalH=hdy.Clone('totalobs')
    totalH.Scale(1-f)
    totalH.Add(fhsig)
    p.add(h=totalH, title='Pseudo-data', color=1,  isData=False,  spImpose=True, isSyst=False)
    
    totalBkgf0=hdy.Clone('totalbkgf0')
    p.add(h=totalBkgf0, title='DY', color=ROOT.kRed+1,  isData=False, spImpose=False, isSyst=False)
    
    totalBkg=hdy.Clone('totalbkg')
    totalBkg.Scale(1-f*dilFactorBkg)
    totalBkg.Add(fhsigmix)
    totalBkg.Scale(totalBkgf0.Integral()/totalBkg.Integral())
    p.add(h=totalBkg, title='DY+mix.signal', color=ROOT.kAzure+1,  isData=True, spImpose=False, isSyst=False)

    extraText='%d #murad\\m_{X}=%d GeV\\signal fraction: %3.2f\\'%(xangle,mass,f)
    p.show(outDir='./',lumi=lumi,noStack=True,saveTeX=False,noRatio=False,extraText=extraText)
    p.reset()


    
