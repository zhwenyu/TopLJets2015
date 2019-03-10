from TopLJets2015.TopAnalysis.Plot import *

comp={"pu":["puup","pudn"],
      "sel":["selup","seldn"],
      "l1":["l1prefireup","l1prefiredn"],
      "qg":["gluonqg","quarkqg"],
      'jer':["JERup","JERdn"],
      'AbsoluteStatJEC':["AbsoluteStatJECup","AbsoluteStatJECdn"],
      'AbsoluteScaleJEC':["AbsoluteScaleJECup","AbsoluteScaleJECdn"],
      'AbsoluteMPFBiasJEC':["AbsoluteMPFBiasJECup","AbsoluteMPFBiasJECdn"],
      "FragmentationJEC":["FragmentationJECup","FragmentationJECdn"],
      "SinglePionECALJEC":["SinglePionECALJECup","SinglePionECALJECdn"],
      "SinglePionHCALJEC":["SinglePionHCALJECup","SinglePionHCALJECdn"],
      "FlavorPureGluonJEC":["FlavorPureGluonJECup","FlavorPureGluonJECdn"],
      "FlavorPureQuarkJEC":["FlavorPureQuarkJECup","FlavorPureQuarkJECdn"],
      "FlavorPureCharmJEC":["FlavorPureCharmJECup","FlavorPureCharmJECdn"],
      "FlavorPureBottomJEC":["FlavorPureBottomJECup","FlavorPureBottomJECdn"],
      "TimePtEtaJEC":["TimePtEtaJECup","TimePtEtaJECdn"],
      "RelativeJEREC1JEC":["RelativeJEREC1JECup","RelativeJEREC1JECdn"],
      "RelativeJEREC2JEC":["RelativeJEREC2JECup","RelativeJEREC2JECdn"],
      "RelativeJERHFJEC":["RelativeJERHFJECup","RelativeJERHFJECdn"],
      "RelativePtBBJEC":["RelativePtBBJECup","RelativePtBBJECdn"],
      "RelativePtEC1JEC":["RelativePtEC1JECup","RelativePtEC1JECdn"],
      "RelativePtEC2JEC":["RelativePtEC2JECup","RelativePtEC2JECdn"],
      "RelativePtHFJEC":["RelativePtHFJECup","RelativePtHFJECdn"],
      "RelativeBalJEC":["RelativeBalJECup","RelativeBalJECdn"],
      "RelativeFSRJEC":["RelativeFSRJECup","RelativeFSRJECdn"],
      "RelativeStatFSRJEC":["RelativeStatFSRJECup","RelativeStatFSRJECdn"],
      "RelativeStatECJEC":["RelativeStatECJECup","RelativeStatECJECdn"],
      "RelativeStatHFJEC":["RelativeStatHFJECup","RelativeStatHFJECdn"],
      "PileUpDataMCJEC":["PileUpDataMCJECup","PileUpDataMCJECdn"],
      "PileUpPtRefJEC":["PileUpPtRefJECup","PileUpPtRefJECdn"],
      "PileUpPtBBJEC":["PileUpPtBBJECup","PileUpPtBBJECdn"],
      "PileUpPtEC1JEC":["PileUpPtEC1JECup","PileUpPtEC1JECdn"],
      "PileUpPtEC2JEC":["PileUpPtEC2JECup","PileUpPtEC2JECdn"],
      "PileUpPtHFJEC":["PileUpPtHFJECup","PileUpPtHFJECdn"]
}

colors=[2,8]

def getProjections(h,c):
    proj=[]
    for ybin in range(h.GetNbinsY()):
        name=h.GetYaxis().GetBinLabel(ybin+1)
        if not name in c: continue
        proj.append( h.ProjectionX(name,ybin+1,ybin+1) )
        proj[-1].SetTitle(name)
        proj[-1].SetLineColor(colors[len(proj)-1])
        proj[-1].SetMarkerColor(colors[len(proj)-1])

    return proj

def drawVariations(h,exp,cat):

    for c in comp:
        projs=getProjections(exp,comp[c])

        p=Plot('%s_%s'%(cat,c))
        p.doPoissonErrorBars=False
        p.add(h,'nominal',1,True,False,False)
        for px in projs:
            if c=='qg':
                px.Scale(h.Integral()/px.Integral())
            p.add(px,px.GetTitle(),px.GetLineColor(),False,False,False)
        p.show('./',41000,True)


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

fIn=ROOT.TFile.Open(sys.argv[1])
for cat in ['HighVPtLowMJJA','LowVPtHighMJJA','HighVPtHighMJJA']:
    h=fIn.Get('%s_tagjetresol'%cat)
    exp=fIn.Get('%s_tagjetresol_exp'%cat)
    drawVariations(h,exp,cat)
