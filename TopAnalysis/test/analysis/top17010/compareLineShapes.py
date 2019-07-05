import ROOT
from TopLJets2015.TopAnalysis.Plot import *

simScenarios={    
    (169.5,1.23) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets_m1695.root',
    (172.5,1.31) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets.root',
    (175.5,1.39) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets_m1755.root',
    }
rwScenarios={
    (169.5,1.23) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/scenario131127/MC13TeV_2016_TTJets.root',
    (172.5,1.0)  : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/scenario917534/MC13TeV_2016_TTJets.root',
    (172.5,2.0)  : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/scenario917634/MC13TeV_2016_TTJets.root',
    (175.5,1.39) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/scenario1704006/MC13TeV_2016_TTJets.root',
    }


def getDist(url):
    fIn=ROOT.TFile.Open(url)
    h=fIn.Get('em_mlb').Clone()
    h.Add(fIn.Get('ee_mlb'))
    h.Add(fIn.Get('mm_mlb'))
    h.SetDirectory(0)
    h.Scale(1./h.Integral())
    h.GetYaxis().SetRangeUser(1e-4,0.005)
    fIn.Close()
    return h


simScenariosH={}
for key in simScenarios:
    simScenariosH[key]=getDist(simScenarios[key])
    simScenariosH[key].SetName('sim%d'%len(simScenariosH))
    simScenariosH[key].SetTitle('(%3.1f,%3.2f)'%key)

rwScenariosH={}
for key in rwScenarios:
    rwScenariosH[key]=getDist(rwScenarios[key])
    rwScenariosH[key].SetName('rw%d'%len(rwScenariosH))
    rwScenariosH[key].SetTitle('(%3.1f,%3.2f)'%key)

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

p=Plot('mlbwidth_closure','13 TeV')
p.ratiorange=[0.95,1.05]
p.savelog=True
p.doPoissonErrorBars=False
p.add(h=simScenariosH[(172.5,1.31)].Clone(),title=simScenariosH[(172.5,1.31)].GetTitle(),color=1,isData=True,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.add(h=rwScenariosH[(172.5,1.0)].Clone(), title=rwScenariosH[(172.5,1.0)].GetTitle(),color=2,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.add(h=rwScenariosH[(172.5,2.0)].Clone(), title=rwScenariosH[(172.5,2.0)].GetTitle(),color=3,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.show('./',lumi=35900,noStack=True)


for name,key in [('m1695',(169.5,1.23)),
                 ('m1755',(175.5,1.39))]:

    p=Plot('%s_closure'%name,'13 TeV')
    p.ratiorange=[0.95,1.05]
    p.savelog=True
    p.doPoissonErrorBars=False
    p.add(h=simScenariosH[key],title='Simulation',color=1,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
    p.add(h=rwScenariosH[key], title='Reweighted',color=2,isData=True,spImpose=False,isSyst=False,doDivideByBinWidth=True)
    p.show('./',lumi=35900,noStack=True)

