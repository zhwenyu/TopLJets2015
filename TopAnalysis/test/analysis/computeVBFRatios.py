import optparse
import ROOT
import sys
from TopLJets2015.TopAnalysis.Plot import *

def normalize(h):
    if h.Integral()<=0 : return
    h.Scale(1./h.Integral(0,h.GetNbinsX()+1))

def getPlotsIn(inF,dirname):
    """gets the data and mc sums in a given directory"""
    data,mc=None,None
    for key in inF.Get(dirname).GetListOfKeys():
        name=key.GetName()
        if 'Graph' in name : continue
        if name==dirname:
            data=key.ReadObj()
            data.SetDirectory(0)
        else:
            h=key.ReadObj()
            if not mc:
                mc=h.Clone(dirname+'_mc')
                mc.SetDirectory(0)
                mc.Reset('ICE')
            mc.Add(h)
    return data,mc

def computeVBFRatios(inUrl,triggerBased=False):
    """opens a plotter and computes the gamma/Z ratios for data and MC
       for different distributions"""

    ratios={}

    inF=ROOT.TFile.Open(inUrl)
    for key in inF.GetListOfKeys():
        name=key.GetName()

        doRatioTo=None
        ratioTitle=''
        if not triggerBased:
            if 'HighPtA_' in name or 'VBFA_' in name: 
                doRatioTo=('A_','MM_')
                ratioTitle='#gamma / Z#rightarrow #mu#mu ratio'
        else:
            if 'HighPtVBFA_' in name:
                doRatioTo=('HighPtVBFA_','HighPtA_')
                ratioTitle='High p_{T} VBF #gamma / High p_{T} #gamma'

        if doRatioTo:
            dataNum,mcNum=getPlotsIn(inF,name)
            dataDen,mcDen=getPlotsIn(inF,name.replace(doRatioTo[0],doRatioTo[1]))

            if triggerBased:
                for h in [dataNum,mcNum,dataDen,mcDen] : normalize(h)

            dataNum.Divide(dataDen)
            mcNum.Divide(mcNum)
            ratios[name]=(dataNum,mcNum)
            for i in xrange(0,2):
                ratios[name][i].GetYaxis().SetTitle(ratioTitle)
    inF.Close()

    return ratios

def showRatios(ratios,outUrl):
    """shows the ratio plots"""

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    os.system('rm %s'%outUrl)
    outDir=os.path.dirname(outUrl)

    for key in ratios:
        p=Plot(key+'_ratio',com='13 TeV')
        p.doPoissonErrorBars=False
        p.ratiorange=(0.38,1.64)
        p.add(h=ratios[key][0],title='Data',color=1,isData=True,spImpose=False,isSyst=False)
        #p.add(h=ratios[key][1],title='MC',color='#e5f5f9',isData=False,spImpose=False,isSyst=False)
        p.show(outDir=outDir,lumi=41400,noStack=False,saveTeX=False,noRatio=True)
        p.appendTo(outUrl)
        p.reset()

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',  dest='input',  help='input plotter [%default]',  default=None,                    type='string')
    parser.add_option('-o',  dest='output', help='output plotter [%default]', default='ratio_plotter.root',    type='string')
    parser.add_option('-t',  dest='triggerBased', help='trigger-based ratio [%default]', default=False, action='store_true')
    (opt, args) = parser.parse_args()

    ratios=computeVBFRatios(opt.input,opt.triggerBased)
    showRatios(ratios,opt.output)

if __name__ == "__main__":
    sys.exit(main())
