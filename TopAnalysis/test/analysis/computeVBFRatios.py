import optparse
import ROOT
import sys
from TopLJets2015.TopAnalysis.Plot import *

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

def computeVBFRatios(inUrl):
    """opens a plotter and computes the gamma/Z ratios for data and MC
       for different distributions"""

    ratios={}

    inF=ROOT.TFile.Open(inUrl)
    for key in inF.GetListOfKeys():
        name=key.GetName()
        if 'HighPtA' in name or 'VBFA' in name:
            dataA,mcA=getPlotsIn(inF,name)
            dataZ,mcZ=getPlotsIn(inF,name.replace('A_','MM_'))
            dataA.Divide(dataZ)
            mcA.Divide(mcZ)
            ratios[name]=(dataA,mcA)
            for i in xrange(0,2):
                ratios[name][i].GetYaxis().SetTitle('#gamma / Z#rightarrow #mu#mu ratio')
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
        p.add(h=ratios[key][1],title='MC',color='#e5f5f9',isData=False,spImpose=False,isSyst=False)
        p.show(outDir=outDir,lumi=41400,noStack=False,saveTeX=False)
        p.appendTo(outUrl)
        p.reset()

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',  dest='input',  help='input plotter [%default]',  default=None,                    type='string')
    parser.add_option('-o',  dest='output', help='output plotter [%default]', default='ratio_plotter.root',    type='string')
    (opt, args) = parser.parse_args()

    ratios=computeVBFRatios(opt.input)
    showRatios(ratios,opt.output)

if __name__ == "__main__":
    sys.exit(main())
