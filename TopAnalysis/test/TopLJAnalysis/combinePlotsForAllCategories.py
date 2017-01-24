from TopLJets2015.TopAnalysis.Plot import *
import sys
import json

"""
"""
def doPlot(plotName,ch,mcsampleOrder):

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    #open the files
    inFiles=[]
    baseDir='~/work/LJets-5TeV'
    if ch=='all' or ch=='mu':  inFiles.append(ROOT.TFile.Open('%s/analysis_mu/plots/plotter.root'%baseDir))
    if ch=='all' or ch=='e':   inFiles.append(ROOT.TFile.Open('%s/analysis_e/plots/plotter.root'%baseDir))

    plotsPerProc={}
    for f in inFiles:
        pName=plotName
        pdir=f.GetDirectory(pName)        
        for key in pdir.GetListOfKeys():
            keyName=key.GetName()
            if 'Graph' in keyName : continue

            h=pdir.Get(keyName)
            title=keyName.replace(pName+'_','')
            if title==keyName : title='Data'
            if 'Multijets' in title:
                for xbin in xrange(1,h.GetNbinsX()+1):
                    val,unc=h.GetBinContent(xbin),h.GetBinError(xbin)
                    h.SetBinError(xbin,ROOT.TMath.Sqrt(val*val+unc*unc))
                h.SetTitle('Multijets')
            if not title in plotsPerProc:
                plotsPerProc[title]=h.Clone(title)
                plotsPerProc[title].Reset('ICE')
                plotsPerProc[title].SetTitle(title)
                plotsPerProc[title].SetDirectory(0)
                plotsPerProc[title].SetFillColor(h.GetFillColor())                        
                plotsPerProc[title].SetLineColor(h.GetLineColor())                        
                plotsPerProc[title].SetMarkerColor(h.GetMarkerColor())   
                if 'drjj' in pName : plotsPerProc[title].GetXaxis().SetTitle("min#DeltaR(j,j')")                
            plotsPerProc[title].Add(h)
                    
    #show
    plot=Plot('%s%s'%(ch,plotName))    
    plot.savelog=True
    plot.wideCanvas=True if plotName=='nbtags' else False
    plot.ratiorange=(0.76,1.24)
    plot.plotformats=['root','pdf','png']
    for key in mcsampleOrder:
        if key in plotsPerProc:
            plot.add(plotsPerProc[key],
                     plotsPerProc[key].GetTitle(),
                     plotsPerProc[key].GetFillColor(),
                     False,
                     False)
    for key in  plotsPerProc:
        if not 'Data' in  plotsPerProc[key].GetTitle() : continue
        plot.add(plotsPerProc[key],
                 plotsPerProc[key].GetTitle(),
                 1,
                 True,
                 False)
    
    plot.finalize()
    plot.ratiorange=(0.32,1.82)
    plot.cmsLabel='#bf{CMS}'
    plot.com='5.02 TeV'
    plot.show(outDir="plots/",lumi=27.4,noRatio=True)
    #raw_input()
                     
def main():

    os.system('mkdir -p plots')

    samplesList=None
    with open(sys.argv[1],'r') as jsonFile :
        samplesList=json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
    mcsampleOrder=[]
    for tag,sample in samplesList:
        isData=sample[1]
        if isData==1: continue
        tag=sample[3]
        if tag in mcsampleOrder: continue
        mcsampleOrder.append( sample[3] )

    plots=sys.argv[2].split(',')
    ch='all'
    if len(sys.argv)>3: ch=sys.argv[3]
    for p in plots : doPlot(p,ch,mcsampleOrder)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
