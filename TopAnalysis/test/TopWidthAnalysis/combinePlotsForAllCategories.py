from TopLJets2015.TopAnalysis.Plot import *
import sys

"""
"""
def doPlot(plotName,chList,url):

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    #open the file
    inF = ROOT.TFile.Open(url)

    plotsPerProc={}
    for ch in chList :
        pName='%s_%s'%(ch,plotName)
        pdir=inF.GetDirectory(pName)
        for key in pdir.GetListOfKeys():
            keyName=key.GetName()
            if 'Graph' in keyName : continue

            h=pdir.Get(keyName)
            title=keyName.replace(pName+'_','')
            if title==keyName : title='Data'
            if not title in plotsPerProc:
                plotsPerProc[title]=h.Clone(title)
                plotsPerProc[title].Reset('ICE')
                plotsPerProc[title].SetTitle(title)
                plotsPerProc[title].SetDirectory(0)
                plotsPerProc[title].SetFillColor(h.GetFillColor())                        
                plotsPerProc[title].SetLineColor(h.GetLineColor())                        
                plotsPerProc[title].SetMarkerColor(h.GetMarkerColor())   
            plotsPerProc[title].Add(h)
                    
    #show
    plot=Plot('%s%s'%(ch,plotName))    
    plot.savelog=True
    plot.wideCanvas=False
    plot.ratiorange=(0.76,1.24)
    plot.plotformats=['root','pdf','png']
    for key in  plotsPerProc:
        isData=True if 'Data' in plotsPerProc[key].GetTitle() else False
        color=1 if isData else plotsPerProc[key].GetFillColor()
        plot.add(plotsPerProc[key],
                 plotsPerProc[key].GetTitle(),
                 color,
                 isData,
                 False)
    plot.finalize()
    plot.mcUnc=0.062
    plot.show(outDir="plots/",lumi=12900)

                     
def main():

    os.system('mkdir -p plots')

    plots=sys.argv[1].split(',')
    chList='EE1b,EE2b,MM1b,MM2b,EM1b,EM2b'.split(',')
    if len(sys.argv)>2: chList=sys.argv[2].split(',')
    plotter='/afs/cern.ch/work/e/ecoleman/public/TopWidth/TopWidth_era2016/analysis/plots/plotter.root'
    if len(sys.argv)>3: plotter=sys.argv[3]
    for p in plots : doPlot(p,chList,plotter)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
