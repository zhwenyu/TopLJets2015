import sys
import ROOT


def main():

    #open files
    plotter        = ROOT.TFile.Open(sys.argv[1])
    syst_plotter   = ROOT.TFile.Open(sys.argv[2])
    target_plotter = ROOT.TFile.Open(sys.argv[3])
    tW_syst_plotter = ROOT.TFile.Open('tW_syst_plotter.root','RECREATE')

    #base configuration expected from the json files
    tWnom='tW'
    tWsysts=["tW m=169.5","tW m=175.5",'tW scale up','tW scale down','tW DS']

    #build the histograms containing the relative variations
    plots = plotter.GetListOfKeys()
    for p in plots:
        pname   = p.GetName()
        nomH    = plotter.Get('%s/%s_%s'%(pname,pname,tWnom))
        nomH.SetDirectory(0)
        targetH = target_plotter.Get('%s/%s_%s'%(pname,pname,tWnom))
        targetH.SetDirectory(0)

        outDir = tW_syst_plotter.mkdir(pname)
        outDir.cd()
        for syst in tWsysts:
            systH=syst_plotter.Get('%s/%s_%s'%(pname,pname,syst)).Clone()
            systH.SetDirectory(outDir)
            systH.Divide(nomH)
            systH.Multiply(targetH)
            systH.Write()
        print pname
        tW_syst_plotter.cd()

    plotter.Close()
    syst_plotter.Close()
    target_plotter.Close()
    tW_syst_plotter.Close()





"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
