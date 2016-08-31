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
        targetH = target_plotter.Get('%s/%s_%s'%(pname,pname,tWnom))

        try:
            nomH.SetDirectory(0)
            targetH.SetDirectory(0)
        except:
            print 'Match not found for',pname
            continue

        outDir = tW_syst_plotter.mkdir(pname)
        outDir.cd()
        for syst in tWsysts:
            systH=syst_plotter.Get('%s/%s_%s'%(pname,pname,syst)).Clone()
            systH.SetDirectory(outDir)
            systH.Divide(nomH)
            systH.Multiply(targetH)
            for xbin in xrange(1,targetH.GetNbinsX()):
                ynom,yvar=targetH.GetBinContent(xbin),systH.GetBinContent(xbin)
                if ynom==0 and yvar==0 : continue
                if ynom!=0 :
                    relVar=yvar/ynom
                    if ROOT.TMath.Abs(1-relVar)>=1: 
                        relVar=0.001 if relVar<1 else 2.0
                        print '[Warning]',yvar,'->',relVar*ynom,'|',ynom,'@ bin=',xbin,'for',pname
                    systH.SetBinContent(xbin,relVar*ynom)
                elif yvar!=0:
                    print '[Warning]',yvar,'->',0,'@ bin=',xbin,'for',pname
                    systH.SetBinContent(xbin,0)
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
