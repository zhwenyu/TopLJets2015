import sys
import ROOT

from optparse import OptionParser
from TopLJets2015.TopAnalysis.Plot import *

def main():

    parser = OptionParser(
        usage="%prog [options] [label=datacard.txt | datacard.txt]",
        epilog=""
        )
    parser.add_option("-i",     type="string", dest="input"  ,   default="shapes.root",   help="file to look for dists in")
    parser.add_option("--cats", type="string", dest="cats",      default="highptEM2b",    help="a list of categories")
    parser.add_option("--altProc",  type="string",  dest="altProc",       default="alternative hypothesis",    help="overlay this alternative hypothesis")
    parser.add_option("--uncs", type="string", dest="uncList",   default="les_EM,jes,jer",   help="a list of uncertainties to look for")
    parser.add_option("--proc", type="string", dest="proc",      default="tbart1p0w",         help="a list of processes to plot")
    parser.add_option("-o",     type="string", dest="outdir" ,   default="./plots",       help="the base filename for the plots")
    parser.add_option("--obs",
                      type="string",
                      dest="obs",
                      default="incmlb",
                      help="observable to consider")
    (opt, args) = parser.parse_args()

    uncList=opt.uncList.split(',')
    systNameMap= {
        "jes"           : "jet E scale",
        "jer"           : "jet E resolution",
        "les_EE"        : "lepton E scale EE",
        "les_EM"        : "lepton E scale EM",
        "les_MM"        : "lepton E scale MM",
        "les"           : "lepton E scale",
        "btag"          : "b-tag efficiencies",
        "ltag"          : "light flavour mistag rate",
        "trigEE"        : "Trigger efficiencies",
        "trigEM"        : "Trigger efficiencies",
        "trigMM"        : "Trigger efficiencies",
        "sel_EE"        : "EE Selection efficiencies",
        "sel_EM"        : "EM Selection efficiencies",
        "sel_MM"        : "MM Selection efficiencies",
        "sel"           : "Selection efficiency",
        "toppt"         : "t#bar{t} top p_{T}",
        "pu"            : "Pileup",
        "ttPSScale"     : "t#bar{t} PS scale",
        "tWQCDScale"    : "tW QCD Scale",
        "tWttInterf"    : "tW/t#bar{t} interference",
        "Mtop"          : "m_{t}",
        "ttMEqcdscale"  : "t#bar{t} #mu_{R}/#mu_{F}",
        "ttPartonShower": "t#bar{t} hadronizer",
        "ttGenerator"   : "t#bar{t} ME generator",
        "ttPDF"         : "NNLOPDF3.0 variation",
        "cflip"         : "Charge flip",
        "ISR"           : "Initial state radiation",
        "FSR"           : "Final state radiation",
        "UE"            : "Underlying event",
        "noSC"          : "No SC",
        "hdamp"         : "HDAMP",
        "Herwig"        : "Herwig",
        "amcnlo"        : "aMC@NLO"
        }

    obsTitleMap={'incmlb':'Inclusive mass(l,b) [GeV]'}

    colors=[ROOT.kRed+1,ROOT.kBlue,
            ROOT.kOrange+8,ROOT.kGreen-4,
            ROOT.kMagenta+2,ROOT.kBlue-9,
            ROOT.kMagenta,ROOT.kGreen+2]

    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    fIn=ROOT.TFile.Open(opt.input)
    for cat in opt.cats.split(',') :

        plot=Plot('%s_%s_%s_%s'%(cat,opt.obs,uncList[0],opt.proc))
        plot.savelog=False
        plot.wideCanvas=False
        plot.frameMin=0.8
        #plot.frameMax=1.2
        plot.ratiorange=(0.76,1.24)
        plot.plotformats=['pdf','png']
        nomH=fIn.Get('%s_%s/%s'%(cat,opt.obs,opt.proc))
        altH=None
        if len(opt.altProc)!=0:
            try:
                altH=fIn.Get('%s_%s/%s'%(cat,opt.obs,opt.altProc))
                altH.Divide(nomH)
                altH.GetXaxis().SetTitle(obsTitleMap[opt.obs])
                altH.GetYaxis().SetTitle("Ratio to reference")
                plot.add(altH,"Alt.",1, False,True,False)
            except:
                altH=None

        isyst=-1
        for syst in opt.uncList.split(','):
            isyst+=1
            upH=fIn.Get('%s_%s_%sUp/%s'%(cat,opt.obs,syst,opt.proc))
            upH.Divide(nomH)
            upH.GetXaxis().SetTitle(obsTitleMap[opt.obs])
            upH.GetYaxis().SetTitle("Ratio to reference")
            plot.add(upH,systNameMap[syst]+" (up)",colors[isyst], False,False,False)
            isyst+=1
            dnH=fIn.Get('%s_%s_%sDown/%s'%(cat,opt.obs,syst,opt.proc))
            dnH.Divide(nomH)
            dnH.GetXaxis().SetTitle(obsTitleMap[opt.obs])
            dnH.GetYaxis().SetTitle("Ratio to reference")
            plot.add(dnH,systNameMap[syst]+" (down)",colors[isyst], False,False,False)
        plot.finalize()
        plot.show(outDir=opt.outdir,lumi=35900,noStack=True)

    fIn.Close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
