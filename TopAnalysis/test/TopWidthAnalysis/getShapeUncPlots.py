import sys
import ROOT

ROOT.gROOT.SetBatch()

import CMS_lumi
import tdrStyle

from optparse import OptionParser
parser = OptionParser(
    usage="%prog [options] [label=datacard.txt | datacard.txt]",
    epilog=""
    )
parser.add_option("-i",     type="string", dest="indir"  ,   default="/afs/cern.ch/work/e/ecoleman/public/TopWidth/TopWidth_era2016/datacards/shapes.root", help="file to look for dists in")
parser.add_option("--wid",  type="string", dest="widList",   default="1",   help="a list of widths to look for in distnames")
parser.add_option("--lfs",  type="string", dest="lfsList",   default="EM",   help="a list of lepton final states to look for in distnames")
parser.add_option("--lbCh", type="string", dest="lbCatList", default="highpt,lowpt",   help="a list of states to look for in distnames")
parser.add_option("--cats", type="string", dest="catList",   default="1b,2b",   help="a list of lepton final states to look for in stats filenames")
parser.add_option("--uncs", type="string", dest="uncList",   default="PDF",   help="a list of uncertainties to look for")
parser.add_option("--proc", type="string", dest="procList",  default="tbart",   help="a list of processes to plot")
parser.add_option("-o",     type="string", dest="outdir" ,   default="./",   help="the base filename for the plots")
parser.add_option("--obs",
        type="string",
        dest="obsList",
        default="incmlb",
        help="a list of observable distributions to consider")

(options, args) = parser.parse_args()

def main():
    tdrStyle.setTDRStyle()

    url = options.indir

    systNameMap= {
            "jes"           : "Jet energy scale",
            "jer"           : "Jet energy resolution",
            "les"           : "Lepton energy scale",
            "ltag"          : "Lepton tagging efficiency",
            "trig"          : "Trigger efficiencies",
            "sel"           : "Selection efficiencies",
            "toppt"         : "Top p_{T}",
            "pu"            : "Pileup",
            "btag"          : "b-tagging efficiencies",
            "ttPartonShower": "t#bar{t} parton shower scale",
            "tWQCDScale"    : "tW QCD Scale",
            "tWttinterf"    : "tW/t#bar{t} interference",
            "Mtop"          : "Top mass",
            "MEmuF"         : "ME: Factorization scale",
            "MEmuR"         : "ME: Renormalization scale",
            "MEtot"         : "ME: Combined variation",
            "Herwig"        : "Hardonizer choice",
            "PDF"           : "NNLOPDF3.0 variation"
    }

    widList   = options.widList.split(',')
    lbCatList = options.lbCatList.split(',')
    lfsList   = options.lfsList.split(',')
    catList   = options.catList.split(',')
    procList  = options.procList.split(',')
    uncList   = options.uncList.split(',')
    obsList   = options.obsList.split(',')

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    fIn=ROOT.TFile.Open(url)

    can=ROOT.TCanvas('c','c',500,500)
    can.SetRightMargin(0.05)
    can.SetTopMargin(0.05)
    can.SetBottomMargin(0.1)
    can.SetLeftMargin(0.12)

    for wid,lbCat,lfs,cat,proc,obs in [(a,b,c,d,e,f)
            for a in widList
            for b in lbCatList
            for c in lfsList
            for d in catList
            for e in procList
            for f in obsList] :

        print "Starting for %s%s%s_%s%s"%(lbCat,lfs,cat,proc,wid)

        nomH=fIn.Get('%s%s%s_%s/%s%sw'%(lbCat,lfs,cat,obs,proc,("%2.1f"%float(wid)).replace('.','p')))
        nomH.SetLineWidth(2)
        nomH.SetLineColor(1)
        nomH.SetFillStyle(0)
        nomH.Draw('hist')
        nomH.GetYaxis().SetTitle("")
        nomH.SetMaximum(nomH.GetMaximum()*1.4)
        nomH.GetXaxis().SetTitle("Mass(lepton, jet) [GeV]")

        leg= ROOT.TLegend(0.60,0.55,0.75,0.9)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.AddEntry(nomH,'Nominal','l')

        colors=[ROOT.kRed+1,ROOT.kBlue,ROOT.kOrange+8,ROOT.kGreen-4,
                ROOT.kBlue-9,ROOT.kGreen+2,
                ROOT.kMagenta+2,ROOT.kMagenta,ROOT.kMagenta-9]
        allGr=[]
        nSysts=len(uncList)
        for isyst in xrange(0,nSysts):
            key=uncList[isyst]
            color=colors[isyst]
            print '%s%s%s_%s_%sUp/%s%sw'%(lbCat,lfs,cat,obs,key,proc,("%2.1f"%float(wid)).replace('.','p'))
            hup=fIn.Get('%s%s%s_%s_%sUp/%s%sw'%(lbCat,lfs,cat,obs,key,proc,("%2.1f"%float(wid)).replace('.','p')))
            hdn=fIn.Get('%s%s%s_%s_%sDown/%s%sw'%(lbCat,lfs,cat,obs,key,proc,("%2.1f"%float(wid)).replace('.','p')))
            allGr.append( ROOT.TGraphAsymmErrors() )
            for xbin in xrange(1,hup.GetNbinsX()+1):
                valUp=hup.GetBinContent(xbin)
                valDn=hdn.GetBinContent(xbin)
                valCen=nomH.GetBinContent(xbin)
                binw=hdn.GetXaxis().GetBinWidth(xbin)
                binc=hdn.GetXaxis().GetBinCenter(xbin)

                dx=(binw*0.9/nSysts)
                if isyst>0:
                    if isyst%2:
                        binc -= dx*(isyst/2+1)
                    else:
                        binc += dx*(isyst/2)

                allGr[-1].SetPoint(xbin,binc,valCen)
                ylo=ROOT.TMath.Min(valUp-valCen,valDn-valCen)
                yhi=ROOT.TMath.Max(valUp-valCen,valDn-valCen)
                allGr[-1].SetPointError(xbin,dx*0.5,dx*0.5,ROOT.TMath.Abs(ylo),ROOT.TMath.Abs(yhi))

            allGr[-1].SetTitle(key)
            allGr[-1].SetName(key)
            allGr[-1].SetMarkerStyle(0)
            allGr[-1].SetMarkerColor(color)
            allGr[-1].SetFillStyle(1001)
            allGr[-1].SetFillColor(color)
            allGr[-1].Draw('2')
            leg.AddEntry(allGr[-1],key if not systNameMap[key] else systNameMap[key],'f')

        leg.Draw()

        CMS_lumi.relPosX = 0.15
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.CMS_lumi(can,4,0)
        can.Modified()
        can.Update()
        can.SaveAs('%s/UncertaintiesPDF_%s%s%s_%s_%s%sw.pdf'%(options.outdir,lbCat,lfs,cat,obs,proc,("%2.1f"%float(wid)).replace('.','p')))
        can.SaveAs('%s/UncertaintiesPDF_%s%s%s_%s_%s%sw.png'%(options.outdir,lbCat,lfs,cat,obs,proc,("%2.1f"%float(wid)).replace('.','p')))

    fIn.Close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
