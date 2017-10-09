import sys
import ROOT

from optparse import OptionParser

ROOT.gROOT.SetBatch(True)

def main():

    parser = OptionParser(
        usage="%prog [options]",
        epilog=""
        )
    parser.add_option("--input", type="string", dest="input" , default="shapes.root", help="file to look for dists in")
    parser.add_option("--cats" , type="string", dest="cats"  , default="EM2bhighpt",  help="a list of categories")
    parser.add_option("-o"     , type="string", dest="outdir", default="./cmpplots",  help="the base filename for the plots")
    parser.add_option("--tag"  , type="string", dest="tag"   , default="exp",         help="the kind of systs to get")
    parser.add_option("--wid"  , type="string", dest="wid"   , default="100",         help="the width to get systs for")
    parser.add_option("--proc" , type="string", dest="procs" , default="t#bar{t}",    help="a list of processes to plot")
    parser.add_option("--systInput", type="string", dest="systinput" , default="shapes.root", help="file to look for dists in")
    (opt, args) = parser.parse_args()

    tag=opt.tag
    wid=opt.wid

    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    fIn=ROOT.TFile.Open(opt.input)

    canv=ROOT.TCanvas("","",600,1000)
    canv.cd()

    canv.SetTopMargin(0.01)
    canv.SetBottomMargin(0.075)
    canv.SetLeftMargin(0.15)
    canv.SetRightMargin(0.15)

    for cat,proc in [(a,b) for a in opt.cats.split(',') for b in opt.procs.split(',')] :
        signalHist = None
        ratioHist  = None
        canv.Clear()

        print "\n\nGetting %s_incmlb_w%s/%s_incmlb_w%s_%s"%(cat,wid,cat,wid,proc)
        signalHist = fIn.Get("%s_incmlb_w%s/%s_incmlb_w%s_%s"%(cat,wid,cat,wid,proc)).Clone()
        signalHist.SetDirectory(0)

        if tag == "sys" :
            systIn = ROOT.TFile(opt.systinput)

            dir = systIn.Get("%s_incmlb_w%s"%(cat,wid))
            ndists = len([x for x in dir.GetListOfKeys() if proc in x.GetName()])

            ratioHist = ROOT.TH2F("ratiohist","ratiohist",signalHist.GetNbinsX(),0,300,ndists,0,ndists)
            iKey=1
            for key in dir.GetListOfKeys() :
                hName=key.GetName()
                if proc not in hName : continue
                hist =dir.Get(hName)

                # fill histogram
                for iBin in range(0,ratioHist.GetNbinsX()+1) :
                    ratioHist.SetBinContent(iBin,iKey,hist.GetBinContent(iBin))

                # set y axis labels
                label=hName.split(" ")[-1]
                if "up"==label or "dn"==label :
                    label = hName.split(" ")[-2] + " " + label
                ratioHist.GetYaxis().SetBinLabel(iKey,label)

                iKey += 1

        else:
            ratioHist = fIn.Get("%s_incmlb_w%s_%s/%s_incmlb_w%s_%s_%s"%(cat,wid,tag,cat,wid,tag,proc)).Clone()
            ratioHist.SetDirectory(0)

        for iBin in range(0,ratioHist.GetNbinsX()+1) :
            nomContents=signalHist.GetBinContent(iBin)
            for iBinY in range(1,ratioHist.GetNbinsY()) :
                newContents = ratioHist.GetBinContent(iBin,iBinY)
                if nomContents == 0 :
                    newContents = 1
                else :
                    newContents /= nomContents

                if newContents < 10e-05 :
                    newContents = 0.001

                ratioHist.SetBinContent(iBin,iBinY,newContents);

        ratioHist.Draw("COLZ")
        ratioHist.GetZaxis().SetTitle("Ratio to Nominal")
        ratioHist.GetZaxis().SetTitleOffset(1.5)
        ratioHist.GetZaxis().SetRangeUser(0.75,1.15)
        ratioHist.GetYaxis().SetTitleOffset(2.25)
        if tag == "gen" :
            ratioHist.GetYaxis().SetRangeUser(0,115)
        #canv.SetLogz()
        canv.Update()
        canv.Modified()

        nproc = proc.replace(' ','')
        nproc = nproc.replace('#','')
        nproc = nproc.replace('{','')
        nproc = nproc.replace('}','')
        nproc = nproc.replace('+','')

        tex = ROOT.TLatex()
        tex.SetTextSize(0.04)
        tex.DrawLatexNDC(0.01,0.01,"%s %s %s"%(cat,wid,proc))

        canv.SaveAs("SystMap_%s_incmlb_w%s_%s_%s.pdf" %(cat,wid,tag,nproc))
        canv.SaveAs("SystMap_%s_incmlb_w%s_%s_%s.png" %(cat,wid,tag,nproc))
        canv.SaveAs("SystMap_%s_incmlb_w%s_%s_%s.C"   %(cat,wid,tag,nproc))
        canv.SaveAs("SystMap_%s_incmlb_w%s_%s_%s.root"%(cat,wid,tag,nproc))


    fIn.Close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
