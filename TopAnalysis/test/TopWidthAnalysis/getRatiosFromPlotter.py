import optparse
import os,sys
import json
import ROOT
import math
import pickle

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default="./plotter.root",       type='string')
    parser.add_option('-d', '--divideFile',  dest='inDiv' ,      help='comparison input plotter',       default="./plotter.root",       type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='mcratio_plotter.root', type='string')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi to print out',              default=41.6,                   type=float)
    parser.add_option(      '--saveLog',     dest='saveLog' ,    help='save log versions of the plots', default=False,             action='store_true')
    parser.add_option('--wids',  dest='wids'  , help='widths to compare to nominal in div', default="4.0",                   type='string')
    parser.add_option('--obs',   dest='obs'   , help='observables to process', default="incmlb,mdrmlb,sncmlb,minmlb,mt2mlb", type='string')
    parser.add_option('--dists', dest='dists' , help='non-observable distributions to process', default="tmass",             type='string')
    parser.add_option('--sigs',  dest='sigs'  , help='signal processes to plot ratios for',     default="t#bar{t}",          type='string')
    parser.add_option('--ptCh',  dest='ptchs' , help='pt categories to consider',               default="highpt,lowpt",      type='string')
    parser.add_option('--ch',    dest='chs'   , help='final states to consider',                default="EE,EM,MM",          type='string')
    parser.add_option('--bcat',  dest='bcats' , help='b categories to consider',                default="1b,2b",             type='string')
    (opt, args) = parser.parse_args()

    wgtFin=ROOT.TFile(opt.inDir)
    divFin=ROOT.TFile(opt.inDiv)


    # get lists for distribution collection
    widList =opt.wids.split(',')
    ptChList=opt.ptchs.split(',')
    chList  =opt.chs.split(',')
    bcatList=opt.bcats.split(',')
    distList=opt.dists.split(',')
    obsList =opt.obs.split(',')
    sigList =opt.sigs.split(',')

    # useful settings (currenty magic)
    divWid="4.0"

    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    outDir='./genvalidation'
    os.system('mkdir -p %s' % outDir)
    os.system('rm %s/%s'%(outDir,opt.outName))

    # collect a list of all distributions to plot
    plots = [("%s_%sw/%s_%sw_%s"%(dist,wid,dist,wid,sig),
        "%s_%sw/%s_%sw_%s"%(dist,divWid,dist,divWid,sig))
                for dist in distList
                for wid in widList
                for sig in sigList]
    plots+= [("%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s"%(ptC,ch,bc,obs,wid,ptC,ch,bc,obs,wid,sig),
        "%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s"%(ptC,ch,bc,obs,divWid,ptC,ch,bc,obs,divWid,sig))
                for ptC in ptChList
                for ch in chList
                for bc in bcatList
                for obs in obsList
                for wid in widList
                for sig in sigList]

    import tdrStyle
    import CMS_lumi

    tdrStyle.setTDRStyle()
    #CMS_lumi.lumiText=
    CMS_lumi.lumiTextSize=0.6
    CMS_lumi.extraText="Simulation Preliminary"
    CMS_lumi.extraOverCmsTextSize=0.5
    CMS_lumi.relPosX=0.18

    for p in plots :
        canvas=ROOT.TCanvas()
        if opt.saveLog : canvas.SetLogy()
        canvas.cd()

        wgtHist=wgtFin.Get(p[0])
        divHist=divFin.Get(p[1])

        divHist.Divide(wgtHist)

        divHist.GetYaxis().SetTitle("Events (4 #times#Gamma_{SM}) / Events (Reweighted)")
        divHist.GetYaxis().SetRangeUser(0.7,1.3)
        divHist.SetTitle("")

        divHist.Draw()

        CMS_lumi.CMS_lumi(canvas,4,0)
        canvas.SaveAs("%s/%s.pdf"%(outDir,p[0].split('/')[1].replace('.','w').replace('#','').replace('{','').replace('}','')))



    print '-'*50
    print 'Plots and summary ROOT file can be found in %s' % outDir
    print '-'*50


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

