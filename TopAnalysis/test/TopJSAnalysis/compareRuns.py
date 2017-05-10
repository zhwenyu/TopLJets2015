import ROOT
import optparse
import os,sys
import json
from collections import OrderedDict
from math import sqrt

debug = True

"""
steer the script
"""
def main():
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/fill',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/result',              type='string')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=36500.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='all', help='flavor [default: %default]')
    (opt, args) = parser.parse_args()

    #read lists of samples
    samplesList=[]
    jsonList = opt.json.split(',')
    for jsonPath in jsonList:
        jsonFile = open(jsonPath,'r')
        samplesList += json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    
    data = {}   # data histograms
    
    for tag,sample in samplesList: 
        isData = sample[1]
        if not isData: continue

        if (not os.path.isfile('%s/%s.root'%(opt.inDir, tag))): continue
        fIn=ROOT.TFile.Open('%s/%s.root'%(opt.inDir, tag))
        if not fIn : continue
        
        h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix_py')

        if any(x in tag for x in ['2016B', '2016C', '2016D']):
            tag = 'BCD'
            color = ROOT.kRed+1
        if any(x in tag for x in ['2016E', '2016F']):
            tag = 'EF'
            color = ROOT.kBlue+1
        if any(x in tag for x in ['2016G', '2016H']):
            tag = 'GH'
            color = ROOT.kGreen+1
        try:
            data[tag].Add(h)
        except:
            data[tag] = h
            data[tag].SetDirectory(0)
            data[tag].SetLineColor(color)
        try:
            data['all'].Add(h)
        except:
            data['all'] = h.Clone('all')
            data['all'].SetDirectory(0)
            data['all'].SetLineColor(ROOT.kBlack)
 
    
    # Plot
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    
    legend = ROOT.TLegend(0.5,0.6,0.85,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    
    #mc.Draw()
    #legend.AddEntry(mc, 'MC', "p")
    isFirst = True
    for tag,hist in data.iteritems():
        print(tag, 'integral', hist.Integral(), 'mean', hist.GetMean())
        normalizeAndDivideByBinWidth(hist)
        if isFirst:
            hist.GetYaxis().SetRangeUser(0, hist.GetMaximum()*1.25)
            hist.Draw('hist')
        else:       hist.Draw('same hist')
        legend.AddEntry(hist, '%s, <N> = %.3f'%(tag, hist.GetMean()), "l")
        isFirst = False
    
    legend.Draw()
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_runs.pdf')

def normalizeAndDivideByBinWidth(hist):
    hist.Scale(1./hist.Integral())
    for i in range(1, hist.GetNbinsX()+1):
        hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
        hist.SetBinError  (i, hist.GetBinError(i)  /hist.GetBinWidth(i))
    return hist
        
"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

