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
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option( '--systJson', dest='systJson', help='json with list of systematics', default='data/era2016/syst_samples.json', type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/fill',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/result',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
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
    
    #read lists of syst samples
    systSamplesList=None
    if opt.systJson:
        jsonFile = open(opt.systJson,'r')
        systSamplesList=json.load(jsonFile,encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()

    allSystVars = ['jec_CorrelationGroupMPFInSitu',
                   'jec_RelativeFSR',
                   'jec_CorrelationGroupUncorrelated',
                   'jec_FlavorPureGluon',
                   'jec_FlavorPureQuark',
                   'jec_FlavorPureCharm',
                   'jec_FlavorPureBottom',
                   'jer',
                   'btag_heavy',
                   'btag_light',
                   'csv_heavy',
                   'csv_light'
                  ]
    varList = []
    for var in allSystVars:
        #varList.append(var+'_up')
        varList.append(var+'_down')
    expSystSamplesList = []
    # TODO: need to reprocess nominal samples
    #for var in varList:
    #    expSystSamplesList.append(['MC13TeV_TTJets_'+var, [832., 0., '', 't#bar{t} '+var]])
    
    data        = None # data histogram
    backgrounds = {}   # background histograms
    nominal     = None # reference unfolding matrix
    systematics = {}   # systematic unfolding matrices
    
    for slist,isSyst in [ (samplesList,False), (systSamplesList,True), (expSystSamplesList,True) ]:
        for tag,sample in slist: 
            if isSyst and not 't#bar{t}' in sample[3] : continue
            isData = sample[1]
            isSig  = not isSyst and sample[3] == 't#bar{t}'
            isBkg  = not (isData or isSig or isSyst)
            xs = sample[0]

            if (not os.path.isfile('%s/%s.root'%(opt.inDir, tag))): continue
            fIn=ROOT.TFile.Open('%s/%s.root'%(opt.inDir, tag))
            if not fIn : continue
            
            if isData:
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix_py')
                try:
                    data.Add(h)
                except:
                    data=h
                    data.SetDirectory(0)
                if debug: print(tag, data.GetEntries())
            
            if isBkg:
                #if not tag in ['MC13TeV_SingleTbar_tW', 'MC13TeV_SingleT_tW', 'MC13TeV_SingleTbar_t', 'MC13TeV_SingleT_t']: continue
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix_py')
                h.Scale(opt.lumi*xs)
                integral = h.Integral()
                integralError = 0
                integralRelError = 999.
                for i in range(1, h.GetNbinsX()+1):
                    integralError = sqrt(integralError**2 + h.GetBinError(i)**2)
                if integral > 0: integralRelError = integralError/integral
                if debug: print(tag, integral, integralError, integralRelError)
                if integralRelError > 1:
                    if debug: print(tag, 'Relative integral error > 1, skipping sample')
                    continue
                try:
                    backgrounds[tag].Add(h)
                except:
                    backgrounds[tag]=h
                    backgrounds[tag].SetDirectory(0)
                if debug:
                    integral = backgrounds[tag].Integral()
                    integralError = 0
                    integralRelError = 0
                    for i in range(1, backgrounds[tag].GetNbinsX()+1):
                        integralError = sqrt(integralError**2 + backgrounds[tag].GetBinError(i)**2)
                    if integral > 0: integralRelError = integralError/integral
                    print(tag, integral, integralError, integralRelError)
            
            if isSig:
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix')
                h.Scale(opt.lumi*xs)
                try:
                    nominal.Add(h)
                except:
                    nominal=h
                    nominal.SetDirectory(0)
                if debug: print(tag, nominal.GetEntries())
                # signal file should contain systematics weights
                for i in range(1, 21):
                    ih=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_wgt'+str(i)+'_responsematrix')
                    itag = tag + '_wgt' + str(i)
                    try:
                        systematics[itag].Add(ih)
                    except:
                        systematics[itag]=ih
                        systematics[itag].SetDirectory(0)
                    if debug: print(itag, systematics[itag].GetEntries())
            
            if isSyst:
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix')
                try:
                    systematics[tag].Add(h)
                except:
                    systematics[tag]=h
                    systematics[tag].SetDirectory(0)
                if debug: print(tag, systematics[tag].GetEntries())

    # TODO: put this into def to call for each systematic
    dataUnfolded, dataFoldedBack, dataBkgSub = unfold('MC13TeV_TTJets', nominal, backgrounds, data)
    systematicUnfolded = {}
    for tag,systematic in systematics.iteritems():
        systematicUnfolded[tag] = unfold(tag, systematic, backgrounds, data)[0]
    
    systUp=[0.]
    systDown=[0.]
    for i in range(1, dataUnfolded.GetNbinsX()+1):
        systUp.append(0.)
        systDown.append(0.)
        for tag,systematic in systematicUnfolded.iteritems():
            diff = systematic.GetBinContent(i) - dataUnfolded.GetBinContent(i)
            if (diff > 0):
                systUp[i] = sqrt(systUp[i]**2 + diff**2)
            else:
                systDown[i] = sqrt(systDown[i]**2 + diff**2)

    dataUnfoldedSys = dataUnfolded.Clone('dataUnfoldedSys')    
    for i in range(1, dataUnfoldedSys.GetNbinsX()+1):
        dataUnfoldedSys.SetBinContent(i, dataUnfolded.GetBinContent(i) + (systUp[i]-systDown[i])/2.)
        dataUnfoldedSys.SetBinError(i, sqrt(dataUnfolded.GetBinError(i)**2 + ((systUp[i]+systDown[i])/2.)**2))
    
    # Plot
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    
    dataUnfolded.SetTitle('')
    dataUnfolded.GetYaxis().SetRangeUser(0., 1.5*dataUnfolded.GetMaximum())
    dataUnfolded.SetLineColor(ROOT.kBlack)
    dataUnfolded.SetMarkerColor(ROOT.kBlack)
    dataUnfolded.SetMarkerStyle(20)
    dataUnfolded.Draw('P X0 E1')
    dataUnfoldedSys.SetLineWidth(2)
    dataUnfoldedSys.SetLineColor(ROOT.kBlack)
    dataUnfoldedSys.SetMarkerColor(ROOT.kBlack)
    dataUnfoldedSys.SetMarkerStyle(1)
    dataUnfoldedSys.Draw('SAME P X0 E1')
    
    nominalGen = normalizeAndDivideByBinWidth(nominal.ProjectionX("nominalGen"))
    nominalGen.GetYaxis().SetRangeUser(0., 1.5*nominalGen.GetMaximum())
    nominalGen.SetLineColor(ROOT.kRed+1)
    nominalGen.SetLineColor(ROOT.kRed+1)
    nominalGen.SetMarkerColor(ROOT.kRed+1)
    nominalGen.SetMarkerStyle(5)
    nominalGen.Draw('SAME P X0 E1')
    
    nominalReco = normalizeAndDivideByBinWidth(nominal.ProjectionY("nominalReco"))
    nominalReco.SetLineColor(ROOT.kRed+1)
    nominalReco.SetLineStyle(2)
    nominalReco.Draw("SAME,HIST")
    
    normalizeAndDivideByBinWidth(dataBkgSub)
    dataBkgSub.SetLineColor(ROOT.kBlack)
    dataBkgSub.SetLineStyle(0)
    dataBkgSub.Draw("SAME,HIST")
    
    dataFoldedBack.SetMarkerColor(ROOT.kMagenta+1);
    dataFoldedBack.SetMarkerStyle(2);
    dataFoldedBack.Draw("SAME P X0 E1");
    
    legend = ROOT.TLegend(0.5,0.6,0.85,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(nominalGen, "MC gen t#bar{t}", "p")
    legend.AddEntry(nominalReco, "MC reco t#bar{t}", "l")
    legend.AddEntry(dataBkgSub, "data (bg-sub)", "l")
    legend.AddEntry(dataUnfolded, "data unfolded", "ep")
    legend.AddEntry(dataFoldedBack, "data folded back", "p")
    legend.Draw()
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_test.pdf')



def unfold(Mtag, M, backgrounds, data):
    # Background subtraction
    dataBkgSub = data.Clone(Mtag+'_dataBkgSub')

    # TODO: background uncertainties
    for tag,background in backgrounds.iteritems():
        dataBkgSub.Add(background, -1.)
    #if debug:
    #    print('dataBkgSub', dataBkgSub.GetEntries(), dataBkgSub.Integral())
    #    for i in range(dataBkgSub.GetNbinsX()+2):
    #        print('dataBkgSub bin content', dataBkgSub.GetBinContent(i), 'bin error', dataBkgSub.GetBinError(i))
    
    sigRecoNonGen = M.ProjectionY("sigRecoNonGen", 0, 0);
    sigReco       = M.ProjectionY("sigReco");
    
    for i in range(dataBkgSub.GetNbinsX()+2):
        sf = sigRecoNonGen.GetBinContent(i) / sigReco.GetBinContent(i)
        dataBkgSub.SetBinContent(i, (1.-sf)*dataBkgSub.GetBinContent(i))
        dataBkgSub.SetBinError  (i, (1.-sf)*dataBkgSub.GetBinError(i))
        M.SetBinContent(0, i, 0.)
    #if debug:
    #    print('dataBkgSub', dataBkgSub.GetEntries(), dataBkgSub.Integral())
    #    for i in range(dataBkgSub.GetNbinsX()+2):
    #        print('dataBkgSub bin content', dataBkgSub.GetBinContent(i), 'bin error', dataBkgSub.GetBinError(i))
    
    # Do unfolding
    unfold = ROOT.TUnfoldDensity(M, ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfold.kRegModeCurvature, ROOT.TUnfold.kEConstraintArea, ROOT.TUnfoldDensity.kDensityModeBinWidth)
    
    if (unfold.SetInput(dataBkgSub) >= 10000):
        print('SetInput bad return value >= 10000')
        return
    
    unfold.DoUnfold(0.);
    
    # Retrieve results
    dataUnfolded = normalizeAndDivideByBinWidth(unfold.GetOutput(Mtag+"_Unfolded"))
    
    #if debug:
    #    for i in range(1, dataUnfolded.GetNbinsX()+1):
    #        print('dataUnfolded', i, dataUnfolded.GetBinContent(i), dataUnfolded.GetBinError(i))
    dataFoldedBack = normalizeAndDivideByBinWidth(unfold.GetFoldedOutput(Mtag+"_FoldedBack"))
    
    return dataUnfolded, dataFoldedBack, dataBkgSub



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

