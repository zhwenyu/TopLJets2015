import ROOT
from array import array
from TopLJets2015.TopAnalysis.storeTools import *
import optparse
import sys

def saveExpectedBtagEff(opt):
    """loops over events and computes expected tagging efficiency per flavour"""

    #prepare data chain
    inputDir=opt.input
    input_list=getEOSlslist(directory=inputDir)
    data=ROOT.TChain('analysis/data') 
    for i in xrange(0,min(5,len(input_list))):
        data.Add(input_list[i])
    print 'Projecting tagging efficiency from ',data.GetEntries(),' events'

    #prepare histograms
    ptBins = [0,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,800,1000]
    preTagH=ROOT.TH1F('preTagH',';p_{T} [GeV];',len(ptBins)-1,array('d',ptBins))
    preTagH.Sumw2()
    tagH=preTagH.Clone('tagH')


    #loop over taggers and flavours
    taggerList=[ x for x in opt.taggers.split(',')]
    flavConds=[('b',"abs(j_hadflav)==5"),
               ('c',"abs(j_hadflav)==4"),
               ('udsg','abs(j_hadflav)!=5 && abs(j_hadflav)!=4')]
    effgrs={}
    for x in taggerList:

        print '\t Starting',x
        tagger,op,opval=x.split(':')        
        cut='j_%s>%s'%('csv' if tagger=='CSVv2' else 'deepcsv',opval)
        
        #count number of tagged jets
        effgrs[(tagger,op)]={}
        for flav,cond in flavConds:
            print '\t computing for',flav,cond
            preTagH.Reset('ICE')
            tagH.Reset('ICE')
            data.Draw("j_pt >> preTagH",cond,'goff')
            data.Draw('j_pt >> tagH',cond + ' && ' + cut,'goff')
            effgrs[(tagger,op)][flav]=ROOT.TGraphAsymmErrors()
            effgrs[(tagger,op)][flav].Divide(tagH,preTagH)

    #save in file
    fOut=ROOT.TFile.Open(opt.output,'RECREATE')
    for key in effgrs:
        fOut.cd()
        outDir=fOut.mkdir('%s_%s'%key)
        outDir.cd()
        for flav in effgrs[key]: effgrs[key][flav].Write(flav)
        fOut.cd()
    fOut.Close()

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',        dest='input',    help='input directory with files', default='/store/cmst3/group/top/RunIIFall17/c29f431/MC13TeV_TTJets')
    parser.add_option('-o', '--out',       dest='output',   help='output directory',           default='data/era2017/btv')
    parser.add_option('--taggers',   dest='taggers',  help='tagger:cut,tagger:cut',      
                      default='DeepCSV:loose:0.1522,DeepCSV:medium:0.4941,DeepCSV:tight:0.8001,CSVv2:loose:0.5803,CSVv2:medium:0.8838,CSVv2:tight:0.9693')
    (opt, args) = parser.parse_args()

    saveExpectedBtagEff(opt)

if __name__ == "__main__":
    sys.exit(main())
