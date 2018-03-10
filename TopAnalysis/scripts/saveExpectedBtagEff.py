import ROOT
from array import array
from TopLJets2015.TopAnalysis.storeTools import *
import optparse
import sys

"""
"""
def saveExpectedBtagEff(opt):

    #prepare histograms
    ptBins = [0,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,800,1000]
    preTagH=ROOT.TH1F('preTagH',';p_{T} [GeV];',len(ptBins)-1,array('d',ptBins))
    preTagH.Sumw2()
    tagH=preTagH.Clone('tagH')

    taggerList=[ x for x in opt.taggers.split(',')]
    flavConds=[('b',"abs(j_hadflav)==5"),
              ('c',"abs(j_hadflav)==4"),
              ('udsg','abs(j_hadflav)!=5 && abs(j_hadflav)!=4'),
              ('pu','abs(j_hadflav)!=5 && abs(j_hadflav)!=4 && abs(j_g)<0')]

    fOut=ROOT.TFile.Open(opt.output+'/expTagEff.root','RECREATE')

    for x in taggerList:

        print '\t Starting',x
        cut='j_%s>%s'%x.split(':')
        inputDir=opt.input
    
        #open a file
        input_list=getEOSlslist(directory=inputDir)
        data=ROOT.TChain('analysis/data') 
        for i in xrange(0,min(5,len(input_list))):
            data.Add(input_list[i])

            print 'Projecting tagging efficiency from ',data.GetEntries(),' events'
            print 'Working point is',tagger,wp
        
        #count number of tagged jets
        effgrs=[]
        for flav,cond in flavConds:
            print '\t computing for',flav,cond
            preTagH.Reset('ICE')
            tagH.Reset('ICE')
            data.Draw("j_pt >> preTagH",cond,'goff')
            data.Draw('j_pt >> tagH',cond + ' && ' + csvWP,'goff')
            effgrs.append(ROOT.TGraphAsymmErrors())
            effgrs[-1].SetName(flav)
            effgrs[-1].Divide(tagH,preTagH)

        #write to file
        fOut.cd()
        outDir=fOut.mkdir(x)
        for gr in effgrs:
            gr.Write()
        fOut.cd()

    fOut.Close()


"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',        dest='input',    help='input directory with files', default='/store/cmst3/group/top/RunIIFall17/c29f431')
    parser.add_option('-o', '--out',       dest='output',   help='output directory',           default='data/era2017')
    parser.add_option(      '--taggers',   dest='taggers',  help='tagger:cut,tagger:cut',      default='deepCSV:0.4941,csv:0.8838')
    (opt, args) = parser.parse_args()

    saveExpectedBtagEff(opt)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
