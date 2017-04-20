#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *
from math import sqrt

flavorMap = {'all': -1, 'light': 1, 'bottom': 5, 'gluon': 0}

def fillHist(hmap, tree, weightindex = 0):
    for event in tree:
        for obs,hmap_reco in hmap.iteritems():
            for reco,hmap_flavor in hmap_reco.iteritems():
                for f,h in hmap_flavor.iteritems():
                    flavor = flavorMap[f]
                    #TODO: weights
                    weight = event.weight[weightindex]
                    #TODO: run updated selection in TOPJetShape.cc
                    
                    for j in range(event.nj):
                        if (event.reco_sel != 1
                            or (flavor > -1 and event.j_flavor[j] != flavor)):
                            valReco = -1
                        else:
                            valReco = eval('event.j_'+obs+'_'+reco)[j]
                        g = event.j_gj[j]
                        if g >= 0: #match
                            if (event.gen_sel != 1
                                or (flavor > -1 and event.gj_flavor[g] != flavor)):
                                valGen = -1
                            else:
                                valGen = eval('event.gj_'+obs+'_'+reco)[g]
                        else: valGen = -1
                        h.Fill(valGen, valReco, weight)
                    for g in range(event.ngj):
                        if event.gj_j[g] < 0: #no match
                            if (event.gen_sel != 1
                                or (flavor > -1 and event.gj_flavor[g] != flavor)):
                                valGen = -1
                            else:
                                valGen = eval('event.gj_'+obs+'_'+reco)[g]
                            h.Fill(valGen, -1, weight)
                    
                    '''
                    reco = -1:
                      - event.gen_sel == 1 and event.reco_sel == -1
                        (all jets, do not throw any gen event based on reco inefficiency)
                      - genjets without reco match
                    gen = -1:
                      - event.gen_sel != 1 and event.reco_sel == 1 (=misreco background)
                      - reco jets without gen match
                    -> loop reco->gen first, put matched into matrix, put unmatched with gen=-1
                       afterwards loop gen->reco, ignore matched, put unmatched with reco=-1
                    If event is not selected in gen fill all jets with gen=-1.
                    If event is not selected in reco fill all jets with reco=-1.
                    '''
                    '''
                    if event.gen_sel == 1 and event.reco_sel == 1:
                        for j in range(event.nj):
                            valReco = eval('event.j_'+obs+'_'+reco)[j]
                            i = event.j_gj[j]
                            if i >= 0: #match
                                valGen = eval('event.gj_'+obs+'_'+reco)[i]
                            else: valGen = -1
                            h.Fill(valGen, valReco)
                        for j in range(event.ngj):
                            if event.gj_j[j] < 0: #no match
                                valGen = eval('event.gj_'+obs+'_'+reco)[j]
                                h.Fill(valGen, -1)
                    if event.gen_sel == 1 and event.reco_sel == -1:
                        for i in range(event.ngj):
                            valGen = eval('event.gj_'+obs+'_'+reco)[i]
                            h.Fill(valGen, -1)
                    if event.gen_sel == -1 and event.reco_sel == 1:
                        for i in range(event.nj):
                            valReco = eval('event.j_'+obs+'_'+reco)[i]
                            h.Fill(-1, valReco)            
                    '''


"""
steer
"""
def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--task',
                            dest='task',   
                            default='fill',
                            help='task: optimize,fill [default: %default]')
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='analysis.root',
                            help='input file, if directory the script will run in batch mode [default: %default]')
    parser.add_option('--obs',
                            dest='obs',   
                            default='mult',
                            help='observable [default: %default]')
    parser.add_option('-o', '--output',
                            dest='output', 
                            default='unfolding/fill/Chunks',
                            help='Output directory [default: %default]')
    parser.add_option('--ri', '--rootinput',
                            dest='rootinput',
                            default='unfolding/optimize/output.root',
                            help='output root file [default: %default]')
    parser.add_option('--skipexisting', 
                            dest='skipexisting',
                            help='skip jobs with existing output files  [%default]',
                            default=True, action='store_true')
    parser.add_option('-q', '--queue', dest='queue',  help='Batch queue to use [default: %default]', default='8nh')
    parser.add_option(      '--only',  dest='only',   help='csv list of samples to process (exact name, applies only for batch mode) [%default]', default=None, type='string')
    parser.add_option('-w', '--weightindex', dest='weightindex', help='weight index (single job) or number of weights to run (batch) [%default]', default=0, type='int')
    (opt, args) = parser.parse_args()

    if opt.input.endswith('.root'):
        rootoutfilepath = opt.output+'/'+os.path.basename(opt.input)
        if (opt.weightindex > 0):
            rootoutfilepath = rootoutfilepath.rsplit('_',1)[0] + "_wgt" + str(opt.weightindex) + "_" + rootoutfilepath.rsplit('_',1)[1]
        if (opt.skipexisting and os.path.isfile(rootoutfilepath)):
            print("skip existing file " + rootoutfilepath)
            return
        
        print("Filling response matrix from " + opt.input)
        rootinfile = ROOT.TFile(opt.rootinput, "READ");
        
        keys = []
        for tkey in rootinfile.GetListOfKeys():
            keys.append(tkey.GetName())
        #observables = ["mult", "width"]
        observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
        
        reco = ["charged"]
        
        flavors = ['all', 'light', 'bottom', 'gluon']
        
        hmap = {}
        for o in observables:
            hmap[o] = {}
            for r in reco:
                if o+'_'+r+'_responsematrix' in keys:
                    hmap[o][r] = {}
                    for f in flavors:
                        obj = rootinfile.Get(o+'_'+r+'_responsematrix')
                        obj.Reset()
                        hmap[o][r][f] = obj.Clone(o+'_'+r+'_'+f+'_responsematrix')
        
        tree = ROOT.TChain('tjsev')
        if opt.input == 'eos':
            opt.input = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_0.root'
        if opt.input == 'eosdata':
            opt.input = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/Data13TeV_SingleMuon_2016G_27.root'
        if (tree.Add(opt.input) == 0): return
        
        fillHist(hmap, tree, opt.weightindex)

        os.system('mkdir -p %s' % opt.output)
        rootoutfile = ROOT.TFile(rootoutfilepath, "RECREATE");    
        rootoutfile.cd()    
        for obs,hmap_reco in hmap.iteritems():
            for reco,hmap_flavor in hmap_reco.iteritems():
                for flavor,h in hmap_flavor.iteritems():
                    h.Write()
                    h.ProjectionX().Write()
                    h.ProjectionY().Write()
    else: # got input directory, go to batch
        cmsswBase = os.environ['CMSSW_BASE']
        workdir   = cmsswBase+'/src/TopLJets2015/TopAnalysis/'
        #parse selection lists
        onlyList=[]
        try:
            onlyList=opt.only.split(',')
        except:
            pass
        #create the tasklist
        inputlist=[]
        if os.path.isdir(opt.input):
            for file_path in os.listdir(opt.input):
                if file_path.endswith('.root'):
                    rootoutfilepath = workdir+opt.output+'/'+os.path.basename(file_path)
                    if (opt.weightindex == 0 and opt.skipexisting and os.path.isfile(rootoutfilepath)): continue
                    #filter tags
                    if len(onlyList)>0:
                        if (not os.path.basename(file_path).rsplit('_',1)[0] in onlyList): continue
                    inputlist.append(os.path.join(opt.input,file_path))
        print(inputlist)
        njobs = len(inputlist) * (opt.weightindex+1)
        print 'Running %d jobs to %s'%(njobs,opt.queue)
        njob = 0
        for inputfile in inputlist:
            for wgt in range(opt.weightindex+1):
                njob+=1
                if (wgt == 12): continue # weight for (mur=1,muf=1)=default
                outputdir=workdir+opt.output
                rootoutfilepath = workdir+opt.output+'/'+os.path.basename(inputfile)
                if (wgt > 0):
                    rootoutfilepath = rootoutfilepath.rsplit('_',1)[0] + "_wgt" + str(wgt) + "_" + rootoutfilepath.rsplit('_',1)[1]
                if (opt.skipexisting and os.path.isfile(rootoutfilepath)): continue
                localRun='python %s/src/TopLJets2015/TopAnalysis/test/TopJSAnalysis/fillUnfoldingMatrix.py --input %s --output %s --weightindex %d'%(cmsswBase,inputfile,outputdir,wgt)
                cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
                print 'Submitting job %d/%d: %s, weight %d -> %s'%(njob,njobs,os.path.basename(inputfile),wgt,os.path.basename(rootoutfilepath))
        #    os.system(cmd)

if __name__ == "__main__":
	sys.exit(main())
