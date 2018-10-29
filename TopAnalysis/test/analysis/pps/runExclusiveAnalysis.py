#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import numpy
import array
import json
import pickle
from collections import OrderedDict

def isValidRunLumi(run,lumi,runLumiList):

    """checks if run is available and lumi section was certified"""

    #no run/lumi to select, all is good by default
    if not runLumiList:
        return True

    #run is available
    if not run in runLumiList: 
        return False

    #lumi is within certified ranges
    for lran in runLumiList[run]:
        try:
            if lumi>=lran[0] and lumi<=lran[1]:
                return True
        except:
            if lumi==lran[0]:
                return True

    #reached this far, nothing found in list
    return False

def getTracksPerRomanPot(tree):

    """loops over the availabe tracks in the event and groups them by roman pot id"""
    tkPerRP={23:[],123:[]}
    try:
        for itk in xrange(0,tree.nRPtk):
            rpid=tree.RPid[itk]
            tkPerRP[rpid].append( tree.RPfarcsi[itk] )
    except:
        pass
    return tkPerRP



def runExclusiveAnalysis(fileList,outFileName,runLumiList,ctrFileList):
    
    """event loop"""

    tree=ROOT.TChain('tree')
    for f in fileList:
        tree.AddFile(f)
    nEntries=tree.GetEntries()  
    print '....analysing',nEntries,'in',len(fileList),'files, with output @',outFileName

    #control
    ctrl_tree=ROOT.TChain('tree')
    for f in ctrFileList:
        ctrl_tree.Add(f)
    nctrl=ctrl_tree.GetEntries()
    rand=ROOT.TRandom2()

    #start histograms
    histos={}
    histos['nvtx'] = {'inc':ROOT.TH1F('nvtx',';Vertex multiplicity;Events',50,0,100)}
    histos['mll'] = {'inc':ROOT.TH1F('mll',';Dilepton invariant mass [GeV];Events',50,20,250)}
    histos['ptll'] = {'inc':ROOT.TH1F('ptll',';Dilepton transverse momentum [GeV];Events',50,0,250)}
    histos['ntk']  = {'inc':ROOT.TH1F('ntk',';Track multiplicity;Events',5,0,5)}
    histos['csi']  = {'inc':ROOT.TH1F('csi',';#xi;Events',50,0,0.3)}
    for name in histos:
        for cat in histos[name]:
            histos[name][cat].SetDirectory(0)
            histos[name][cat].Sumw2()
            
    def fillHisto(val,weight,key):
        name,cat=key
        if not cat in histos[name]:
            histos[name][cat]=histos[name]['inc'].Clone('%s_%s'%key)
            histos[name][cat].SetDirectory(0)
            histos[name][cat].Reset('ICE')
        histos[name][cat].Fill(val,weight)

    #loop over events in the tree and fill histos
    irand=0
    for i in xrange(0,nEntries):

        tree.GetEntry(i)

        if not isValidRunLumi(tree.run,tree.lumi,runLumiList):
            continue

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nEntries))))
        
        if tree.isSS : continue
        dilId=abs(int(tree.l1id*tree.l2id))
        if dilId==11*11 : cat='ee'
        if dilId==11*13 : cat='em'
        if dilId==13*13 : cat='mm'
        
        fillHisto(tree.nvtx,tree.evwgt,key=('nvtx',cat))
        fillHisto(tree.mll,tree.evwgt,key=('mll',cat))
        fillHisto(tree.llpt,tree.evwgt,key=('ptll',cat))

        rptks = getTracksPerRomanPot(tree)
        for rpid in rptks:
            fillHisto(len(rptks[rpid]),tree.evwgt,key=('ntk','%d'%rpid))
            for csi in rptks[rpid]:
                fillHisto(csi,tree.evwgt,key=('csi','%d'%rpid))


        #get an event to mix RP information
        while True:
            #irand=rand.Integer(nctrl)
            #ctrl_tree.GetEntry( irand )
            ctrl_tree.GetEntry( irand )
            irand=irand+1 if irand<nctrl else 0
            if not isValidRunLumi(ctrl_tree.run,ctrl_tree.lumi,runLumiList):
                continue
            break

        ctrl_rptks = getTracksPerRomanPot(ctrl_tree)
        for rpid in ctrl_rptks:
            fillHisto(len(ctrl_rptks[rpid]),tree.evwgt,key=('ntk','mix%d'%rpid))
            for csi in ctrl_rptks[rpid]:
                fillHisto(csi,tree.evwgt,key=('csi','%d'%rpid))

    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    fOut.cd()
    for name in histos:
        for cat in histos[name]:
            if histos[name][cat].GetEntries()==0 : 
                histos[name][cat].Delete()
                continue
            histos[name][cat].SetDirectory(fOut)
            histos[name][cat].Write()
    fOut.Close()

 
def runExclusiveAnalysisPacked(args):

    """wrapper for parallel execution"""

    try:
        runExclusiveAnalysis(*args)
    except Exception as e:
        print 50*'<'
        print "  Problem with", args[1], "continuing without"
        print e
        print 50*'<'
        return False
    
"""
Create analysis tasks
"""
def runAnalysisTasks(opt):

    def getGroupedTasks(injson,inDir):

        #read samples to process
        task_dict={}
        with open(injson,'r') as cachefile:
            samples=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
            for x in samples:
                task_dict[x[0]]=[]

        #group chunks matching the same name
        for file_path in os.listdir(inDir):
            file_name,ext=os.path.splitext(file_path)
            if ext != '.root' : continue

            #check if file tag is already in the list of samples to process
            lastTkn=file_name.rfind('_')
            tag=file_name[0:lastTkn]
            if not tag in task_dict: continue
            task_dict[tag].append( os.path.join(inDir,file_path) )
        
        return task_dict


    task_dict=getGroupedTasks(opt.json,     opt.input)
    ctrl_task=getGroupedTasks(opt.ctrlJson, opt.input)
    incCtrlSampleList=[]
    for x in ctrl_task: incCtrlSampleList += ctrl_task[x]
    ctrl_task={x.split('_')[1]:ctrl_task[x] for x in ctrl_task.keys()}

    #parse json file with list of run/lumi sections
    with open(opt.RPin,'r') as cachefile:
        runLumi=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
        runLumi={int(x[0]):x[1] for x in runLumi}

    #create the tasks
    import multiprocessing as MP
    pool = MP.Pool(opt.jobs)
    task_list=[]
    for x in task_dict.keys():
        runLumiList=None
        ctrlSampleList=incCtrlSampleList
        if 'Data' in x:
            runLumiList=runLumi        
            #era=x.split('_')[1]
            #ctrlSampleList=ctrl_task[era]
            
                
        task_list.append( (task_dict[x],
                           os.path.join(opt.output,x)+'.root',
                           runLumiList,
                           ctrlSampleList) )

    pool.map(runExclusiveAnalysisPacked,task_list)


def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                      dest='input',   
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/Chunks',
                      help='input directory with the files [default: %default]')
    parser.add_option('--jobs',
                      dest='jobs', 
                      default=8,
                      type=int,
                      help='# of jobs to process in parallel the trees [default: %default]')
    parser.add_option('--json',
                      dest='json', 
                      default='pps_samples.json',
                      type='string',
                      help='json with the files to process')
    parser.add_option('--ctrJson',
                      dest='ctrlJson', 
                      default='pps_jetht_samples.json',
                      type='string',
                      help='json with control samples for background estimation (event mixing)')
    parser.add_option('--RPin',
                      dest='RPin', 
                      default='combined_RPIN_CMS.json',
                      type='string',
                      help='json with the runs/lumi sections in which RP are in')
    parser.add_option('-o', '--output',
                      dest='output', 
                      default='analysis',
                      help='Output directory [default: %default]')
    (opt, args) = parser.parse_args()
    
    os.system('mkdir -p %s' % opt.output)
    runAnalysisTasks(opt)
        

if __name__ == "__main__":
    sys.exit(main())
