import os
import sys
import optparse
import ROOT
import pickle
import json
import re
import commands
from TopLJets2015.TopAnalysis.storeTools import *
from TopLJets2015.TopAnalysis.batchTools import *

"""
Wrapper to be used when run in parallel
"""
def RunMethodPacked(args):

    method,inF,outF,channel,charge,flav,runSysts,systVar,era,tag,debug=args
    print 'Running ',method,' on ',inF
    print 'Output file',outF
    print 'Selection ch=',channel,' charge=',charge,' flavSplit=',flav,' systs=',runSysts
    print 'Normalization applied from tag=',tag
    print 'Corrections will be retrieved for era=',era

    try:
        cmd='analysisWrapper --era %s --normTag %s --in %s --out %s --method %s --charge %d --channel %d --flav %d --systVar %s'\
            %(era, tag, inF, outF, method, charge, channel, flav, systVar)
        if runSysts : cmd += ' --runSysts'
        if debug : cmd += ' --debug'
        print(cmd)
        os.system(cmd)

    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inF)
        print 50*'<'
        return False
    return True

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-m', '--method',      dest='method',      help='method to run [%default]',                   default='TOP-16-006::RunTop16006',  type='string')
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file [%default]',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)  [%default]',  default='analysis', type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process  [%default]',             default=None,       type='string')
    parser.add_option(      '--skip',        dest='skip',        help='csv list of samples to skip  [%default]',             default=None,       type='string')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--systVar',     dest='systVar',     help='single systematic variation  [%default]',   default='nominal',       type='string')
    parser.add_option(      '--debug',       dest='debug',      help='debug mode  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--flav',        dest='flav',        help='split according to heavy flavour content  [%default]',   default=0,          type=int)
    parser.add_option(      '--ch',          dest='channel',     help='channel  [%default]',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge  [%default]',                                     default=0,          type=int)
    parser.add_option(      '--era',         dest='era',         help='era to use for corrections/uncertainties  [%default]',   default='era2016',       type='string')
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag  [%default]',                    default=None,       type='string')
    parser.add_option('-q', '--queue',       dest='queue',       help='submit to this queue  [%default]',                       default='local',    type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel  [%default]',                                default=0,    type='int')
    parser.add_option(      '--skipexisting',dest='skipexisting',help='skip jobs with existing output files  [%default]',                            default=False,      action='store_true')
    (opt, args) = parser.parse_args()

    #parse selection lists
    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass
    skipList=[]
    try:
        skipList=opt.skip.split(',')
    except:
        pass
    #parse list of systematic variations
    varList=[]
    if opt.systVar == 'all':
        allSystVars = ['jec_CorrelationGroupMPFInSitu', 'jec_CorrelationGroupInterCalibration',
                       'jec_CorrelationGroupUncorrelated', 'jec_FlavorPureGluon', 'jec_FlavorPureQuark',
                       'jec_FlavorPureCharm', 'jec_FlavorPureBottom', 'jer',
                       'btag_heavy', 'btag_light', 'csv_heavy', 'csv_light']
        for var in allSystVars:
            varList.append(var+'_up')
            varList.append(var+'_down')
    else:
        try:
            varList=opt.systVar.split(',')
        except:
            pass
    print 'Running following variations: ', varList

    #prepare output if a directory
    if not '.root' in opt.output :
        print opt.output
        if not '/store/' in opt.output:
            os.system('mkdir -p %s/Chunks'%opt.output)
        else:
            os.system('eos mkdir %s'%opt.output)
            os.system('eos mkdir %s/Chunks'%opt.output)
    #correct location of corrections to be used using cmsswBase, if needed
    cmsswBase=os.environ['CMSSW_BASE']
    if not cmsswBase in opt.era : opt.era=cmsswBase+'/src/TopLJets2015/TopAnalysis/data/'+opt.era

    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        inF=opt.input
        if '/store/' in inF and not 'root:' in inF : inF='root://eoscms//eos/cms'+opt.input        
        for systVar in varList:
            outF=opt.output
            if systVar != 'nominal': outF=opt.output[:-5]+'_'+systVar+'.root'
            task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,opt.flav,opt.runSysts,systVar,opt.era,opt.tag,opt.debug) )
    else:

        inputTags=getEOSlslist(directory=opt.input,prepend='')
        for baseDir in inputTags:

            tag=os.path.basename(baseDir)
            if tag=='backup' : continue

            #filter tags
            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if itag in tag:
                        processThisTag=True
                        break
                if not processThisTag : continue
            if len(skipList)>0:
                processThisTag=True
                for itag in skipList:
                    if itag in tag:
                        processThisTag=False
                        break
                if not processThisTag : continue

            input_list=getEOSlslist(directory='%s/%s' % (opt.input,tag) )
            nexisting = 0
            
            for ifile in xrange(0,len(input_list)):
                inF=input_list[ifile]
                for systVar in varList:
                    outF=os.path.join(opt.output,'Chunks','%s_%d.root' %(tag,ifile))
                    if systVar != 'nominal': outF=os.path.join(opt.output,'Chunks','%s_%s_%d.root' %(tag,systVar,ifile))
                    if (opt.skipexisting and os.path.isfile(outF)):
                        nexisting += 1
                        continue
                    task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,opt.flav,opt.runSysts,systVar,opt.era,tag,opt.debug) )
            if (opt.skipexisting and nexisting): print '--skipexisting: skipping %d of %d tasks as files already exist'%(nexisting,len(input_list))

    #run the analysis jobs
    if opt.queue=='local':
        print 'launching %d tasks in %d parallel jobs'%(len(task_list),opt.njobs)
        if opt.njobs == 0:
            for args in task_list: RunMethodPacked(args)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        print 'launching %d tasks to submit to the %s queue'%(len(task_list),opt.queue)        

        FarmDirectory                      = opt.output+"/FARM"
        if '/store' in FarmDirectory : FarmDirectory = './FARM%s'%os.path.basename(opt.output)
        os.system('mkdir -p %s'%FarmDirectory)

        jobNb=0
        for method,inF,outF,channel,charge,flav,runSysts,systVar,era,tag,debug in task_list:
            jobNb+=1
            cfgfile='%s/job_%s.sh'%(FarmDirectory,os.path.splitext(os.path.basename(outF))[0])
            logfile='%s/job_%s.log'%(FarmDirectory,os.path.splitext(os.path.basename(outF))[0])
            cfg=open(cfgfile,'w')
            cfg.write('WORKDIR=`pwd`\n')
            cfg.write('echo "Working directory is ${WORKDIR}"\n')
            cfg.write('cd %s\n'%cmsswBase)
            cfg.write('eval `scram r -sh`\n')
            cfg.write('cd ${WORKDIR}\n')
            localOutF=os.path.basename(outF)
            runOpts='-i %s -o ${WORKDIR}/%s --charge %d --ch %d --era %s --tag %s --flav %d --method %s --systVar %s'\
                    %(inF, localOutF, charge, channel, era, tag, flav, method, systVar)
            if runSysts : runOpts += ' --runSysts'
            if debug :    runOpts += ' --debug'
            cfg.write('python %s/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py %s &> run.log\n'%(cmsswBase,runOpts))
            if '/store' in outF:
                cfg.write('xrdcp ${WORKDIR}/%s root://eoscms//eos/cms/%s\n'%(localOutF,outF))
                cfg.write('rm ${WORKDIR}/%s'%localOutF)
            elif outF!=localOutF:
                cfg.write('mv -v ${WORKDIR}/%s %s\n'%(localOutF,outF))
                cfg.write('mv -v ${WORKDIR}/run.log %s\n'%(logfile))
            cfg.close()
            os.system('chmod u+x %s'%cfgfile)
            print 'Submitting job %d/%d'%(jobNb, len(task_list))
            os.system('bsub -q %s %s -R "pool>30000"'%(opt.queue,
                                                       os.path.abspath(cfgfile)) )
        


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
