import os
import sys
import optparse
import ROOT
import pickle
from collections import OrderedDict
import json
import re
import commands

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)  [%default]',  default='/store/group/phys_top/TOP-17-015/', type='string')
    parser.add_option('--outList',           dest='outList',       help='file list to copy to output [%default]',     default='MiniEvents.root',    type='string')    
    parser.add_option('-q', '--queue',       dest='queue',       help='if not local send to batch with condor. queues are now called flavours, see http://batchdocs.web.cern.ch/batchdocs/local/submit.html#job-flavours   [%default]',     default='workday',    type='string')    
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel  [%default]',                  default=0,    type='int')
    parser.add_option(      '--opts',        dest='cmsswOpts',   help='cmssw options [%default]',             default="nevents=50000 asfsr=0.1365",       type='string')

    parser.add_option(      '--farmappendix',        dest='farmappendix',        help='Appendix to condor FARM directory [%default]',             default='FARMGEN',       type='string')
    (opt, args) = parser.parse_args()

    #read env variable
    cmsswBase=os.environ['CMSSW_BASE']

    #prepare output if a directory
    if not '/store/' in opt.output:
        os.system('mkdir -p %s/Chunks'%opt.output)
    else:
        os.system('eos mkdir %s'%opt.output)
        os.system('eos mkdir %s/Chunks'%opt.output)
   

    #prepare jobs
    FarmDirectory = '%s/%s'%(cmsswBase,opt.farmappendix)
    os.system('mkdir -p %s'%FarmDirectory)
    print 'Preparing %d tasks to submit to the batch'%opt.njobs
    print 'Executables and condor wrapper are stored in %s'%FarmDirectory
    with open ('%s/condor.sub'%FarmDirectory,'w') as condor:
        
        condor.write('executable = {0}/$(cfgFile).sh\n'.format(FarmDirectory))
        #condor.write('output     = {0}/output_$(cfgFile).out\n'.format(FarmDirectory))
        #condor.write('error      = {0}/output_$(cfgFile).err\n'.format(FarmDirectory))
        #condor.write('log        = {0}/output_$(cfgFile).log\n'.format(FarmDirectory))
        condor.write('+JobFlavour = "{0}"\n'.format(opt.queue))

        for i in xrange(0,opt.njobs):

            cfgFile='job_%d'%i

            condor.write('cfgFile=%s\n'%cfgFile)
            condor.write('queue 1\n')
            condor.write('max_retries = 2\n')
                
            with open('%s/%s.sh'%(FarmDirectory,cfgFile),'w') as cfg:

                cfg.write('#!/bin/bash\n')
                cfg.write('trap "exit" INT\n')
                cfg.write('WORKDIR=`pwd`\n')
                cfg.write('echo "Working directory is ${WORKDIR}"\n')
                cfg.write('cd %s\n'%cmsswBase)
                cfg.write('eval `scram r -sh`\n')
                cfg.write('cd ${WORKDIR}\n')
                cfg.write('cmsRun %s/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzerGEN_cfg.py %s initialSeed=%d\n'%(cmsswBase,opt.cmsswOpts,i))
                for x in opt.outList.split(','):
                    base,ext=os.path.splitext(x)
                    cfg.write('xrdcp ${WORKDIR}/%s %s \n'%(x,os.path.join('root://eoscms//eos/cms/',opt.output,base+'_%d'%i+ext)))
                    cfg.write('rm ${WORKDIR}/%s'%x)

                os.system('chmod u+x %s/%s.sh'%(FarmDirectory,cfgFile))

        print 'Submitting jobs to condor, flavour "%s"'%(opt.queue)
        #os.system('condor_submit %s/condor.sub'%FarmDirectory)
        

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
