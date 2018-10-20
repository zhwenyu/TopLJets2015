import os
import sys
import optparse
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',      dest='inDir',       help='input directory with files',  default=None,       type='string')
    parser.add_option('-o', '--outDir',     dest='outDir',      help='output directory with files', default=None,       type='string')
    parser.add_option(      '--only',       dest='only',        help='only this tag',               default=None    ,   type='string')
    parser.add_option(      '--farm',       dest='farm',        help='farm tag',               default=None    ,   type='string')
    parser.add_option('-q', '--queue',      dest='queue',       help='batch queue [%default]',                 default='workday',  type='string')
    (opt, args) = parser.parse_args()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir

    cmsswBase=os.environ['CMSSW_BASE']
    FARMDIR='%s/INTEGRITYFARM'%cmsswBase
    if opt.farm: FARMDIR += opt.farm
    os.system('mkdir -p %s'%FARMDIR)
    print 'Submissions scripts @ ',FARMDIR

    onlyList=opt.only.split(',') if opt.only else []

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')

    with open('%s/checkInteg.sh'%FARMDIR,'w') as shell:
        shell.write('#!/bin/bash\n')
        shell.write('WORKDIR=`pwd`\n')
        shell.write('echo "Working directory is ${WORKDIR}"\n')
        shell.write('cd %s\n'%cmsswBase)
        shell.write('eval `scram r -sh`\n')
        shell.write('cd ${WORKDIR}\n')
        shell.write('echo "CMSSW_BASE=${CMSSW_BASE}"\n')
        shell.write('echo "Preparing output directory"\n')
        shell.write('mkdir -p /eos/cms/%s\n'%opt.outDir)
        shell.write('echo "Running integrity checker"\n')
        shell.write('python ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/checkProductionIntegrity.py -i ${1} -o ${2} --nocheck --only ${3}\n')
    os.system('chmod u+x %s/checkInteg.sh'%FARMDIR)

    condor=open('%s/condor.sub'%FARMDIR,'w')
    condor.write('executable = %s/checkInteg.sh\n'%FARMDIR)
    condor.write('output = %s/job_$(ProcId).out\n'%FARMDIR)
    condor.write('error  = %s/job_$(ProcId).err\n'%FARMDIR)
    condor.write('log    = %s/job_$(ProcId).log\n'%FARMDIR)
    condor.write('+JobFlavour ="%s"\n'%opt.queue)

    for dset in dset_list:
        dsetname=dset.split('/')[-1]

        pub_list=getEOSlslist(directory=dset,prepend='')
        for pubDir in pub_list:

            if not 'crab' in pubDir:
                print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
                continue
            pub=pubDir.split('/crab_')[-1]

            if len(onlyList)>0:
                found=False
                for tag in onlyList:
                    if not tag in pub: continue
                    found=True
                    break
                if not found : 
                    print 'Skipping %s, not in process only list'%pub
                    continue

            condor.write('arguments = %s %s %s\n'%(opt.inDir,opt.outDir,pub))
            condor.write('queue 1\n')

    condor.close()
    os.system('condor_submit %s/condor.sub'%FARMDIR)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

