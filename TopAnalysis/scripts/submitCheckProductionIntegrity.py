import os
import sys
import optparse
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
from subprocess import Popen, PIPE

"""
steer the script
"""
def main():

    eos_cmd = '/afs/cern.ch/project/eos/installation/cms/bin/eos.select'

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',      dest='inDir',       help='input directory with files',               default=None,   type='string')
    parser.add_option('-o', '--outDir',     dest='outDir',      help='output directory with files',              default=None,   type='string')
    parser.add_option('-q', '--queue',      dest='queue',       help='batch queue',                              default='2nd',  type='string')
    (opt, args) = parser.parse_args()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir

    cmsswBase=os.environ['CMSSW_BASE']
    FARMDIR='%s/src/TopLJets2015/TopAnalysis/FARM'%cmsswBase
    os.system('mkdir -p %s'%FARMDIR)

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    for dset in dset_list:
        dsetname=dset.split('/')[-1]

        pub_list=getEOSlslist(directory=dset,prepend='')
        for pubDir in pub_list:

            if not 'crab' in pubDir:
                print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
                continue
            pub=pubDir.split('/crab_')[-1]
            
            cfgfile='%s/mergejob_%s.sh'%(FARMDIR,pub)
            with open(cfgfile,'w') as cfg:
                cfg.write('WORKDIR=`pwd`\n')
                cfg.write('echo "Working directory is ${WORKDIR}"\n')
                cfg.write('cd %s\n'%cmsswBase)
                cfg.write('eval `scram r -sh`\n')
                cfg.write('cd ${WORKDIR}\n')
                cfg.write('echo "Preparing output directory"\n')
                cfg.write('%s mkdir %s\n'%(eos_cmd,opt.outDir))
                cfg.write('echo "Running integrity checker"\n')
                cfg.write('python ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/checkProductionIntegrity.py -i %s -o %s --nocheck --mount --only %s'%(opt.inDir,opt.outDir,pub))
                
            os.system('chmod u+x %s'%cfgfile)
            cmd='bsub -q %s %s' % (opt.queue,cfgfile)
            os.system(cmd)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

