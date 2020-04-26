import ROOT
import os
import sys
import argparse
import itertools
from generateBinnedWorkspace import VALIDLHCXANGLES

KINEMATICS = '(((cat==169 || cat==121) && l1pt>30 && l2pt>20 && bosonpt>40) || (cat==22 && bosonpt>95))'
RPSEL      = 'csi1>0.035 && csi2>0.035'
CATEGS     = [list(itertools.product( ['protonCat==%d'%protonCat],
                                      ['']+['xangle==%d'%x for x in VALIDLHCXANGLES],
                                      ['']+['nvtx<20','nvtx>=20','nch<15','njets==0']))
              for protonCat in [1,2,3,4] ]
CATEGS= [filter(None,list(y)) for x in CATEGS for y in x]
OPTIMLIST=[ ' && '.join( [KINEMATICS]+[RPSEL]+x ) for x in CATEGS]

def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p05',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--injectMass',
                        dest='injectMass',
                        default=None,
                        help='mass to inject in pseudo-data [%default]')
    parser.add_argument('--injectMu',
                        dest='injectMu',
                        default=1.0,
                        help='signal strength of the mass to inject in pseudo-data [%default]')
    parser.add_argument('--signed',
                        dest='signed',
                        default=False,
                        help='used rapidity-signed missing mass [%default]',
                        action='store_true')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='ppvx_analysis',
                        help='Output directory [default: %default]')
    parser.add_argument('--just',
                        dest='just',
                        default=None,
                        help='Run only this optimization points (CSV list) [default: %default]')
    parser.add_argument('--unblind',
                        dest='unblind', 
                        default=False,
                        action='store_true',
                        help='Use non-mixed data in the final fit [default: %default]')
    opt=parser.parse_args(args)

    #build a list of the points to run
    if opt.just: opt.just=[int(x) for x in opt.just.split(',')]

    for ipt,ana in enumerate(OPTIMLIST):

        if opt.just and not ipt in opt.just: continue


        workDir='%s/optim_%d'%(opt.output,ipt)
        print workDir

        os.system('mkdir -p ' + workDir)
        with open('%s/optimJob.sh'%workDir,'w') as script:
            script.write('#!/bin/bash\n')

            #initialization
            script.write('\n')
            script.write('input=%s\n'%opt.input)
            script.write('cuts="%s"\n'%ana)
            if opt.injectMass:
                script.write('injectMassOpt="--injectMass %s"\n'%opt.injectMass)
            else:
                script.write('injectMassOpt=""\n')
            if opt.signed:
                script.write('extraOpt=--signed\n')
            else:
                script.write('extraOpt=""\n')
            if opt.unblind:
                print 'You\'re submitting unblinded jobs - how sure are you of what you did?'
                script.write('extraOpt="${extraOpt} --unblind"\n')

            script.write('output=%s\n'%os.path.abspath(workDir))     
            script.write('\n')

            #environment
            script.write('echo "Setting up environment"\n')
            script.write('cd %s/src\n'%os.environ['CMSSW_BASE'])
            script.write('eval `scram r -sh`\n')
            script.write('cd -\n')
            script.write('\n')

            #create datacard
            script.write('echo "Running datacard creation"\n')
            script.write('python ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/generateBinnedWorkspace.py -i ${input} -o ${output} --cuts "${cuts}" ${injectMassOpt} ${extraOpt}\n')
            script.write('\n')

            #combine cards
            script.write('echo "Combining datacards"\n')
            script.write('cd $output\n')        
            script.write('for v in z zmm zee g; do\n')
            script.write('\ta=(`ls shapes-parametric.datacard_${v}*.dat`)\n')
            script.write('\tcombstr=""\n')
            script.write('\tncards=${#a[@]}\n')
            script.write('\tfor (( i=0; i<${ncards}; i++ ));\n')
            script.write('\tdo\n')
            script.write('\t\tcombstr="${combstr} ${v}cat${i}=${a[$i]}"\n')
            script.write('\tdone\n') 
            script.write('\tcombineCards.py $combstr > ${v}_datacard.dat\n')
            script.write('done\n')
            script.write('\n')

    #submit optimization points to crab
    print 'Will submit %d optimization scan points'%len(OPTIMLIST)
    with open('%s/zxstatana_scan.sub'%opt.output,'w') as condor:
        condor.write("executable  = %s/optim_$(ProcId)/optimJob.sh\n"%os.path.abspath(opt.output))
        condor.write("output       = zxstatana_scan.out\n")
        condor.write("error        = zxstatana_scan.err\n")
        condor.write("log          = zxstatana_scan.log\n")
        condor.write("+JobFlavour = \"tomorrow\"\n")
        condor.write("request_cpus = 4\n")
        condor.write("queue %d\n"%len(OPTIMLIST))
    os.system('condor_submit %s/zxstatana_scan.sub'%opt.output)
    

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
