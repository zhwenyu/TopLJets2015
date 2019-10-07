import ROOT
import os
import sys
import itertools

KINEMATICS = [('bosonpt>20', 'bosonpt>95'),
              ('bosonpt>30', 'bosonpt>100'),
              ('bosonpt>40', 'bosonpt>110'),
              ('bosonpt>50', 'bosonpt>120'),
              ('bosonpt>60', 'bosonpt>130')]
RPSEL      = ['csi1>0.05  && csi2>0.05',
              'csi1>0.055 && csi2>0.055',
              'csi1>0.06  && csi2>0.06']
CATEGS     = ['nvtx>0',
              'nvtx<15,nvtx>=15',
              'nvtx<20,nvtx>=20',
              'nvtx<25,nvtx>=25']


baseDir=sys.argv[1]

ipt=0
for ana in list(itertools.product(KINEMATICS, RPSEL,CATEGS)):
    ipt+=1

    kin,rpsel,cats=ana

    workDir='%s/optim_%d'%(baseDir,ipt)
    print workDir

    os.system('mkdir -p ' + workDir)
    with open('%s/optimJob.sh'%workDir,'w') as script:
        script.write('#!/bin/bash\n')

        #initialization
        script.write('\n')
        script.write('input=/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p05/\n')
        script.write('preselZ="%s && %s"\n'%(kin[0],rpsel))
        script.write('preselGamma="%s && %s"\n'%(kin[1],rpsel))        
        script.write('categs="%s"\n'%(cats))
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
        script.write('python ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/generateBinnedWorkspace.py -i ${input} -o ${output} --preselZ "${preselZ}" --preselGamma "${preselGamma}" --categs "${categs}"\n')
        script.write('\n')

        #combine cards
        script.write('echo "Combining datacards"\n')
        script.write('cd $output\n')        
        script.write('combineCards.py zmm120=shapes-parametric.datacard_zmm_a120.dat zmm130=shapes-parametric.datacard_zmm_a130.dat zmm140=shapes-parametric.datacard_zmm_a140.dat zmm150=shapes-parametric.datacard_zmm_a150.dat zmm120=shapes-parametric.datacard_zmm_a120.dat zee130=shapes-parametric.datacard_zee_a130.dat zee140=shapes-parametric.datacard_zee_a140.dat zee150=shapes-parametric.datacard_zee_a150.dat > z_datacard.dat\n')
        script.write('combineCards.py g120=shapes-parametric.datacard_g_a120.dat g130=shapes-parametric.datacard_g_a130.dat g140=shapes-parametric.datacard_g_a140.dat g150=shapes-parametric.datacard_g_a150.dat > g_datacard.dat\n')
        script.write('\n')
        
        #run statistical analysis
        script.write('echo "Running combine"\n')
        script.write('for m in 800 900 1000 1080 1200 1320 1400 1500; do\n')
        script.write('\t for b in Z g; do\n')
        script.write('\t\t text2workspace.py ${b}_datacard.dat -o ${b}_workspace.root\n')
        script.write('\t\t combine ${b}_workspace.root -n PP${b}X -M AsymptoticLimits -m ${m}\n')
        script.write('\t\t combine ${b}_workspace.root -n PP${b}X -M Significance     -m ${m} --expectSignal=1 --setParameters mu_outfidsig=1\n')
        script.write('\t\t combine ${b}_workspace.root -n PP${b}X -M FitDiagnostics   -m ${m} --expectSignal=1 --setParameters mu_outfidsig=1\n')
        script.write('\t done\n')
        script.write('done\n')
        script.write('cd -\n')        

#submit optimization points to crab
print 'Will submit %d optimization scan points'%ipt
with open('zxstatana_scan.sub','w') as condor:
    condor.write("executable  = %s/optim_$(point)/optimJob.sh\n"%os.path.abspath(baseDir))
    condor.write("output       = zxstatana_scan.out\n")
    condor.write("error        = zxstatana_scan.err\n")
    condor.write("log          = zxstatana_scan.log\n")
    condor.write("+JobFlavour = \"tomorrow\"\n")
    condor.write("request_cpus = 4\n")
    for i in range(ipt):
        condor.write("point=%d\n"%(i+1))
        condor.write("queue 1\n")
os.system('condor_submit zxstatana_scan.sub')
