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
inputDir='eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p05'
if len(sys.argv)>2:
    inputDir=sys.argv[2]

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
        script.write('input=%s\n'%inputDir)
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
        script.write('for v in z g; do\n')
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
        
        #run statistical analysis
        script.write('echo "Running combine"\n')
        script.write('for b in Z g; do\n')
        script.write('\t for m in 800 900 1000 1080 1200 1320 1400 1500; do\n')
        script.write('\t\t pfix=${b}_m${m}\n')
        script.write('\t\t text2workspace.py ${b}_datacard.dat -m ${m} -o ${pfix}_workspace.root\n')
        script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X -M AsymptoticLimits -m ${m}\n')
        script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X -M Significance     -m ${m} --expectSignal=1 --setParameters mu_outfidsig=1\n')
        script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X -M FitDiagnostics   -m ${m} --expectSignal=1 --setParameters mu_outfidsig=1\n')
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
