## command python ....py path to fit_type_dir

import sys, os

tag_list = [
#'data',
 'tbart',
# 'tbart0p5w',
# 'tbart169p5',
# 'tbart175p5',
# 'tbart2l2nu',
# 'tbart2l2nufxfx',
# 'tbart4w',
# 'tbarterdon',
# 'tbartfsrup',
# 'tbarthdampdn',
# 'tbarthdampup',
# 'tbartisrdn',
# 'tbartuedn'
]

job = 'fitdil_inc_condor3.sub'
script=open('/afs/cern.ch/user/w/wenyu/afswork/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/%s'%(job), 'w+')
script.write('executable  = /afs/cern.ch/work/w/wenyu/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/runFitWrapper.sh \n')
script.write('output      = %s.out \n'%(job))
script.write('error       = %s.err \n'%(job))
script.write('log         = %s.log \n'%(job))
script.write('+JobFlavour = "workday" \n')

num = 0
baseDir=sys.argv[1]
for dist in os.listdir(baseDir): # scenarios dir
    if '.pdf' in dist: continue
    elif '.png' in dist: continue
    distDir=os.path.join(baseDir,dist)
    for tag in tag_list :
        rootname = 'fitresults_'+ tag + '.root'
        if os.path.isfile(os.path.join(distDir,rootname)): continue
	else:
	    script.write('arguments   = /afs/cern.ch/user/w/wenyu/afswork/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/%s/runFit_%s.sh \n'%(distDir, tag))
	    script.write('queue 1 \n')
	    num = num +1

script.close()
print "jobs to resubmit = %s jobs"%(num)


