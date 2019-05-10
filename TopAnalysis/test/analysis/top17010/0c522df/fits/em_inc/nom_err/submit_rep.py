
wrtF = open('submit_rep.sub', 'w')
wrtF.write('executable  = /afs/cern.ch/user/w/wenyu/afswork/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/0c522df/fits/em_inc/nom_err/runFit_tbartfsrup_rep.sh')
wrtF.write('\n')
wrtF.write('output      = submit_rep.sub.out ')
wrtF.write('\n')
wrtF.write('error       = submit_rep.sub.err ')
wrtF.write('\n')
wrtF.write('log         = submit_rep.sub.log ')
wrtF.write('\n')
wrtF.write('requirements = (OpSysAndVer =?= "SLCern6") ')
wrtF.write('\n')
wrtF.write('+JobFlavour = "workday" ')
wrtF.write('\n')
wrtF.write('arguments   = $(num) ')
wrtF.write('\n')
wrtF.write('queue num from ( ')
wrtF.write('\n')

for i in range(3, 999):
  wrtF.write(str(i))
  wrtF.write('\n')

wrtF.write(')')

wrtF.close()  
