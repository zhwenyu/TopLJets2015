#!/usr/bin/env python

import sys
import os
import ROOT

def getEntries(url,tname):
    if not os.path.isfile(url): 
        raise NameError('No file for %s'%tname)

    f=ROOT.TFile.Open(url)
    if f.IsZombie():
        raise RuntimeError('Zombie file for %s'%tname)
    n=f.Get(tname).GetEntriesFast()
    f.Close()

    return n

def runPredPacked(args):
    url,tag=args
    os.system('sh test/analysis/pps/wrapPUDiscrTrain.sh {0}/pudiscr {1}'.format(url,tag))

runLocally=False
dontrun=False
url=sys.argv[1]

toCheck=[]
for f in os.listdir(url):
    fullf=os.path.join(url,f)
    if not os.path.isfile(fullf) : 
        continue
    try:
        norig=getEntries(fullf,'tree')
        nfriend=getEntries(os.path.join(url,'pudiscr',f),'pudiscr')        
        if norig!=nfriend : raise ValueError('%d!=%d'%(norig,nfriend))
    except Exception as e:
        toCheck.append( (fullf,e) )

if dontrun:

    for x in toCheck:
        print x

else:

    if runLocally:
        import multiprocessing as MP
        pool = MP.Pool(8)
        task_list=[]
        for x in toCheck:
            task_list.append( (url,x[0]) )
        pool.map( runPredPacked,task_list)

    else:
        with open('runpred_condor.sub','w') as c:
            c.write("executable  = {0}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapPUDiscrTrain.sh\n".format(os.environ['CMSSW_BASE']))
            c.write("output      = runpred_condor.out\n")
            c.write("error       = runpred_condor.err\n")
            c.write("log         = runpred_condor.log\n")
            c.write("arguments   = {0}/pudiscr $(chunk)\n".format(url))
            for x in toCheck: 
                c.write("chunk={0}\n".format(x[0]))
                c.write("queue 1\n")
            os.system('condor_submit runpred_condor.sub')
