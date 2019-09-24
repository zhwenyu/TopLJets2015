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
    n=0
    #n=f.Get(tname).GetEntriesFast()
    f.Close()

    return n

def runAnaPacked(args):
    url,tag=args
    os.system('sh test/analysis/pps/wrapAnalysis.sh 1 {0}/analysis {0}/Chunks {1} {0}/mixing/mixbank.pck'.format(url,tag))

runMode=2
url=sys.argv[1]

toCheck=[]
nTot=0
for f in os.listdir(url+'/Chunks'):
    #if 'MC13TeV' in f : continue
    nTot+=1
    fullf=os.path.join(url+'/analysis/Chunks',f)
    if not os.path.isfile(fullf) : 
        print 'Missing in action',f
        toCheck.append(f)
    else:
        try:
            nentries=getEntries(fullf,'data')                    
        except Exception as e:
            print 'Corrupted or empty file for',f,e
            toCheck.append(f)

if runMode==0:
    for x in toCheck:
        print x
elif runMode==1:
    print 'Running locally',len(toCheck),'/',nTot,'jobs'
    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map( runAnaPacked,[ (url,x) for x in toCheck ])
elif runMode==2:
    print 'Submitting',len(toCheck),'/',nTot,'jobs'
    cmssw_base=os.environ['CMSSW_BASE']
    with open('zxana_recover.sub','w') as cache:
        cache.write("executable  = %s/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh\n"%cmssw_base)
        cache.write("output      = zxana_recover.out\n")
        cache.write("error       = zxana_recover.err\n")
        cache.write("log         = zxana_recover.log\n")
        for x in toCheck:
            cache.write("arguments   = 1 {predout} {predin} {filein} {mix_file}\n".format(predout=url+'/analysis',
                                                                                              predin=url+'/Chunks',
                                                                                              filein=x,
                                                                                              mix_file=url+'/mixing/mixingback.pck'))
            cache.write("queue 1\n")
        #os.system('condor_submit zxana_recover.sub')
