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

dontrun=False
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

if dontrun:

    for x in toCheck:
        print x

else:

    print 'Submitting',len(toCheck),'/',nTot,'jobs'
    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map( runAnaPacked,[ (url,x) for x in toCheck ])

