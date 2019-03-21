import pickle
import os
import sys
from random import shuffle
from collections import defaultdict

def runLocalPacked(args):
    baseDir,f=args
    os.system('sh test/analysis/pps/wrapAnalysis.sh 0 /tmp/ {0}/Chunks {1}'.format(baseDir,f.replace('.pck','.root')))
    os.system('mv -v /tmp/Chunks/{1} {0}/mixing/Chunks/{1}'.format(baseDir,f))

rpData = defaultdict(list)

baseDir=sys.argv[1]
fList=[f.replace('.root','.pck') for f in os.listdir(baseDir+'/Chunks') if 'Data' in f]

toCheck=[]
for f in fList:

    if 'ZeroBias' in f : continue
    if 'DoubleEG' in f : continue
    if 'Photon' in f : continue

    fullf=os.path.join(baseDir+'/mixing/Chunks',f)

    #check if file was produced
    if not os.path.isfile(fullf):
        toCheck.append(f)
        continue

    try:
        with open(fullf,'r') as cache:
            a=pickle.load(cache)
        for key in a: 
            newKey=(key[0],int(key[1]),int(key[2]))
            rpData[newKey] += a[key]
    except Exception as e:
        toCheck.append(f)

print 'Total number of events for mixing'
for key in rpData:
    print key,len(rpData[key])

with open('mixbank.pck','w') as cache:
    pickle.dump(rpData,cache)
os.system('mv -v mixbank.pck %s/mixing/mixbank.pck'%baseDir)

if len(toCheck)>0:
    print 'Re-running locally',len(toCheck),'jobs'
    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map( runLocalPacked,[ (baseDir,x) for x in toCheck ])

