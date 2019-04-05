#!/usr/bin/env python
import os
import sys

def runPacked(args):
    tag,cmd=args
    print 'Creating cards for',tag
    os.system(cmd)


cmd_list=[]
baseDir=sys.argv[1]
motherDir=os.path.dirname(baseDir)
toCheck=['tbart.datacard.dat','pseudodata.tbart.shapes.root']
base_anchor='{0},%s/{1}/MC13TeV_2016_TTJets.root'%motherDir
base_cmd='test/analysis/top17010/prepareDataCard.py -d {0} -t %s/templates dataDef=sig,%s/plots/plotter.root,{0}/{0}_t#bar{{t}} dataDef=data,%s/plots/plotter.root,{0}/{0}  -s {1} -o %s/datacards --systs test/analysis/top17010/systs_dict.json'%(motherDir,motherDir,motherDir,motherDir)



for dist in os.listdir(baseDir):
    distDir=os.path.join(baseDir,dist)
    for a in os.listdir(distDir):
        allFound=True
        for f in toCheck:
            if os.path.isfile(os.path.join(distDir,a,f)): continue
            allFound=False
            break
        if allFound: continue

        anchor=base_anchor.format(a,a if a!='nom' else '')
        cmd_list.append( (dist+a,base_cmd.format(dist+'_mlb',anchor)) )

print 'Print  i have found %d directories which need to re-run'%len(cmd_list)
print cmd_list
from multiprocessing import Pool
pool = Pool(8)
pool.map(runPacked,cmd_list)
        
