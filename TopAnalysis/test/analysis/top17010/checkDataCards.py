#!/usr/bin/env python
import os
import sys

def runPacked(args):
    tag,cmd=args
    print 'Creating cards for',tag
    print cmd
    os.system(cmd)


cmd_list=[]
baseDir=sys.argv[1]
motherDir=os.path.dirname(baseDir)
tagList=['','fsrup','fsrdn','erdon','gluonmove','qcdBased','scenario1179712','scenario131125','scenario1704005','scenario655418','scenario917534','scenario917634']
toCheck=['tbart%s.datacard.dat'%tag for tag in tagList]+['data.datacard.dat']
toCheck+=[ 'pseudodata.tbart%s.shapes.root'%tag for tag in tagList]+['data.shapes.root']
base_anchor='{0},%s/{1}/MC13TeV_2016_TTJets.root'%motherDir
#base_cmd='test/analysis/top17010/prepareDataCard.py -d {0} -t %s/templates dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}4w dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}0.5w dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}175.5 dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}169.5 dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}fsrup dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}fsrdn dataDef=data,%s/plots/plotter.root,{0}/{0} -s {1} -o %s/datacards --systs test/analysis/top17010/systs_dict.json'%(motherDir,motherDir,motherDir,motherDir,motherDir,motherDir,motherDir,motherDir,motherDir)

dataDefList='dataDef=sig,{dir}/syst_plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}}gluonmove dataDef=sig,{dir}/syst_plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}}qcdBased dataDef=sig,{dir}/syst_plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}}erdon dataDef=sig,{dir}/syst_plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}}fsrup dataDef=sig,{dir}/syst_plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}}fsrdn dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}},scenario1704005/MC13TeV_2016_TTJets.root,{{dist}} dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}},scenario1179712/MC13TeV_2016_TTJets.root,{{dist}} dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}},scenario917634/MC13TeV_2016_TTJets.root,{{dist}} dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}},scenario917534/MC13TeV_2016_TTJets.root,{{dist}} dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}},scenario655418/MC13TeV_2016_TTJets.root,{{dist}} dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}},scenario131125/MC13TeV_2016_TTJets.root,{{dist}} dataDef=sig,{dir}/plotter.root,{{dist}}/{{dist}}_t#bar{{{{t}}}} dataDef=data,{dir}/plotter.root,{{dist}}/{{dist}}'.format(dir=motherDir+'/plots')
base_cmd='test/analysis/top17010/prepareDataCard.py -d {dist} -t %s/templates %s -s {anchor} -o %s/datacards --systs test/analysis/top17010/systs_dict.json'%(motherDir,dataDefList,motherDir)

anchors=[]
for dist in os.listdir(baseDir):
    distDir=os.path.join(baseDir,dist)
    for a in os.listdir(distDir):
        anchors.append(a)
anchors=list(set(anchors))
print '%s anchors found'%len(anchors)

for dist in os.listdir(baseDir):
    distDir=os.path.join(baseDir,dist)
    d_anchors=os.listdir(distDir)
    for a in anchors:

        allFound=True
        error=0
        if not a in d_anchors:
            allFound=False
            error=(1,a)
        else:
            for f in toCheck:
                if os.path.isfile(os.path.join(distDir,a,f)): continue
                allFound=False
                error=(1,f)
                break
        if allFound: continue
        print len(cmd_list),dist+a,error
        anchor=base_anchor.format(a,a if a!='nom' else '')
        cmd_list.append( (dist+a,base_cmd.format(dist=dist+'_mlb',anchor=anchor)) )

print 'I have found %d directories which need to re-run'%len(cmd_list)
#print cmd_list
from multiprocessing import Pool
pool = Pool(8)
pool.map(runPacked,cmd_list)
        
