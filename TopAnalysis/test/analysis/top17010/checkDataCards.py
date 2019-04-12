#!/usr/bin/env python
import os
import sys

#toCheck=['tbart.datacard.dat','pseudodata.tbart.shapes.root', 'ttbar.shapes.root']
#base_anchor='{0},/eos/cms/store/cmst3/group/top/TOP17010/0c522df/{1}/MC13TeV_2016_TTJets.root'
#base_cmd='${{CMSSW_BASE}}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/prepareDataCard.py -d {0} -t /eos/cms/store/cmst3/group/top/TOP17010/0c522df/templates dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}4w dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}0.5w dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}hdampup dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}hdampdn dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}uedn dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}175.5 dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}169.5 dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}erdon dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}isrdn dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}fsrup dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}2l2nufxfx dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}2l2nu dataDef=sig,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/plotter.root,{0}/{0}_t#bar{{t}} dataDef=data,/eos/cms/store/cmst3/group/top/TOP17010/0c522df/plots/plotter.root,{0}/{0}  -s {1} -o /eos/cms/store/cmst3/group/top/TOP17010/0c522df/datacards --systs /afs/cern.ch/user/w/wenyu/afswork/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/systs_dict.json'


def runPacked(args):
    tag,cmd=args
    print 'Creating cards for',tag
    print cmd
    os.system(cmd)


cmd_list=[]
baseDir=sys.argv[1]  # /eos/cms/store/cmst3/group/top/TOP17010/0c522df/datacards
motherDir=os.path.dirname(baseDir)
tagList=['','0p5w','4w','169p5','175p5','fsrup','fsrdn']
toCheck=['tbart%s.datacard.dat'%tag for tag in tagList]
toCheck+=[ 'pseudodata.tbart%s.shapes.root'%tag for tag in tagList]
base_anchor='{0},%s/{1}/MC13TeV_2016_TTJets.root'%motherDir
base_cmd='test/analysis/top17010/prepareDataCard.py -d {0} -t %s/templates dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}4w dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}0.5w dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}175.5 dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}169.5 dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}fsrup dataDef=sig,%s/plots/syst_plotter.root,{0}/{0}_t#bar{{t}}fsrdn dataDef=data,%s/plots/plotter.root,{0}/{0} -s {1} -o %s/datacards --systs test/analysis/top17010/systs_dict.json'%(motherDir,motherDir,motherDir,motherDir,motherDir,motherDir,motherDir,motherDir,motherDir)

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
        if not a in d_anchors:
            allFound=False
        else:
            for f in toCheck:
                if ( os.path.isfile(os.path.join(distDir,a,f)) and os.path.getsize(os.path.join(distDir,a,f)) >0 ): continue  # check exist and size
                allFound=False
                break
        if allFound: continue

        anchor=base_anchor.format(a,a if a!='nom' else '')
        cmd_list.append( (dist+a,base_cmd.format(dist+'_mlb',anchor)) )

print 'Print  i have found %d directories which need to re-run'%len(cmd_list)
print cmd_list
#exit(0)

from multiprocessing import Pool
pool = Pool(8)
pool.map(runPacked,cmd_list)
        
