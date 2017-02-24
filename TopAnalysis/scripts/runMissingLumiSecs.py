#!/usr/bin/env python

import das_client
from sys import argv
import os
import json
def getListOfFiles(inputPD,runSel):
      das_query_string = 'das_client.py --query="file dataset=%s  run in [%s] | grep file.name" --limit=0' % (inputPD,runSel)
      das_output = os.popen(das_query_string).read() #query DAS via bash
      das_output = das_output.replace("\"\"","'")
      das_output = das_output.replace("root","root',")
      das_output = das_output.replace("\n","") #don't want newline in string
      das_output = das_output[:-1]#remove last comma in list of secondary files
      return das_output

job=argv[1]

print 'Starting',job

jsonF='%s/results/lumisToProcess.json'%job
runSel=[]
try:
      with open(jsonF) as missingLumis :
            data = json.load(missingLumis)
            runSel=[int(x) for x in data.keys()]
except:
      print 'No missing lumis - yay!'
if len(runSel)==0 : exit
print '\t %d runs with missing lumi sections'%len(runSel)

cfg=job.replace('crab_','')+'_cfg.py'
newCfg=cfg.replace('_cfg','_ext_cfg')
newCfgFile=open(newCfg,'w')
for line in open(cfg,'r'):
      newLine=line      
      if 'requestName' in newLine : newLine=newLine[:-2]+"_ext\"\n"
      if 'lumiMask' in newLine    : newLine='config.Data.lumiMask = \"%s\"\n'%os.path.abspath(jsonF)
      if 'unitsPerJob' in newLine : newLine='config.Data.unitsPerJob = 25\n'
      newCfgFile.write( newLine )
newCfgFile.close()
print '\t new cfg to process missing lumis @',newCfg
                 
os.system('crab submit -c %s' % newCfg)
print '\t jobs have been submitted'
