#!/bin/bash
HOME=`pwd`
CMSSW=/afs/cern.ch/user/p/psilva/work/CMSSW_9_4_11
PPS=${CMSSW}/src/TopLJets2015/TopAnalysis/test/analysis/pps
samples_json=${PPS}/samples.json
RPout_json=${PPS}/golden_noRP.json
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
python $PPS/runExclusiveAnalysis.py --step $1 --json ${samples_json} --RPout ${RPout_json} -o $2 -i $3 --only $4

