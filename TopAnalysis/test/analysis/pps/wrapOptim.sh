#!/bin/bash
HOME=`pwd`
CMSSW=/afs/cern.ch/user/p/psilva/work/CMSSW_9_4_11
PPS=${CMSSW}/src/TopLJets2015/TopAnalysis/test/analysis/pps
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
python $PPS/optimizeSR.py ${1} ${2} ${3}