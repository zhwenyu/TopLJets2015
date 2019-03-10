#!/bin/bash
HOME=`pwd`
CMSSW=/afs/cern.ch/user/p/psilva/work/CMSSW_10_3_0_pre4/
PPS=/afs/cern.ch/user/p/psilva/work/CMSSW_9_4_11/src/TopLJets2015/TopAnalysis/test/analysis/pps
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
python $PPS/trainPUdiscriminators.py -m $PPS/pu_models.pck -o ${1} -i ${2}

