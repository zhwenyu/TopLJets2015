#!/bin/bash
HOME=`pwd`
CMSSW=/afs/cern.ch/user/p/psilva/work/CMSSW_10_3_0_pre4/
PY=/afs/cern.ch/user/p/psilva/work/CMSSW_9_4_11/src/TopLJets2015/TopAnalysis/test/analysis/pps
PPS=/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/train_results/
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
python $PY/trainPUdiscriminators.py -m $PPS/pu_models_Zmm.pck -o ${1} -i ${2}

