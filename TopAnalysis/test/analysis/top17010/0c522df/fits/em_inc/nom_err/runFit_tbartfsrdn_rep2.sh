#!/bin/bash

OUTNUM=$1
#setup environment
cd /afs/cern.ch/work/w/wenyu/work/topwidth/CMSSW_10_3_0_pre4
eval `scram r -sh`
source /afs/cern.ch/user/b/bendavid/work/cmspublic/pythonvenv/tensorflowfit_h5py/bin/activate
cd -

#make local copy of the datacard
cat /afs/cern.ch/work/w/wenyu/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/0c522df/datacards1/em/nom/tbartfsrdn.datacard.dat > datacard.dat
sed -i "s|pseudodata_999|pseudodata_${OUTNUM}|g" datacard.dat

#convert to HDF5 and run TF-based fits
text2hdf5.py datacard.dat
combinetf.py datacard.dat.hdf5 -o fitresults_tbartfsrdn_${OUTNUM}.root --doImpacts
cp -v fitresults*root /afs/cern.ch/user/w/wenyu/afswork/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/0c522df/fits/em_inc/nom_err
#python /afs/cern.ch/user/w/wenyu/afswork/work/topwidth/TopLJets2015/TopAnalysis/test/analysis/top17010/nuisanceScan.py datacard.dat
#cp -v fitresults_fixedgroups.pck /afs/cern.ch/user/w/wenyu/afswork/work/topwidth/CMSSW_9_4_10/src/TopLJets2015/TopAnalysis/test/analysis/top17010/0c522df/fits/em_inc/nom_err/fitresults_tbartfsrdn_fixedgroups.pck
