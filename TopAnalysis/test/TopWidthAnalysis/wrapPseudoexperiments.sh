#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
WORKDIR=${SCRIPTPATH}/../

#configure environment
cd /afs/cern.ch/user/${USER: 0 : 1}/${USER}/CMSSW_7_4_7/src/ 
#export SCRAM_ARCH=$ARCH
eval `scram r -sh`

#run with the arguments passed
cd $1
echo ${PWD}
$2

echo $2 
