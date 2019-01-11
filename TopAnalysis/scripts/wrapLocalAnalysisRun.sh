#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
WORKDIR=`dirname ${SCRIPTPATH}`

#configure environment
cd $WORKDIR
eval `scram r -sh`
cd -

#run with the arguments passed
echo $@
#$*
