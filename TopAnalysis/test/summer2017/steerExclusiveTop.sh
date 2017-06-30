#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerExclusiveTop.sh <SEL/MERGE/PLOT/WWW>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output"
    echo "        PLOT         - make plots"
    echo "        WWW          - move plots to web-based are"
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N
queue=2nd
githash=b312177
lumi=35922
lumiUnc=0.025
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/summer2017/ExclusiveTop
wwwdir=~/www/ExclusiveTop


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    SEL )
        #to run locally use "--njobs 8 -q local" instead of "-q condor"
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            --only test/summer2017/exctop_samples.json --exactonly \
            -o ${outdir} -q condor \
            --era era2016 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
	;;

    MERGE )
	./scripts/mergeOutputs.py ${outdir};
	;;
    PLOT )
	commonOpts="-i ${outdir} --puNormSF puwgtctr -j test/summer2017/exctop_samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc} --noStack"
	python scripts/plotter.py ${commonOpts}; 
	;;
    WWW )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
esac
