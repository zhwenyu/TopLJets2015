#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerTOPUE.sh <TESTSEL/SEL/MERGE/PLOTSEL/WWWSEL>";
    echo "        TESTSEL      - test locally the selection script on a single file";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
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
summaryeosdir=/store/cmst3/group/top/summer2017/WbChargeAsymm
outdir=test/summer2017/WbChargeAsymm
wwwdir=~/www/WbChargeAsymm


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_W1Jets/MergedMiniEvents_0.root
	analysisWrapper \
	    --in ${file} \
	    --out w1jets_test.root \
	    --era ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016 \
	    --method WbChargeAsymmetry::RunWbChargeAsymmetry \
	    --ch 0;
	;;

    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            --only test/summer2017/wb_samples.json --exactonly \
            -q ${queue} -o ${summaryeosdir} \
            --era era2016 -m WbChargeAsymmetry::RunWbChargeAsymmetry --ch 0 --runSysts;
	;;

    MERGE )
	mkdir -p ${outdir};
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount eos;
	./scripts/mergeOutputs.py eos/cms${summaryeosdir} True ${outdir};	
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse umount eos;
	;;
    PLOT )
	commonOpts="-i ${outdir} --puNormSF puwgtctr -j test/summer2017/wb_samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts}; 
	;;
    WWW )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
esac
