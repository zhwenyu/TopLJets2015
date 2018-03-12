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

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=workday
githash=c29f431
eosdir=/store/cmst3/group/top/RunIIFall17/${githash}
lumi=41367
lumiUnc=0.025
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/ExclusiveTop
wwwdir=~/www/ExclusiveTop


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	python scripts/runLocalAnalysis.py \
            -i ${eosdir}/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root \
            -o MC13TeV_TTJets.root \
            --njobs 1 -q local \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
        ;;
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            --only data/era2017/top_samples.json --exactonly \
            -o ${outdir} \
            -q ${queue} \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
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
