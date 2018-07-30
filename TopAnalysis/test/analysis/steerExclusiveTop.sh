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
queue=tomorrow
githash=fbc74ae
eosdir=/store/cmst3/group/top/RunIIFall17/${githash}
lumi=41367
lumiUnc=0.025
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/ExclusiveTop
wwwdir=~/www/ExclusiveTop
inputfileTESTSEL=/store/cmst3/group/top/RunIIFall17/fbc74ae/Data13TeV_SingleMuon_2017D/MergedMiniEvents_2_ext0.root

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	python scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} \
            -o TESTSEL.root \
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
    SELDATA )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
	    --only data/era2017/top_samples_Dataonly.json --exactonly \
            -o ${outdir} -q ${queue} \
	    --skipexisting \
            --era era2017 \
	    -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;	
	;;
    SELMC )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts \
	    -o ${outdir} \
            -q ${queue} \
            --skipexisting \
            --only data/era2017/top_samples_MConly.json --exactonly;
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir};
	;;
    PLOT )
	commonOpts="-i ${outdir} --puNormSF puwgtctr -j data/era2017/top_samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts}; 
	;;
    WWW )
	mkdir -p ${wwwdir}
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}
	cp test/index.php ${wwwdir}
	;;
esac
