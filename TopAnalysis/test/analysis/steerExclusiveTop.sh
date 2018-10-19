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
githash=f93b8d8
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
outdir=/store/cmst3/user/psilva/ExclusiveAna
samples_json=test/analysis/pps_samples.json
zx_samples_json=test/analysis/zx_samples.json
wwwdir=~/www/ExclusiveAna
inputfileTag=MC13TeV_2017_ZH
inputfileTESTSEL=/store/cmst3/group/top/RunIIReReco/f93b8d8/${inputfileTag}/MergedMiniEvents_0_ext0.root
lumi=41367
lumiUnc=0.025


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	python scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} --tag ${inputfileTag} \
            -o testsel.root --genWeights genweights_${githash}.root\
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
        ;;
    SEL )
        baseOpt="-i ${eosdir} --genWeights genweights_${githash}.root"
        baseOpt="${baseOpt} -o ${outdir} -q ${queue} --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts"
        baseOpt="${baseOpt} --exactonly"
	python scripts/runLocalAnalysis.py ${baseOpt} --only ${samples_json};
	python scripts/runLocalAnalysis.py ${baseOpt} --only ${zx_samples_json};
	;;
    MERGE )
	./scripts/mergeOutputs.py /eos/cms/${outdir} True;
	;;
    PLOT )
	commonOpts="-i /eos/cms/${outdir} --puNormSF puwgtctr -j ${samples_json} --signalJson ${zx_samples_json} -l ${lumi} --mcUnc ${lumiUnc} "
	python scripts/plotter.py ${commonOpts} -O plots; 
	;;
    WWW )
	mkdir -p ${wwwdir}
	cp plots/*.{png,pdf} ${wwwdir}
	cp test/index.php ${wwwdir}
	;;
esac
