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
githash=5fb8f4f
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
lumi=41367
lumiUnc=0.025
outdir=/eos/cms/store/cmst3/user/psilva/ExclusiveAna
wwwdir=~/www/ExclusiveTop
inputfileTESTSEL=/store/cmst3/group/top/RunIIReReco/5fb8f4f/Data13TeV_2017D_DoubleMuon/MergedMiniEvents_0_ext0.root

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	python scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} -o testsel.root \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
        ;;
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} --genWeights genweights_${githash}.root \
            -o ${outdir} -q ${queue} \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
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
