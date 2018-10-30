#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerExclusiveTop.sh <SEL/MERGE/PLOT/WWW>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output"
    echo "        PLOT         - make plots"
    echo "        WWW          - move plots to web-based are"
    echo "        PREPAREANA   - prepare bank of events for event mixing from the summary trees"
    echo "        ANA          - run analysis on the summary trees"
    echo "        PLOTANA      - plot analysis results"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=tomorrow
githash=f93b8d8
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
outdir=/store/cmst3/user/psilva/ExclusiveAna
samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/samples.json
jetht_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/jetht_samples.json
zx_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
wwwdir=~/www/ExclusiveAna
inputfileTag=MC13TeV_2017_ZH
inputfileTESTSEL=/store/cmst3/group/top/RunIIReReco/f93b8d8/${inputfileTag}/MergedMiniEvents_0_ext0.root
#inputfileTag=Data13TeV_2017B_DoubleMuon
#inputfileTESTSEL=/store/cmst3/group/top/RunIIReReco/f93b8d8/${inputfileTag}/MergedMiniEvents_0_ext0.root
lumi=41367
ppsLumi=37500
lumiUnc=0.025


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )

        #-i /store/cmst3/group/top/psilva/1e783f4/JetHT/crab_Data13TeV_2017C_JetHT/181018_194446/0000/MiniEvents_50.root \
        #-i /store/cmst3/group/top/psilva/1e783f4/DoubleEG/crab_Data13TeV_2017E_DoubleEG/181018_194206/0000/MiniEvents_7.root 
        #
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} --tag ${inputfileTag} \
            -o testsel.root --genWeights genweights_${githash}.root \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
        ;;
    SEL )
        baseOpt="-i ${eosdir} --genWeights genweights_${githash}.root"
        baseOpt="${baseOpt} -o ${outdir} -q ${queue} --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts"
        baseOpt="${baseOpt} --exactonly"
	python scripts/runLocalAnalysis.py ${baseOpt} --only ${samples_json};
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} --only ${zx_samples_json};
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} --only ${jetht_samples_json};
	;;
    MERGE )
	./scripts/mergeOutputs.py /eos/cms/${outdir} True;
	;;
    PLOT )
	commonOpts="-i /eos/cms/${outdir} --puNormSF puwgtctr -j ${samples_json} -l ${lumi} --mcUnc ${lumiUnc} "
	python scripts/plotter.py ${commonOpts} -O plots; 
	;;
    WWW )
	mkdir -p ${wwwdir}/presel
	cp plots/*.{png,pdf} ${wwwdir}/presel
	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel
	;;

    PREPAREANA )
        python runExclusiveAnalysis.py --step 0;        
        ;;

    ANA )
        python runExclusiveAnalysis.py --step 1;
        ;;

    PLOTANA )
        commonOpts="-i analysis -j ${samples_json} -l ${ppsLumi} --mcUnc ${lumiUnc} "
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} -O plots/ana;
        python compareBackgroundEstimation.py;
        mv *.{png,pdf} plots/ana/;
        ;;

    WWWANA )
        mkdir -p ${wwwdir}/ana
        cp plots/ana/*.{png,pdf} ${wwwdir}/ana
        cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/ana
        ;;

esac
