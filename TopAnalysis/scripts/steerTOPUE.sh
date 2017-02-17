#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerTOPUE.sh <SEL/PLOTSEL/WWWSEL>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N


queue=2nw
githash=f016290
lumi=36460
lumiSpecs="" #--lumiSpecs EE:11391"
lumiUnc=0.026
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
summaryeosdir=/store/cmst3/group/top/winter2017/TopUE
outdir=/afs/cern.ch/work/${myletter}/${whoami}/TopUE_ReReco2016
wwwdir=~/www/TopUE_ReReco2016/


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	analysisWrapper --in MC13TeV_TTJets.root \
	    --out ue_test.root \
	    --era ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016 \
	    --method TOP-UE::RunTopUE \
	    --ch 0;
	;;

    FULLSEL ) 
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era era2016 -m TOP-UE::RunTopUE --ch 0 --only TTJ;
	;;

    SEL )
	samplesToProcess=(Double,MuonEG DY,Single,W,ZZ _TT)
	for s in ${samplesToProcess[@]}; do
	    
	    echo -e "${RED} Submitting ${s} ${NC}"
	    python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era era2016 -m TOP-UE::RunTopUE --ch 0 --runSysts --only ${s} --babySit;

	    echo -e "${RED} Merging ${s} ${NC}"
	    ./scripts/mergeOutputs.py ${outdir} True;	

	    echo -e "${RED} Moving to store ${s} ${NC}"
	    a=(`ls ${outdir}/Chunks/*.root`)
	    for i in ${a[@]}; do
		baseName=`basename ${i}`
		xrdcp ${i} root://eoscms//eos/cms/${summaryeosdir}/${baseName};
		rm ${i};
	    done
	done
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir} True;
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc};	
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	eosprefix=root://eoscms//eos/cms
	echo "Computing resolutions"
	python test/TopUEAnalysis/runUEanalysis.py -i ${eosprefix}/${summaryeosdir}/MC13TeV_TTJets_0.root --step 0 --ptThr 1.0,0.9;

	#echo "Defining analysis configuration"
	#python test/TopUEAnalysis/runUEanalysis.py -i ${eosprefix}/${summaryeosdir}/MC13TeV_TTJets_dilpowheg_0.root --step 1;
	#python test/TopUEAnalysis/runUEanalysis.py --step 3;
	;;

esac
