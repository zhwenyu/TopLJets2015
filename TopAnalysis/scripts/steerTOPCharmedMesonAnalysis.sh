#!/bin/bash

WHAT=$1; 
ERA=$2
if [ "$#" -ne 2 ]; then 
    echo "steerTOPCharmedMesonAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=2nd
githash=8db9ad6
lumi=12870
lumiUnc=0.062
eosdir=/store/cmst3/user/psilva/LJets2016/${githash}
summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}_ichepv2
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84
	lumiUnc=0.027
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}
	;;
esac

outdir=~/work/TopCharmedMesons_${ERA}
wwwdir=~/www/TopCharmedMesons_${ERA}

RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era ${ERA} -m TopCharmedMesonAnalysis::RunTopCharmedMesonAnalysis --ch 0 --only Data;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir};	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog --mcUnc ${lumiUnc};	
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
esac
