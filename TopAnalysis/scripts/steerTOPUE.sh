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
githash=b32c02e
lumi=36460
lumiSpecs="" #--lumiSpecs EE:11391"
lumiUnc=0.027
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
summaryeosdir=/store/cmst3/group/top/winter2017/TopUE
outdir=./UEanalysis/
wwwdir=~/www/TopUE_ReReco2016/


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	analysisWrapper --in root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b32c02e/MC13TeV_TTJets/MergedMiniEvents_0.root \
	    --out ue_test.root \
	    --era ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016 \
	    --method TOP-UE::RunTopUE \
	    --ch 0;
	;;

    FULLSEL ) 
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${summaryeosdir} --era era2016 -m TOP-UE::RunTopUE --ch 0;
	;;

    MERGE )
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount eos;
	./scripts/mergeOutputs.py eos/cms${summaryeosdir} True ${outdir};
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse umount eos;
	;;
    PLOTSEL )
	commonOpts="-i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts} --only mll --outName mll_plotter.root;	
     	python scripts/runDYRinRout.py --in ${outdir}/plots/mll_plotter.root --categs ""  --out ${outdir}/plots/;
	python scripts/plotter.py ${commonOpts} --procSF DY:${outdir}/plots/.dyscalefactors.pck --only njets --rebin 4 --saveTeX --outName count_plotter.root;
	python scripts/plotter.py ${commonOpts} --procSF DY:${outdir}/plots/.dyscalefactors.pck;
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	eosprefix=root://eoscms//eos/cms
	echo "Computing resolutions"
	base="${eosprefix}/${summaryeosdir}/Chunks/MC13TeV_TTJets"
	#python test/TopUEAnalysis/runUEanalysis.py -i ${base}_0.root,${base}_1.root,${base}_2.root,${base}_3.root,${base}_4.root, --step 0 --ptThr 1.0,0.9;

	echo "Defining analysis configuration"
	python test/TopUEAnalysis/runUEanalysis.py --step 1;
	
	echo "Filling the histograms"
	#python test/TopUEAnalysis/runUEanalysis.py -i ${summaryeosdir}/Chunks --step 2 -q ${queue};
	;;
    MERGEANA )
	./scripts/mergeOutputs.py UEanalysis/analysis_0_0 True 
	;;
    PLOTANA )
	#python scripts/plotter.py -i UEanalysis/analysis_0_0 -j data/era2016/samples.json      -l ${lumi} --saveLog --mcUnc ${lumiUnc};	
	#python scripts/plotter.py -i UEanalysis/analysis_0_0 -j data/era2016/syst_samples.json -l ${lumi} --saveLog --mcUnc ${lumiUnc} --silent --outName syst_plotter.root;	
	python test/TopUEAnalysis/compareAtRecoLevel.py UEanalysis/analysis_0_0/plots/plotter.root UEanalysis/analysis_0_0/plots/syst_plotter.root
	#python test/TopUEAnalysis/UETools.py -o UEanalysis/analysis_0_0/plots/
	;;
    WWWANA )
	mkdir -p ${wwwdir}/rawana
        cp UEanalysis/analysis_0_0/plots/*.{png,pdf} ${wwwdir}/rawana
	cp  UEanalysis/analysis_0_0/plots/plotter/*.{png,pdf} ${wwwdir}/rawana
        cp test/index.php ${wwwdir}/rawana
	;;

esac
