#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerTOP17010.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW/HYPOTEST> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        DY           - estimate DY scale factor from data";
    echo "        PLOT_P1      - runs the plotter tool on the analysis outputs";
    echo "        PLOT_P2      - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo "        HYPOTEST     - create the datacards, steering scripts for hypothesis testing and submit to batch";
    echo "        PLOTHYPOTEST - create summaries of the hypothesis tests";
    echo "        WWWHYPOTEST  - move summaries to afs-web-based area"
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N


queue=2nd
githash=b312177
lumi=35922
lumiSpecs="" #"--lumiSpecs EE:11391"
lumiUnc=0.025
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
summaryeosdir=/store/cmst3/group/top/TOP-17-010/
COMBINERELEASE=~/scratch0/CMSSW_7_4_7/src/
outdir=/afs/cern.ch/work/${myletter}/${whoami}/TOP-17-010/
wwwdir=~/www/TOP-17-010/


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    TEST )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q local -o /tmp/`whoami` --era era2016 -m TOP-17-010::RunTop17010 --ch 0 --runSysts --only TTJets;
        ;;
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${summaryeosdir} --era era2016 -m TOP-17-010::RunTop17010 --ch 0 --runSysts;
	;;
    MERGESEL )
	mkdir -p ${outdir}
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount eos;
	./scripts/mergeOutputs.py eos/cms${summaryeosdir} True ${outdir};	
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse umount eos;
	;;
    PLOTSEL )
        commonOpts="-i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts}
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${outdir}/analysis/Chunks -q ${queue};	
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    BKG )
	python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json  -l ${lumi} ${lumiSpecs} --onlyData --only mll -o dy_plotter.root; 
	;;
    DY )
	python scripts/runDYRinRout.py --in ${outdir}/analysis/plots/dy_plotter.root --categs 1b,2b --out ${outdir}/analysis/plots/;
	;;
    PLOT_P1 )
		#python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --only count --saveTeX -o count_plotter2.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck > count2.out & 
		#python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --only njets,ptlb -o njets_plotter2.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck > njets2.out &
    	#python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck -o plotter4.root > plotter4.out &
		#python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/syst_samples.json -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --silent -o syst_plotter2.root > syst2.out &
		#stdcmd="python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --only count --saveTeX -o count_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck;"
	    #	bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd}; 
		#stdcmd="python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --only njets,ptlb -o njets_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck; "
	    #	bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd}; 
    	#stdcmd="python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck &"
	    #	bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd}; 
		#stdcmd="python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/syst_samples.json -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --silent -o syst_plotter.root;"
	    #	bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd}; 
	;;
	PLOT_P2 )
	mv ${outdir}/analysis/plots/syst_plotter.root ${outdir}/analysis/plots/syst_plotter_orig.root
	python test/TopWidthAnalysis/createtWShapeUncs.py ~/work/TopWidth_era2015/analysis/plots/plotter.root ~/work/TopWidth_era2015/analysis/plots/syst_plotter.root  ${outdir}/analysis/plots/plotter.root;
	hadd ${outdir}/analysis/plots/syst_plotter.root tW_syst_plotter.root ${outdir}/analysis/plots/syst_plotter_orig.root 
	mv tW_syst_plotter.root ${outdir}/analysis/plots/tW_syst_plotter.root;	       
	
	#combined plots
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb        EE1b,EE2b,MM1b,MM2b,EM1b,EM2b    ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb        EE1b,MM1b,EM1b                   ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb        EE2b,MM2b,EM2b                   ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w lowptEE1b,lowptMM1b,lowptEM1b    ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w lowptEE2b,lowptMM2b,lowptEM2b    ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w highptEE1b,highptMM1b,highptEM1b ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w highptEE2b,highptMM2b,highptEM2b ${outdir}/analysis/plots/plotter.root
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	mkdir -p ${wwwdir}/comb
        cp plots/*.{root,png,pdf} ${wwwdir}/comb 
        cp test/index.php ${wwwdir}/comb
	;;
    HYPOTEST ) 
	mainHypo=1.0
	CATS=(
	    "lowptEE1b,lowptEE2b,highptEE1b,highptEE2b,lowptEM1b,lowptEM2b,highptEM1b,highptEM2b,lowptMM1b,lowptMM2b,highptMM1b,highptMM2b"
	    "lowptEE1b,lowptEE2b,highptEE1b,highptEE2b,lowptMM1b,lowptMM2b,highptMM1b,highptMM2b"
	    "lowptEM1b,lowptEM2b,highptEM1b,highptEM2b"
	)
	TAGS=("inc" "ll" "em")
	altHypo=(0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0)
	data=(-1.0 1.0 4.0)	
	#still to be debugged
        #cmd="${cmd} --addBinByBin 0.3" 
	for h in ${altHypo[@]}; do
	    for d in ${data[@]}; do
		for i in ${!TAGS[*]}; do

		    icat=${CATS[${i}]}
		    itag=${TAGS[${i}]}
		    
		    cmd="python test/TopWidthAnalysis/runHypoTestDatacards.py"
		    cmd="${cmd} --combine ${COMBINERELEASE}"
		    cmd="${cmd} --mainHypo=${mainHypo} --altHypo ${h} --pseudoData=${d}"
		    cmd="${cmd} -s tbart,tW --replaceDYshape"
		    cmd="${cmd} --dist incmlb"		    
		    cmd="${cmd} -i ${outdir}/analysis/plots/plotter.root"
		    cmd="${cmd} --systInput ${outdir}/analysis/plots/syst_plotter.root"
		    cmd="${cmd} -c ${icat}"
		    cmd="${cmd} --rebin 2"            
		    if [ "$h" == "2.2" ]; then
			if [ "$d" == "-1.0" ]; then
			    echo "    validation will be included"
			    cmd="${cmd} --doValidation"
			fi
		    fi

		    echo "Submitting ($mainHypo,$h,$d,$itag,$icat)"		
		    stdcmd="${cmd} -o ${outdir}/datacards_${itag}/"
		    bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd};
		    if [ "$itag" == "inc" ]; then
			if [ "$d" == "1.0" ]; then
			    echo "    injecting pseudo-data from nloproddec"
			    nlocmd="${cmd} --pseudoDataFromWgt nloproddec -o ${outdir}/datacards_${itag}_nloproddec"
			    bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${nlocmd};   
			    echo "    injecting pseudo-data from widthx4"
			    width4cmd="${cmd} --pseudoDataFromSim=t#bar{t}_widthx4 -o ${outdir}/datacards_${itag}_widthx4"
			    bsub -q ${queue} sh ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${width4cmd};
			fi
		    fi
		done
            done
	done
	;;
    PLOTHYPOTEST )
	summaryScript=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/summarizeHypoTestResults.py;
	cd ${COMBINERELEASE}/
	eval `scramv1 r -sh`
	cd -
	TAGS=("inc" "ll" "em" "inc_nloproddec" "inc_widthx4")
	for i in ${!TAGS[*]}; do
	    itag=${TAGS[${i}]}
	    python ${summaryScript} -i ${outdir}/datacards_${itag}/ --doCLs --recreateCLsSummary --doNuisances --doFitSummary	
	done	
	;;
    WWWHYPOTEST )
	TAGS=("_ll" "_em" "_inc" "_inc_nloproddec" "_inc_widthx4")
        for i in ${!TAGS[@]}; do
	    itag=${TAGS[${i}]}
            mkdir -p ${wwwdir}/hypo${itag}
            cp ${outdir}/datacards${itag}/*.{png,pdf} ${wwwdir}/hypo${itag};
            cp ${outdir}/datacards${itag}/hypotest_1.0vs2.2_data/*.{png,pdf} ${wwwdir}/hypo${itag};
            cp test/index.php ${wwwdir}/hypo${itag}
	done
	;;

esac
