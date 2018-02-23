#!/bin/bash

WHAT=$1; 
if [ "$#" -lt 1 ]; then 
    echo "steerTOP17010.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW/HYPOTEST> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        DY           - estimate DY scale factor from data";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo "        HYPOTEST     - create the datacards, steering scripts for hypothesis testing and submit to batch";
    echo "        PLOTHYPOTEST - create summaries of the hypothesis tests";
    echo "        WWWHYPOTEST  - move summaries to afs-web-based area"
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi


queue=workday
lumi=35922
lumiSpecs=""
lumiUnc=0.025
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/b312177
dataeosdir=/store/cmst3/group/top/ReReco2016/be52dbe_03Feb2017
summaryeosdir=/store/cmst3/group/top/TOP-17-010-final-v2/
AnalysisDir=/eos/cms/${summaryeosdir}/analysis
ChunksDir=${AnalysisDir}/Chunks
COMBINERELEASE=${HOME}/scratch0/CMSSW_7_4_7/src/
outdir=/afs/cern.ch/work/${myletter}/${whoami}/TOP-17-010-final-v2/
wwwdir=${HOME}/www/TOP-17-010/final-v2


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    TEST )
        file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_SingleTbar_tW/MergedMiniEvents_0_ext0.root
        #file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root
	python scripts/runLocalAnalysis.py -i ${file} \
            -q local -o /tmp/`whoami`/MC13TeV_SingleT_tW_test.root --era era2016 -m TOP-17-010::RunTop17010 --ch 0 --runSysts;
        python scripts/runTopWidthAnalysis.py -i /tmp/`whoami`/MC13TeV_SingleT_tW_test.root -o /tmp/`whoami`/Chunks -q local;
        ;;
    SEL )
        commonOpts="-q ${queue} -o ${summaryeosdir} --era era2016 -m TOP-17-010::RunTop17010 --ch 0 --runSysts --skipexisting";
	python scripts/runLocalAnalysis.py -i ${eosdir} ${commonOpts}     --only MC --farmappendix TOP17010MC;
	python scripts/runLocalAnalysis.py -i ${dataeosdir} ${commonOpts} --only Data --farmappendix TOP17010Data;
	;;
    MERGESEL )
	mkdir -p ${outdir}
	./scripts/mergeOutputs.py /eos/cms${summaryeosdir} True ${outdir};	
	;;
    PLOTSEL )
        commonOpts="-i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc} --doDataOverMC"
	python scripts/plotter.py ${commonOpts} 
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    TESTANA )
	python scripts/runTopWidthAnalysis.py -i root://eoscms//eos/cms/${summaryeosdir}/Chunks/MC13TeV_TTJets_12.root -o ${outdir}/analysis/Chunks -q local;
        ;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${ChunksDir} -q ${queue} --only MC;
        python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${ChunksDir} -q ${queue} --only Data13TeV_Single,Data13TeV_Double --farm TOP17010DataANA;
        python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${ChunksDir} -q ${queue} --only MuonEG --farm TOP17010DataMuEGANA;
	;;
    CHECKANA )
        for FARM in TOP17010ANA TOP17010DataANA TOP17010DataMuEGANA; do
            python scripts/checkAnalysisIntegrity.py ${CMSSW_BASE}/${FARM} /eos/cms/${summaryeosdir}/analysis/Chunks;
        done
        ;;
    MERGE )
	./scripts/mergeOutputs.py /eos/cms/${summaryeosdir}/analysis;
	;;
    BKG )
        opts="-j data/era2016/samples.json  -l ${lumi} ${lumiSpecs} --onlyData --doDataOverMC --mcUnc ${lumiUnc}"
	python scripts/plotter.py -i ${AnalysisDir} ${opts} --only mll -o dy_plotter.root; 
	python scripts/runDYRinRout.py --in ${AnalysisDir}/plots/dy_plotter.root --categs 1b,2b --out ${AnalysisDir}/plots/ > ${AnalysisDir}/plots/dysf.dat;
        opts="${opts} --procSF DY:${AnalysisDir}/plots/.dyscalefactors.pck --doDataOverMC --mcUnc ${lumiUnc}"
	python scripts/plotter.py -i ${AnalysisDir} ${opts} --only mll,evcount,dphilb,drlb,met,njets,ptlb,incmlb_w100 -o dysf_plotter.root; 
	python scripts/plotter.py -i ${AnalysisDir} ${opts} --only count --saveTeX -o count_plotter.root;
	;;
    PLOT )
        opts="-l ${lumi} ${lumiSpecs} --procSF DY:${AnalysisDir}/plots/.dyscalefactors.pck --mcUnc ${lumiUnc} --silent"
        python scripts/plotter.py -i ${AnalysisDir}  -j data/era2016/samples.json       ${opts};
        #python scripts/plotter.py -i ${AnalysisDir}  -j data/era2016/syst_samples.json  ${opts} -o syst_plotter.root;
        ;;
    COMBPLOT )
	#combined plots
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb:ptlbinc            EE1b,EE2b,MM1b,MM2b,EM1b,EM2b     ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb:ptlb1b             EE1b,MM1b,EM1b                    ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb:ptlb2b             EE2b,MM2b,EM2b                    ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py drlb:drlb1b             EE1b,MM1b,EM1b                    ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py drlb:drlb2b             EE2b,MM2b,EM2b                    ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py drlb:drlbinc            EE1b,EE2b,MM1b,MM2b,EM1b,EM2b     ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlb1blowpt  EE1blowpt,MM1blowpt,EM1blowpt     ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlb2blowpt  EE2blowpt,MM2blowpt,EM2blowpt     ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlb1bhighpt EE1bhighpt,MM1bhighpt,EM1bhighpt  ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlb2bhighpt EE2bhighpt,MM2bhighpt,EM2bhighpt  ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlblowpt    EE1blowpt,MM1blowpt,EM1blowpt,EE2blowpt,MM2blowpt,EM2blowpt  ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlbhighpt   EE1bhighpt,MM1bhighpt,EM1bhighpt,EE2bhighpt,MM2bhighpt,EM2bhighpt  ${AnalysisDir}/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100:mlbinc      EE1blowpt,MM1blowpt,EM1blowpt,EE2blowpt,MM2blowpt,EM2blowpt,EE1bhighpt,MM1bhighpt,EM1bhighpt,EE2bhighpt,MM2bhighpt,EM2bhighpt  ${AnalysisDir}/plots/plotter.root
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${AnalysisDir}/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	mkdir -p ${wwwdir}/comb
        cp plots/*.{root,png,pdf} ${wwwdir}/comb 
        cp test/index.php ${wwwdir}/comb
	;;
    HYPOTEST ) 
        mainHypo=(20 100 400)        
	altHypo=(20 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 350 400)        
	data=(-1 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 350 400)
        queue=2nd       
        for mh in ${mainHypo[@]}; do
	    for h in ${altHypo[@]}; do
	        for d in ${data[@]}; do		    		       
		    cmd="python test/TopWidthAnalysis/runHypoTestDatacards.py"
		    cmd="${cmd} --combine ${COMBINERELEASE}"
		    cmd="${cmd} --mainHypo ${mh} --altHypo ${h} --pseudoData ${d}"		    
                    cmd="${cmd} -s tbart"
		    cmd="${cmd} --dist incmlb"		    
		    cmd="${cmd} --nToys 2000"
		    cmd="${cmd} -i /eos/cms/${summaryeosdir}/plotter/plotter.root"
		    cmd="${cmd} --systInput /eos/cms/${summaryeosdir}/plotter/syst_plotter.root"
                    cats="EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,EE1blowpt,EE2blowpt,EE1bhighpt,EE2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt"
		    cmd="${cmd} -c ${cats}"
		    cmd="${cmd} --rebin 2"            
                    cmd="${cmd} --doValidation"
		    stdcmd="${cmd} -o ${outdir}/datacards_${mh}/"
		    echo "Submitting ($mh,$h,$d)"		
                    echo ${stdcmd}
		    bsub -q ${queue} sh ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd};
                done
            done
	done
	;;
    PLOTHYPOTEST )
	summaryScript=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/summarizeHypoTestResults.py;
	cd ${COMBINERELEASE}/
	eval `scramv1 r -sh`
	cd -
	TAGS=("inc") # "ll" "em" "inc_nloproddec" "inc_widthx4")
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
            cp ${outdir}/datacards${itag}/hypotest_100vs220_data/*.{png,pdf} ${wwwdir}/hypo${itag};
            cp test/index.php ${wwwdir}/hypo${itag}
	done
	;;
    IMPACTS )
        echo "Starting nuisance impact estimation..."
        cd /afs/cern.ch/work/e/ecoleman/CMSSW_7_4_7/
        eval `scramv1 r -sh`
        cd ${outdir}/$2/

        combineTool.py -M Impacts -d workspace.root -m 172.5 --doInitialFit --robustFit 1
        combineTool.py -M Impacts -d workspace.root -m 172.5 --robustFit 1 --doFits --job-mode lxbatch --task-name nuis-impacts --sub-opts='-q 1nw'
        ;;
    FINALIMP )
        echo "Finalizing nuisance impacts..."
        cd /afs/cern.ch/work/e/ecoleman/CMSSW_7_4_7/
        eval `scramv1 r -sh`
        cd ${outdir}/$2/
       
        combineTool.py -M Impacts -d workspace.root -m 172.5 -o impacts.json
        plotImpacts.py -i impacts.json -o impacts
        ;;
esac
