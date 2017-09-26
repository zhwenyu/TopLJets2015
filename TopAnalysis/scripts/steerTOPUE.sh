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


queue=workday
lumi=35922
lumiSpecs="" #--lumiSpecs EE:11391"
lumiUnc=0.025
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/b312177
dataeos=/store/cmst3/group/top/ReReco2016/be52dbe_03Feb2017
summaryeosdir=/store/cmst3/group/top/TopUE
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/UEanalysis/
wwwdir=~/www/TOP-17-015


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root
        #file==root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets2l2nu_amcatnlo/MergedMiniEvents_1_ext0.root
	#file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/Data13TeV_MuonEG_2016D/MergedMiniEvents_0.root
	#analysisWrapper \
	#    --in ${file} \
	#    --out ue_test.root \
	#    --era ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016 \
	#    --method TOP-UE::RunTopUE \
	#    --runSysts \
	#    --ch 0;
        for step in 1 2; do
	    python test/TopUEAnalysis/runUEanalysis.py -i ue_test.root --step ${step} --ptThr 1.0,0.9  --obs chmult --slice ptll=0,9999. --reg ptll=awa -o ./UEanalysis_test;
        done
	#python test/TopUEAnalysis/runUEanalysis.py --step 1 -o ./UEanalysis_test;
	#python test/TopUEAnalysis/runUEanalysis.py -i ue_test.root      --step 2 -q local -o ./UEanalysis_test;
	#python test/TopUEAnalysis/showFastFinalDistributions.py UEanalysis_test/analysis/Chunks/ue_test.root --cfg ./UEanalysis_test/analysisaxiscfg.pck
	;;

    SEL )
        commonOpts="-q ${queue} -o ${summaryeosdir}      --era era2016 -m TOP-UE::RunTopUE --ch 0 --runSysts";
	python scripts/runLocalAnalysis.py -i ${eosdir}  --farmappendix TopUEMC ${commonOpts} --only MC;
	python scripts/runLocalAnalysis.py -i ${dataeos} --farmappendix TopUEMC ${commonOpts} --only Data;
	;;

    MERGESEL )
	mkdir -p ${outdir}
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount eos;
	./scripts/mergeOutputs.py eos/cms${summaryeosdir} True ${outdir};	
	/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse umount eos;
	;;
    PLOTSEL )
	commonOpts="-i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts} --only mll --outName mll_plotter.root;	
     	python scripts/runDYRinRout.py --in ${outdir}/plots/mll_plotter.root --categs "0t,1t,"  --out ${outdir}/plots/ > ${outdir}/plots/dy.dat;
	python scripts/plotter.py ${commonOpts} --procSF DY:${outdir}/plots/.dyscalefactors.pck --only njets --rebin 7 --saveTeX --outName count_plotter.root;
	python scripts/plotter.py ${commonOpts} --procSF DY:${outdir}/plots/.dyscalefactors.pck;
	python scripts/plotter.py ${commonOpts} --only nbtags,rho,nvtx,0t,1t
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	eosprefix=root://eoscms//eos/cms

	base="${eosprefix}/${summaryeosdir}/Chunks/MC13TeV_TTJets"
	baseFiles=${base}_0.root,${base}_1.root,${base}_2.root,${base}_3.root,${base}_4.root

	echo "Preparing analysis configuration based on ${baseFiles}"
        obs=("chmult" "chavgpt" "chavgpz" "chfluxz" "chflux" "sphericity" "aplanarity" "C" "D")
        analyses=(
            "" 
            "--reg ptttbar=awa" 
            "--reg ptttbar=tow" 
            "--reg ptttbar=tra" 
            "--reg ptll=awa" 
            "--reg ptll=tow" 
            "--reg ptll=tra"
            "--slice nj=0,1 " 
            "--slice nj=0,1 --reg ptttbar=awa" 
            "--slice nj=0,1 --reg ptttbar=tow" 
            "--slice nj=0,1 --reg ptttbar=tra" 
            "--slice nj=0,1 --reg ptll=awa" 
            "--slice nj=0,1 --reg ptll=tow" 
            "--slice nj=0,1 --reg ptll=tra"
            "--slice nj=1,2 " 
            "--slice nj=1,2 --reg ptttbar=awa" 
            "--slice nj=1,2 --reg ptttbar=tow" 
            "--slice nj=1,2 --reg ptttbar=tra" 
            "--slice nj=1,2 --reg ptll=awa" 
            "--slice nj=1,2 --reg ptll=tow" 
            "--slice nj=1,2 --reg ptll=tra"
            "--slice nj=2,999" 
            "--slice nj=2,999 --reg ptttbar=awa" 
            "--slice nj=2,999 --reg ptttbar=tow" 
            "--slice nj=2,999 --reg ptttbar=tra" 
            "--slice nj=2,999 --reg ptll=awa" 
            "--slice nj=2,999 --reg ptll=tow" 
            "--slice nj=2,999 --reg ptll=tra"
            "--slice ptttbar=0,20 " 
            "--slice ptttbar=0,20 --reg ptttbar=awa" 
            "--slice ptttbar=0,20 --reg ptttbar=tow" 
            "--slice ptttbar=0,20 --reg ptttbar=tra" 
            "--slice ptttbar=20,60" 
            "--slice ptttbar=20,60--reg ptttbar=awa" 
            "--slice ptttbar=20,60--reg ptttbar=tow" 
            "--slice ptttbar=20,60--reg ptttbar=tra" 
            "--slice ptttbar=60,120" 
            "--slice ptttbar=60,120 --reg ptttbar=awa" 
            "--slice ptttbar=60,120 --reg ptttbar=tow" 
            "--slice ptttbar=60,120 --reg ptttbar=tra" 
            "--slice ptttbar=120,9999" 
            "--slice ptttbar=120,9999 --reg ptttbar=awa" 
            "--slice ptttbar=120,9999 --reg ptttbar=tow" 
            "--slice ptttbar=120,9999 --reg ptttbar=tra" 
            "--slice ptll=0,60"
            "--slice ptll=0,60 --reg ptll=awa" 
            "--slice ptll=0,60 --reg ptll=tow" 
            "--slice ptll=0,60 --reg ptll=tra"
            "--slice ptll=60,120"
            "--slice ptll=60,120 --reg ptll=awa" 
            "--slice ptll=60,120 --reg ptll=tow" 
            "--slice ptll=60,120 --reg ptll=tra"
            "--slice ptll=120,9999"
            "--slice ptll=120,9999 --reg ptll=awa" 
            "--slice ptll=120,9999 --reg ptll=tow" 
            "--slice ptll=120,9999 --reg ptll=tra"
        )
        for o in "${obs[@]}"; do
            for a in "${analyses[@]}"; do
                options="--ptThr 0.9,0.9 --obs ${o} ${a}"
                echo $options
	        python test/TopUEAnalysis/runUEanalysis.py -i ${baseFiles}     --step 1             ${options} -o ./UEanalysis;
            done
        done

        queue=workday
	echo "Filling the histograms for unfolding"
        a=(`ls UEanalysis/*/*/analysiscfg.pck`)
        for i in ${a[@]}; do
            echo ${i};
            continue
            python test/TopUEAnalysis/runUEanalysis.py -i ${summaryeosdir}/Chunks --step 2 -q ${queue} -o `dirname ${i}`;
        done
        
	;;
    MERGEANA )
        a=(`ls UEanalysis/*/*/analysiscfg.pck`)
        for i in ${a[@]}; do
            dir=`dirname ${i}`
            if [ -f "${dir}/condor.sub" ]; then
                echo "Checking results for ${dir}"
                python scripts/checkAnalysisIntegrity.py ${dir}/FARM-UEANA/ ${dir}/Chunks/
	        ./scripts/mergeOutputs.py ${dir} True 
	        commonOpts="-l ${lumi} --mcUnc ${lumiUnc} --procSF DY:${outdir}/plots/.dyscalefactors.pck";
	        python scripts/plotter.py -i ${dir} -j data/era2016/samples.json      ${commonOpts} --only _0;
                python scripts/plotter.py -i ${dir} -j data/era2016/samples.json      ${commonOpts} --silent;
	        python scripts/plotter.py -i ${dir} -j data/era2016/syst_samples.json ${commonOpts} --silent --outName syst_plotter.root;	
            fi
        done
	;;
    PLOTANA )

        mkdir -p UEanalysis/analysis/plots/rawana
	commonOpts="-l ${lumi} --saveLog --mcUnc ${lumiUnc} --procSF DY:${outdir}/plots/.dyscalefactors.pck";
	python scripts/plotter.py -i UEanalysis/analysis -j data/era2016/samples.json      ${commonOpts};
	python scripts/plotter.py -i UEanalysis/analysis -j data/era2016/syst_samples.json ${commonOpts} --silent --outName syst_plotter.root;	
        python test/TopUEAnalysis/UETools.py -o UEanalysis/analysis/plots/ -i  UEanalysis/analysis/MC13TeV_TTJets.root
	python test/TopUEAnalysis/showFinalRecoDistribution.py \
            UEanalysis/analysis/plots/plotter.root \
            UEanalysis/analysis/plots/syst_plotter.root \
            --out UEanalysis/analysis/plots/reco;
        #python test/TopUEAnalysis/showFinalUnfoldedDistribution.py \
        #    UEanalysis/analysis/plots/plotter.root \
        #    UEanalysis/analysis/plots/syst_plotter.root \
        #    --out UEanalysis/analysis/plots/rawana \
        #    --reco;
	
	;;
    WWWANA )
	mkdir -p ${wwwdir}/rawana        
        #cp UEanalysis/analysis/plots/*.{png,pdf,dat} ${wwwdir}/rawana
	#cp  UEanalysis/analysis/plots/rawana/*.{png,pdf,dat} ${wwwdir}/rawana
        #cp test/index.php ${wwwdir}/rawana
        mkdir -p ${wwwdir}/reco
        cp UEanalysis/analysis/plots/reco/*{png,pdf} ${wwwdir}/reco
        cp test/index.php ${wwwdir}/reco
	;;
    UNFOLD )
        commonOpts="--plotter UEanalysis/analysis/plots/plotter.root --syst UEanalysis/analysis/plots/syst_plotter.root -d UEanalysis/analysis/Chunks/"
        #python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 0

        for ovar in chmult chflux chfluxz chavgpt chavgpz sphericity aplanarity C D; do                  

            #inclusive
            #python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 1 --histo ${ovar}_None_0_inc;            
            #python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 2 --histo ${ovar}_None_0_inc;
            #python test/TopUEAnalysis/showUnfoldSummary.py -i UEanalysis/unfold/${ovar}_None_0_inc.root;                                   
        
            #differential
            for svar in nj chmult; do

                if [ "$svar" = "$ovar" ]
                then
                    continue
                fi

                bins=`seq 1 3`;
                if [ "$svar" = "chmult" ] 
                then
                    bins=`seq 1 9`;
                fi

                for bin in $bins; do
                    echo $ovar $svar $bin
                    python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 1 --histo ${ovar}_${svar}_${bin}_inc;            
                    python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 2 --histo ${ovar}_${svar}_${bin}_inc;
                    python test/TopUEAnalysis/showUnfoldSummary.py -i UEanalysis/unfold/${ovar}_${svar}_${bin}_inc.root;                                   
                done
            done
        done
        ;;
    PLOTUNFOLD )
        python test/TopUEAnalysis/showFinalUnfoldedDistribution.py \
            UEanalysis/unfold UEanalysis/analysis/plots/plotter.root \
            UEanalysis/analysis/plots/syst_plotter.root \
            --out UEanalysis/unfold/;
        ;;
    WWWUNFOLD )
        mkdir -p ${wwwdir}/unfold;
        cp UEanalysis/unfold/*.{png,pdf,dat} ${wwwdir}/unfold;
        cp test/index.php ${wwwdir}/unfold;
        ;;

esac
