#!/bin/bash

WHAT=$1; 
TAGANA=$2
if [ "$#" -lt 1 ]; then 
    echo "steerTOPUE.sh <OPTIONL>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge selection outputs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        PREPAREANA   - prepares configuration files to run the analysis";
    echo "        SUBMITANA    - submits a specific analysis (tag must be provided as extra argument e.g. chmult/inc)";
    echo "        MERGEANA     - checks and merges analysis outputs (tag must be provided as extra argument e.g. chmult/inc)";
    echo "        MVANA        - use to move to another area mounted under store (tag must be provided as extra argument e.g. chmult/inc)";
    echo "        UNFOLDANA    - unfold results in a given directory (directory must be provided as extra argument e.g. store/TOP-17-015/chmult/inc)";
    echo "        WWWANA       - move summary plots to web area (directory must be provided as extra argument e.g. store/TOP-17-015/chmult/inc)";
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
markuseos=/eos/user/m/mseidel/ReReco2016/b312177_merged/
efeeos=/store/group/phys_top/efe/output_ue_root_alpha_s_isr
efeeos2=/store/group/phys_top/efe/
summaryeosdir=/store/cmst3/group/top/TOP-17-015
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/UEanalysis/
wwwdir=~/www/TOP-17-015


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
	file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root
        file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets2l2nu_noSC/MergedMiniEvents_0_ext0.root
        #file==root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets2l2nu_amcatnlo/MergedMiniEvents_1_ext0.root
	#file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/Data13TeV_MuonEG_2016D/MergedMiniEvents_0.root
	outFile=ue_test_nosc.root
        #analysisWrapper \
	#    --in ${file} \
	#    --out ue_test.root \
	#    --era ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016 \
	#    --method TOP-UE::RunTopUE \
	#    --runSysts \
	#    --ch 0;
        for step in 1; do
            continue
	    #python test/TopUEAnalysis/runUEanalysis.py -i ue_test.root --step ${step} --ptThr 1.0,0.9  --obs chmult --slice ptll=0,9999. --reg ptll=awa -o ./UEanalysis_test;
            python test/TopUEAnalysis/runUEanalysis.py -i ue_test.root  --ptThr 0.9,0.9 --step ${step} --obs chmult -o ./UEanalysis_test;
        done
	#python test/TopUEAnalysis/runUEanalysis.py --step 1 -o ./UEanalysis_test;
	#python test/TopUEAnalysis/runUEanalysis.py -i ue_test.root      --step 2 -q local -o ./UEanalysis_test;
	#python test/TopUEAnalysis/showFastFinalDistributions.py UEanalysis_test/analysis/Chunks/ue_test.root --cfg ./UEanalysis_test/analysisaxiscfg.pck
	;;

    SEL )
        commonOpts="-q ${queue} -o ${summaryeosdir}      --era era2016 -m TOP-UE::RunTopUE --ch 0 --runSysts";
	python scripts/runLocalAnalysis.py -i ${eosdir}  --farmappendix TopUEMC ${commonOpts} --only MC13TeV;
	python scripts/runLocalAnalysis.py -i ${dataeos} --farmappendix TopUEMC ${commonOpts} --only Data;        
        python scripts/runLocalAnalysis.py -i ${markuseos} --farmappendix TopUEMC ${commonOpts} --only asfsr,herwig7,sherpa;
        python scripts/runLocalAnalysis.py -i ${efeeos} --farmappendix TopUEMCASISR  ${commonOpts} --ignore "~";
        python scripts/runLocalAnalysis.py -i ${efeeos2} --farmappendix TopUEMCMPICR ${commonOpts} --only mpi_off,cr_off --ignore "~";        
        #python scripts/runLocalAnalysis.py -i ${efeeos3} --farmappendix TopUEMCMPICR ${commonOpts} --ignore "~";        
	;;
    CHECKSELINTEG )
        for farm in TopUEMC TopUEMCASISR TopUEMCMPICR; do
            python scripts/checkAnalysisIntegrity.py ${CMSSW_BASE}/FARMTOP-17-015${farm} /eos/cms/${summaryeosdir}/Chunks
        done
        ;;

    MERGESEL )
	mkdir -p ${outdir}
	./scripts/mergeOutputs.py /eos/cms${summaryeosdir} True ${outdir};	
	;;

    PLOTSEL )
	commonOpts="-i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts} --only mll --outName mll_plotter.root;	
     	python scripts/runDYRinRout.py --in ${outdir}/plots/mll_plotter.root --categs "0t,1t,"  --out ${outdir}/plots/ > ${outdir}/plots/dy.dat;
	python scripts/plotter.py ${commonOpts} --procSF DY:${outdir}/plots/.dyscalefactors.pck --only njets --rebin 7 --saveTeX --outName count_plotter.root;
	python scripts/plotter.py ${commonOpts} --procSF DY:${outdir}/plots/.dyscalefactors.pck;
	python scripts/plotter.py ${commonOpts} --only nbtags,rho,nvtx,0t,1t
	;;
    PLOTSELPAPER )
        python scripts/plotter.py -i ${outdir} -outName paper_plotter.root \
            --puNormSF puwgtctr  -j data/era2016/samples.json \
            -l ${lumi} --procSF DY:${outdir}/plots/.dyscalefactors.pck --noRatio --noUncs\
            --cmsLabel "#bf{CMS}" --only mll_EM,ptll_EM,njets_EM --formats "pdf,png,root";
        cp -v ${outdir}/plots/{mll,ptll,njets}_EM.* ${wwwdir}/sel/
        ;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;

    PREPAREANA )
	eosprefix=root://eoscms//eos/cms

	base="${eosprefix}/${summaryeosdir}/Chunks/MC13TeV_TTJets"
	baseFiles=${base}_0.root,${base}_1.root,${base}_2.root,${base}_3.root,${base}_4.root

	echo "Preparing analysis configuration based on ${baseFiles} - this will take a long time..."
        obs=("C" "D" "sphericity" "aplanarity" "chmult" "chavgpt" "chavgpz" "chfluxz" "chflux" "maxRap" "rapDist")
        analyses=(
            "" 
            "--slice nj=0,1" 
            "--slice nj=1,2" 
            "--slice nj=2,999" 
            "--slice nj=0,1   --reg ptll=tra" 
            "--slice nj=0,1   --reg ptll=tow" 
            "--slice nj=0,1   --reg ptll=awa" 
            "--slice nj=1,2   --reg ptll=tra" 
            "--slice nj=1,2   --reg ptll=tow" 
            "--slice nj=1,2   --reg ptll=awa" 
            "--slice nj=2,999 --reg ptll=tra" 
            "--slice nj=2,999 --reg ptll=tow" 
            "--slice nj=2,999 --reg ptll=awa" 
            "--slice mll=0,60" 
            "--slice mll=60,120" 
            "--slice mll=120,200" 
            "--slice mll=200,9999" 
            "--reg ptll=awa" 
            "--reg ptll=tow" 
            "--reg ptll=tra"
            "--slice ptll=0,20"
            "--slice ptll=0,20 --reg ptll=awa" 
            "--slice ptll=0,20 --reg ptll=tow" 
            "--slice ptll=0,20 --reg ptll=tra"
            "--slice ptll=20,40"
            "--slice ptll=20,40 --reg ptll=awa" 
            "--slice ptll=20,40 --reg ptll=tow" 
            "--slice ptll=20,40 --reg ptll=tra"
            "--slice ptll=40,80"
            "--slice ptll=40,80 --reg ptll=awa" 
            "--slice ptll=40,80 --reg ptll=tow" 
            "--slice ptll=40,80 --reg ptll=tra"
            "--slice ptll=80,120"
            "--slice ptll=80,120 --reg ptll=awa" 
            "--slice ptll=80,120 --reg ptll=tow" 
            "--slice ptll=80,120 --reg ptll=tra"
            "--slice ptll=120,9999"
            "--slice ptll=120,9999 --reg ptll=awa" 
            "--slice ptll=120,9999 --reg ptll=tow" 
            "--slice ptll=120,9999 --reg ptll=tra"
        )
        obs=("chavgpz")
        analyses=("--slice ptll=80,120"
            "--slice ptll=40,80 --reg ptll=tow"
            "--slice ptll=40,80 --reg ptll=tra"
        )
        
        obs=("chavgpt")
        analyses=("--slice ptll=0,20"
            "--slice ptll=20,40"
            "--slice ptll=0,20 --reg ptll=tow"
            "--slice ptll=0,20 --reg ptll=tra"
            "--slice ptll=0,20 --reg ptll=awa"
            "--reg ptll=tow"
            "--reg ptll=tra"
        )

        for o in "${obs[@]}"; do
            for a in "${analyses[@]}"; do
                options="--ptThr 0.9,0.9 --obs ${o} ${a}"
                if [[ $a == *"--reg"* ]]; then
                    if [ "$o" == "sphericity" ] || [ "$o" == "aplanarity" ] || [ "$o" == "C" ] || [ "$o" == "D" ]; then
                        echo "Skipping ${a} for ${o} as this is an inclusive observable";
                        continue
                    fi
                fi    
                if [[ $a == *"chmult"* ]]; then
                    if [ "$o" == "chmult" ];then
                        echo "Skipping ${a} for ${o} as this is the variable being sliced"
                    fi
                fi

                #configure histograms (binning etc.)
	        python test/TopUEAnalysis/runUEanalysis.py -i ${baseFiles} --step 1  ${options} -o ./UEanalysis;

                #create condor
                dir=`cat lastUE.dat`
                python test/TopUEAnalysis/runUEanalysis.py -i ${summaryeosdir}/Chunks --step 2 -q ${queue} -o ${dir} --dryRun;
                
                #submit
                cd ${dir}
                condor_submit condor.sub;
                cd -
            done
        done
        
	;;
    SUBMITANA )
        queue=longlunch
	echo "Filling the histograms for unfolding"
        cd UEanalysis/${TAGANA};
        condor_submit condor.sub;
        cd -
        ;;
    SUBMITSPECIALANA )
        obs=("sphericity" "aplanarity" "C" "D" "chmult" "chavgpt" "chavgpz" "chfluxz" "chflux")
        for i in ${obs[@]}; do
            a=(`ls store/TOP-17-015/${i}`)
            #a=("nj=0,1_ptll=awa" "nj=0,1_ptll=tow" "nj=0,1_ptll=tra" "nj=1,2_ptll=awa" "nj=1,2_ptll=tow" "nj=1,2_ptll=tra" "nj=2,999_ptll=awa" "nj=2,999_ptll=tow" "nj=2,999_ptll=tra")
            for j in ${a[@]}; do
                dir=store/TOP-17-015/${i}/${j};
                if [ -d ${dir} ]; then
                    echo "Preparing analysis cfg for $dir"
                    mkdir -p UEanalysis/${i}/${j};
                    cp -v ${dir}/analysis*.pck UEanalysis/${i}/${j};

                    echo "Creating jobs for special MC (gen only)"
                    python test/TopUEAnalysis/runUEanalysis.py -i /eos/cms/store/cmst3/group/top/TopUE_extra/Chunks --only asfsr --step 2 -q ${queue} -o UEanalysis/${i}/${j} --dryRun;
                    cd UEanalysis/${i}/${j};
                    condor_submit condor.sub;
                    cd -;
                    #sleep 30s;
                fi
            done
        done
        ;;

    CHECKANA )
        dir=UEanalysis/${TAGANA}
        python scripts/checkAnalysisIntegrity.py ${dir}/FARM-UEANA/ ${dir}/Chunks/ #> ${dir}/ana_integ.dat
        ;;

    MVANA)
        mkdir -p store/TOP-17-015/${TAGANA}
        rsync -axu --remove-source-files --delete-after --progress UEanalysis/${TAGANA}/* store/TOP-17-015/${TAGANA}/
        rm -rf UEanalysis/${TAGANA}/
        rm -rf store/TOP-17-015/${TAGANA}/FARM-UEANA/
        ;;

    MERGEANA )
        dir=${TAGANA}
        echo "Checking results for ${dir}"
	./scripts/mergeOutputs.py ${dir} True 
	commonOpts="-l ${lumi} --mcUnc ${lumiUnc} --procSF DY:${outdir}/plots/.dyscalefactors.pck";
	python scripts/plotter.py -i ${dir} -j data/era2016/samples.json      ${commonOpts} --only _0;
        python scripts/plotter.py -i ${dir} -j data/era2016/samples.json      ${commonOpts} --silent;
	python scripts/plotter.py -i ${dir} -j data/era2016/syst_samples.json ${commonOpts} --silent --outName syst_plotter.root;	            
	;;


    UNFOLDANA )
        dir=$TAGANA
        commonOpts="-o ${dir}/unfold --plotter ${dir}/plots/plotter.root --syst ${dir}/plots/syst_plotter.root -d ${dir}/Chunks/"            
        python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 0;
        python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 1;
        python test/TopUEAnalysis/runUEUnfolding.py ${commonOpts} -s 2;
        python test/TopUEAnalysis/showUnfoldSummary.py -i ${dir}/unfold/unfold_summary.root;
        python test/TopUEAnalysis/showFinalDistributions.py \
            --cfg ${dir}/analysiscfg.pck --cmsLabel "#bf{CMS}"\
            ${dir}/unfold/unfold_summary.root \
            ${dir}/plots/plotter.root \
            ${dir}/plots/syst_plotter.root;        
        ;;

    WWWANA )

        tks=(`echo $TAGANA | tr "/" "\n"`)
        ntks=${#tks[@]}
        tag="${tks[$ntks-2]}_${tks[$ntks-1]}"
        tag=${tag//\=/_eq_}
        tag=${tag//,/_}
        tag=${tag//./p}
        wwwdir="${wwwdir}/ana_${tks[$ntks-2]}"
        mkdir -p ${wwwdir}
        a=(`ls $TAGANA/unfold/*.{png,pdf,dat}`)
        b=(`ls $TAGANA/*.{png,pdf}`)
        a+=( "${a[@]}" "${b[@]}" )
        for i in ${a[@]}; do            
            cp ${i} ${wwwdir}/${tag}_`basename ${i}`;
        done
        cp test/index.php ${wwwdir}/

	;;

    PROFILEANA )
        dir=$TAGANA
        tks=(`echo $TAGANA | tr "/" "\n"`)
        ntks=${#tks[@]}
        ana="${tks[$ntks-1]}"
        for s in 1 2; do
            python test/TopUEAnalysis/showFinalProfiles.py -i ${TAGANA} -s ${s} --doPull --cmsLabel "#bf{CMS}";
        done
        if [[ $dir = *"chavgp"* ]]; then
            python test/TopUEAnalysis/showFinalProfiles.py -i ${TAGANA} -s 3 --doPull --cmsLabel "#bf{CMS}";
        fi
        a=(`ls ${dir}/*ueprofile*.{png,pdf}`)
        for i in ${a[@]}; do
            cp -v ${i} ${wwwdir}/ana_${ana}/${ana}_`basename ${i}`;
        done
        ;;
esac
