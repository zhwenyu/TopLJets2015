#!/bin/bash

WHAT=$1;
if [ "$#" -ne 1 ]; then
    echo "steerTOPJS.sh <SEL/PLOTSEL/WWWSEL>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    exit 1;
fi

export LSB_JOB_REPORT_MAIL=N


queue=8nh
githash=b312177
lumi=36460
lumiSpecs="" #--lumiSpecs EE:11391"
lumiUnc=0.027
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
summaryeosdir=/eos/user/${myletter}/${whoami}/analysis/TopJetShapes/${githash}
outdir=${summaryeosdir}
wwwdir=/eos/user/${myletter}/${whoami}/www/cms/TopJS/


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        scram b -j 8 && python scripts/runLocalAnalysis.py -i ${eosdir}/MC13TeV_TTJets/MergedMiniEvents_0.root --tag MC13TeV_TTJets -o analysis.root --era era2016 -m TOPJetShape::RunTopJetShape --debug;
        ;;

    NORMCACHE )
        python scripts/produceNormalizationCache.py -i ${eosdir} -o data/era2016/genweights.root;
        ;;

    FULLSEL ) #TODO -o ${summaryeosdir} for all
        cd batch;
        python ../scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --skip Data13TeV_Double,Data13TeV_MuonEG,MC13TeV_TTJets2l2nu,MC13TeV_TTJets_cflip,MC13TeV_TTJets_m,MC13TeV_TTJets_width;
        #python ../scripts/runLocalAnalysis.py -i ${eosdir}_syst -q ${queue} -o ${summaryeosdir}_syst --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --only MC13TeV_TTJets_isr,MC13TeV_TTJets_fsr,MC13TeV_TTJets_hdamp,MC13TeV_TTJets_herwig,MC13TeV_TTJets_ue;
        python ../scripts/runLocalAnalysis.py -i ${eosdir}_qcd -q ${queue} -o ${summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting;
        cd -;
        ;;

    FULLSELSYST )
        cd batch;
        python ../scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${summaryeosdir}_expsyst --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --only MC13TeV_TTJets --systVar all --exactonly;
        cd -;
        ;;

    MERGE )
        python scripts/mergeOutputs.py ${summaryeosdir} True;
        #TODO Not needed if everything goes in one directory
        #python scripts/mergeOutputs.py ${summaryeosdir}_syst True;
        #cp ${summaryeosdir}_syst/*.root ${summaryeosdir}
        #python scripts/mergeOutputs.py ${summaryeosdir}_qcd True;
        #cp ${summaryeosdir}_qcd/*.root ${summaryeosdir}
        ;;

    PLOTSEL )
        rm -r plots
        commonOpts="-i ${summaryeosdir} -j data/era2016/samples.json,data/era2016/qcd_samples.json --systJson data/era2016/syst_samples.json -l ${lumi} --mcUnc ${lumiUnc} --rebin 1"
        python scripts/plotter.py ${commonOpts} --outDir plots;
        ;;

    TESTPLOTSEL )
        commonOpts="-i ${summaryeosdir} -j data/era2016/samples.json,data/era2016/qcd_samples.json --systJson data/era2016/syst_samples.json -l ${lumi}"
        python scripts/plotter.py ${commonOpts} --outDir plots/test --only L4_1l4j2b2w_njets,js_tau32_puppi;
        ;;

    WWWSEL )
        rm -r ${wwwdir}/sel
        mkdir -p ${wwwdir}/sel
        cp plots/*.{png,pdf} ${wwwdir}/sel
        cp test/index.php ${wwwdir}/sel
        ;;

    BINNING )
        python test/TopJSAnalysis/optimizeUnfoldingMatrix.py -i eos --obs all
        ;;

    WWWBINNING )
        rm -r ${wwwdir}/binning
        mkdir -p ${wwwdir}/binning
        cp unfolding/optimize/*.{png,pdf} ${wwwdir}/binning
        cp test/index.php ${wwwdir}/binning
        ;;

    FILL )
        cd batch;
        python ../test/TopJSAnalysis/fillUnfoldingMatrix.py -i /eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/;
        cd -;
        ;;
    
    FILLWEIGHTS )
        cd batch;
        python ../test/TopJSAnalysis/fillUnfoldingMatrix.py -i /eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/ --only MC13TeV_TTJets --nweights 20 -q 8nh;
        cd -
        ;;
    
    MERGEFILL )
        ./scripts/mergeOutputs.py unfolding/fill True
        ;;
        
    TOYUNFOLDING )
        for OBS in mult width ptd ptds ecc tau21 tau32 tau43 zg zgxdr zgdr ga_width ga_lha ga_thrust c1_02 c1_05 c1_10 c1_20 c2_02 c2_05 c2_10 c2_20 c3_02 c3_05 c3_10 c3_20
        do
          mkdir -p unfolding/toys/farm/${OBS}
          cp test/TopJSAnalysis/testUnfold0Toys.C unfolding/toys/farm/${OBS}
          root -l -b -q "unfolding/toys/farm/${OBS}/testUnfold0Toys.C++(\"${OBS}\", 1000)"&
        done
        ;;
        
        

esac
