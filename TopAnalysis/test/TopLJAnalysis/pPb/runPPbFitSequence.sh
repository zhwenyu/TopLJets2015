#!/bin/bash

WHAT=$1;

case $WHAT in

    TESTSEL)

	#for sample in MC8.16TeV_TTbar_pPb_Pohweg MC8.16TeV_WJets_pPb MC8TeV_WJets_pp; do 
	for sample in MC8.16TeV_TTbar_pPb_Pohweg; do
	    for wjjOrder in drjj mjj sumpt; do
		for thadOrder in dr dm2tlep dm2pdg; do
	            python prepareWorkspace.py  -d ${sample} --wjjOrder ${wjjOrder} --thadOrder ${thadOrder};
	            mv plots/${sample}/controlplots.root plots/${sample}/controlplots_${wjjOrder}_${thadOrder}.root;
                    python ./checkMatchingEfficiency.py plots/${sample}/controlplots_${wjjOrder}_${thadOrder}.root #> plots/${sample}/algoeff_${wjjOrder}_${thadOrder}.tex 
		done
	    done

	    rm plots/${sample}/*.{png,pdf}

	    python checkMatchingDistributions.py "min#DeltaR(jj)":plots/${sample}/controlplots_drjj_dm2tlep.root \
	        "max#Sigmap_{T}(j)":plots/${sample}/controlplots_sumpt_dm2tlep.root \
	        "min|M(jj)-M(W)|":plots/${sample}/controlplots_mjj_dm2tlep.root
	    mkdir plots/${sample}/w
	    mv *.{png,pdf} plots/${sample}/w


	    python checkMatchingDistributions.py "min|m_{thad}-m_{tlep}|":plots/${sample}/controlplots_drjj_dm2tlep.root \
		"min|m_{thad}_m_{t}|":plots/${sample}/controlplots_sumpt_dm2pdg.root \
		"min#DeltaR(b,W)":plots/${sample}/controlplots_mjj_dr.root
	    mkdir plots/${sample}/top
	    mv *.{png,pdf} plots/${sample}/top
	done
	;;

    SEL)
	common="--wjjOrder drjj --thadOrder dm2tlep"
        for sample in MC8.16TeV_TTbar_pPb MC8.16TeV_TTbar_pPb_Pohweg MC8.16TeV_WJets_pPb; do        
	    python prepareWorkspace.py  -d ${sample} ${common};
	done
	;;

    SELDATA )
     	common="--wjjOrder drjj --thadOrder dm2tlep"
        for sample in Data8.16TeV_pPb_nonsubiso; do #Data8.16TeV_pPb
	    python prepareWorkspace.py  -d ${sample} ${common} --jerProf plots/MC8.16TeV_TTbar_pPb/controlplots.root;
        done
	;;
    
    PARAM)

        for sample in MC8.16TeV_TTbar_pPb MC8.16TeV_TTbar_pPb_Pohweg; do

	    python parameterizeMCShapes.py workspace_${sample}.root;
	    outDir=plots/$sample/workspace/
	    mkdir ${outDir}
	    b=(`ls *.{png,pdf}`)
	    for i in ${b[@]}; do
		f=${i/./_ttbar.}
		mv -v ${i} ${outDir}/${f};
	    done
	done
        
	python parameterizeMCShapes.py workspace_MC8.16TeV_WJets_pPb.root pdf_workspace_MC8.16TeV_TTbar_pPb.root;
	a=(`ls *.{png,pdf}`)
	for i in ${a[@]}; do
	    f=${i/./_w.}
	    mv -v ${i} plots/MC8.16TeV_WJets_pPb/${f};
	done
	;;

    FITS)
        sigRef=MC8.16TeV_TTbar_pPb #_Pohweg
        dataRef=Data8.16TeV_pPb_nonsubiso
	for f in `seq 0 0`; do
	    python runDataFit.py --fitType ${f} -s pdf_workspace_${sigRef}.root -i workspace_${dataRef}.root -o plots/${dataRef};
	    python runDataFit.py --fitType ${f} -o plots/${dataRef}/ --verbose 9 --finalWorkspace finalworkspace.root; #--impacts
	done
	;;

    FINALCOMP )
        python doPRLplots.py -i finalfitworkspace_0.root;
        for var in mjj mthad mtlep; do
            for cat in 1l4j2b 1l4j1b1q 1l4j2q; do
                python mergePLRplotCategs.py $var $cat;
            done
        done
	;;
esac
