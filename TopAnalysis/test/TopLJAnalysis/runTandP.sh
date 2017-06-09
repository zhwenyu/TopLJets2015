#!/bin/bash                                                                                                                                                                                                
WHAT=$1;

case $WHAT in
    FITS )
        f=(BWResCBExp) #CB2Exp)
        d=(True False)

        for coll in pPb_noembed; do #pPb
            for i in ${f[@]}; do 
                for j in ${d[@]}; do                        
                    cmsRun runTandP_HLT.py      data=${j} func=${i} tree=${coll}_8.16TeV & 
	            cmsRun runTandP_IDISO.py    data=${j} func=${i} tree=${coll}_8.16TeV &  
	            cmsRun runTandP_HLTIDISO.py data=${j} func=${i} tree=${coll}_8.16TeV & 
                done       
            done
        done
        ;;
    
    PLOTS )
        
        for coll in pPb; do
            #python drawTandPResults.py hltfits_${coll}_8.16TeV_1_BWResCBExp.root:"Data",hltfits_${coll}_noembed_8.16TeV_0_BWResCBExp.root:"MC (no embed.)",hltfits_${coll}_8.16TeV_0_BWResCBExp.root:"MC"           ~/www/LJets-${coll}/HLT_mu
            #python drawTandPResults.py idisofits_${coll}_8.16TeV_1_BWResCBExp.root:"Data",idisofits_${coll}_noembed_8.16TeV_0_BWResCBExp.root:"MC (no embed)",idisofits_${coll}_8.16TeV_0_BWResCBExp.root:"MC"       ~/www/LJets-${coll}/IDISO_mu
            #python drawTandPResults.py hltidisofits_${coll}_8.16TeV_1_BWResCBExp.root:"Data",hltidisofits_${coll}_noembed_8.16TeV_0_BWResCBExp.root:"MC (no embed)",hltidisofits_${coll}_8.16TeV_0_BWResCBExp.root:"MC" ~/www/LJets-${coll}/HLTIDISO_mu

            #python drawTandPResults.py hltfits_${coll}_noembed_8.16TeV_0_BWResCBExp.root:'Data',hltfits_${coll}_noembed_8.16TeV_1_BWResCBExp.root:'MC (no embed.)'  ~/www/LJets-${coll}/HLT_mu;
            python drawTandPResults.py idisofits_${coll}_noembed_8.16TeV_1_BWResCBExp.root:"Data",idisofits_${coll}_noembed_8.16TeV_0_BWResCBExp.root:"MC (no embed.)"       ~/www/LJets-${coll}/IDISO_mu;
            python drawTandPResults.py hltidisofits_${coll}_noembed_8.16TeV_1_BWResCBExp.root:"Data",hltidisofits_${coll}_noembed_8.16TeV_0_BWResCBExp.root:"MC (no embed.)" ~/www/LJets-${coll}/HLTIDISO_mu;
        done
        
        ;;
esac