import os

def getEraConfiguration(era,isData):
    globalTags = {
        'era2016':('80X_mcRun2_asymptotic_2016_TrancheIV_v7', '80X_dataRun2_2016SeptRepro_v6'),
        'era2017':('94X_mc2017_realistic_v10',                '94X_dataRun2_ReReco_EOY17_v2')
        }
    jecTags    = {
        'era2016':('Summer16_23Sep2016V4',  'Summer16_23Sep2016AllV4'),
        'era2017':('Fall17_17Nov2017_V4',   'Fall17_17Nov2017_V4')
        }

    idx=1 if isData else 0
    globalTag = globalTags[era][idx]
    jecTag = jecTags[era][idx]
    jecTag_pf = 'DATA' if isData else 'MC'
    #fixme once there are JEC for 2017 data
    if era=='era2017': jecTag_pf='MC'
    jecDB  = 'jec_DATA.db'  if isData else 'jec_MC.db'

    if not os.path.isfile(jecDB) and not os.path.islink(jecDB):
        print 'Creating symbolic link to DB correction DB'
        os.system('ln -s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s_%s.db %s'%(era,jecTag,jecTag_pf,jecDB))

    return globalTag, jecTag, jecDB
    
