import os

def getEraConfiguration(era,isData):
    
    """ defines global tags, JEC/R corrections, etc. depending on the era """

    globalTags = {
        'era2016':('94X_mcRun2_asymptotic_v3', '94X_dataRun2_v10'),
        'era2017':('94X_mc2017_realistic_v14', '94X_dataRun2_v6')
        }
    jecFiles    = {
        'era2016':('Summer16_07Aug2017_V11_MC',  'Summer16_07Aug2017All_V11_DATA'),
        'era2017':('Fall17_17Nov2017_V6_MC',     'Fall17_17Nov2017BCDEF_V6_DATA')
        }
    jerFiles    = {
        'era2016':('Summer16_25nsV1_MC',         'Summer16_25nsV1_DATA'),
        'era2017':('Summer16_25nsV1_MC',         'Summer16_25nsV1_DATA'),
        }

    globalTag = globalTags[era][isData]
    jecFile   = jecFiles[era][isData]
    jecTag    = '_'.join( jecFile.split('_')[0:-1] )
    jecDB     = 'jec_DATA.db'  if isData else 'jec_MC.db'
    jerFile   = jerFiles[era][isData]
    jerTag    = '_'.join( jerFile.split('_')[0:-1] )
    jerDB     = 'jer_DATA.db'  if isData else 'jer_MC.db'

    if not os.path.isfile(jecDB) and not os.path.islink(jecDB):
        print 'Creating symbolic link to JEC correction DB'
        os.system('ln -s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jecFile,jecDB))
    if not os.path.isfile(jerDB) and not os.path.islink(jerDB):
        print 'Creating symbolic link to JER correction DB'
        os.system('ln -s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jerFile,jerDB))
    print 'JEC tag: ',jecTag,'to be read from',jecDB
    print 'JER tag: ',jerTag,'to be read from',jerDB

    return globalTag, jecTag, jecDB, jerTag, jerDB
    
