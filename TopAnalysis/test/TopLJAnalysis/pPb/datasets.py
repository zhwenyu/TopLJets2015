#!/usr/bin/env python2.7

from collections import OrderedDict

def defineDataset(data):

    if 'Data8.16TeV_pPb' in data:
        dataset=OrderedDict()
        if '_e' in data:
            dataset['l+jets']         = (['data/Data_EJets_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'],       '#4575b4', [11],   False)
            dataset['non-iso l+jets'] = (['data/DataNonIso_EJets.root'], '#91bfdb', [1100], False)
        if '_mu' in data:
            dataset['l+jets']         = (['data/Data_MuJets_RecheckWithMergedFiles.root'],       '#4575b4', [13],   False)   #data/Data_MuJets.root
            dataset['non-iso l+jets'] = (['data/DataNonIso_MuJets.root'], '#91bfdb', [1300], False)
    if 'MC8.16TeV_TTbar_pPb' in data:
        dataset=OrderedDict()
        
        dataset['l+jets']         = (['data/PYQUEN_TTJets_PPb_MuJets.root','data/PYQUEN_TTJets_PbP_MuJets.root',
                                      'data/PYQUEN_TTJets_PPb_EJets.root','data/PYQUEN_TTJets_PbP_EJets.root'],
                                      '#4575b4', [11,13],   True)
        if '_e' in data:
            dataset['l+jets']=(dataset['l+jets'][0][2:4],dataset['l+jets'][1],[11],True)
        if '_mu' in data:
            dataset['l+jets']=(dataset['l+jets'][0][0:2],dataset['l+jets'][1],[13],True)


    if 'MC8.16TeV_WJets_pPb' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['data/PYQUEN_WJets_PbP_MuJets.root','data/PYQUEN_WJets_PPb_MuJets.root'],
                                    '#4575b4', [11,13],   False)

    if 'MC8TeV_WJets_pp' in data:
        dataset=OrderedDict()
        dataset['l+jets'] = ( [#'data_pp/MC8TeV_WJets.root',
                               #'data_pp/MC8TeV_W1Jets.root',
                               #'data_pp/MC8TeV_W2Jets.root',
                               'data_pp/MC8TeV_W3Jets.root',
                               ],
                               '#4575b4',
                               [11,13],#[11,13],
                               False )
    if data=='Data8TeV_pp':
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['data_pp/Data8TeV_SingleMu2012A.root'
                                #'data_pp/Data8TeV_SingleMu2012B.root',
                                #'data_pp/Data8TeV_SingleMu2012C.root',
                                #'data_pp/Data8TeV_SingleMu2012D.root'
                                ],
                                 '#4575b4',
                                 [11,13],
                                 False )
        dataset['non-iso l+jets'] = ( dataset['l+jets'][0],'#fc8d59',[11000,13000], False )

    if data=='MC8TeV_TTJets':
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['data_pp/MC8TeV_TTJets_MSDecays_172v5.root'],
                                 '#4575b4',
                                 [13], #[11,13]
                                 True )
        #dataset['non-iso l+jets'] = ( dataset['l+jets'][0], '#91bfdb', [11000,13000],False )

        #dataset['non-iso l+jets'] = ( dataset['l+jets'][0], '#91bfdb', [11000,13000],False )

    return dataset
