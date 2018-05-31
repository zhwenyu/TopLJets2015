#!/usr/bin/env python

import ROOT
import os
import numpy as np
from array import array

baseDir='/eos/cms/store/cmst3/group/top/psilva/9bea816/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_MC13TeV_TTJets/180519_130226/0000/'
t=ROOT.TChain('analysis/data')
for x in os.listdir(baseDir):
    t.AddFile( os.path.join(baseDir,x) )


uncList={11:[],13:[]}
for i in xrange(0,t.GetEntries()):
    t.GetEntry(i)
    for l in xrange(0,t.nl):
        pt=t.l_pt[l]
        aeta=abs(t.l_eta[l])
        pdgId=abs(t.l_id[l])

        if pt<20 or aeta>2.4 : 
            continue

        pid=t.l_pid[l]
        relIso=t.l_relIso[l]
        passId=False
        if abs(pdgId)==11:
            passId=True if ((pid>>7)&0x1)==1 else False
        else:
            passId=True if ((pid>>4) &0x1)==1 else False
            if relIso>0.15: passId=False

        if not passId: continue
        uncList[pdgId].append([pt,aeta,t.l_scaleUnc_1[l]/pt,t.l_scaleUnc_2[l]/pt,t.l_scaleUnc_2[l]/pt])

scaleUncH={}
for pdgId in uncList:        
    vals=np.array(uncList[pdgId])
    nq=[10,20,40,60,80,99]
    q=np.percentile(vals,nq, axis=0)

    ptRange=np.insert(q[:,0], 0, 20., axis=0)
    etaRange=np.insert(q[:,1],0, 0., axis=0)

    for i in xrange(0,3):
        scaleUncH[(pdgId,i)]=ROOT.TH2F('scaleUnc_%d_%d'%(pdgId,i),
                                       'scaleUnc_%d_%d'%(pdgId,i),
                                       len(ptRange)-1,  array('d',ptRange),
                                       len(etaRange)-1, array('d',etaRange))
        
    for j in xrange(0,len(ptRange)-1):
        for k in xrange(0,len(etaRange)-1):
            ptMin=ptRange[j]
            ptMax=ptRange[j+1]
            etaMin=etaRange[k]
            etaMax=etaRange[k+1]
            cVals=vals[vals[:,0]>ptMin]
            cVals=cVals[cVals[:,0]<ptMax]
            cVals=cVals[cVals[:,1]>etaMin]
            cVals=cVals[cVals[:,1]<etaMax]
            if len(cVals)==0: continue
            avg=np.average(cVals,axis=0)    
            for i in xrange(0,3):
                scaleUncH[(pdgId,i)].SetBinContent(j+1,k+1,avg[i+2])

fOut=ROOT.TFile.Open('ParametrizedLeptonScaleUnc.root','recreate')
for key in scaleUncH:
    scaleUncH[key].SetDirectory(fOut)
    scaleUncH[key].Write()
fOut.Close()


            
