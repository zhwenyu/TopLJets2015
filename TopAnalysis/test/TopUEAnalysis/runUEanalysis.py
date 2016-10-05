#!/usr/bin/env/python

import glob
import sys
import os
import optparse
import ROOT
import numpy as np
import array as array

"""
0 - tow(ards), 1 trans(verse), 2 away
"""
def getRegionFor(dphi) :
    if ROOT.TMath.Abs(dphi) < ROOT.TMath.Pi()/3.     : return 0
    elif ROOT.TMath.Abs(dphi) < 2*ROOT.TMath.Pi()/3. : return 1
    return 2

"""
converts index to name
"""
def getRegionName(idx):
    if idx==0 : return 'tow'
    elif idx==1: return 'trans'
    return 'away'

"""
return the most appropriate bin for a given value, taking into account the range available
"""
def getBinForVariable(h,val):
    xmin,xmax=h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
    if val>xmax : return h.GetNbinsX()
    if val<xmin : return 0
    return h.GetXaxis().FindBin(val)

"""
opens analysis files and reads out the histograms defining the base binning
"""
def readSlicesAndObservables(url):
    
    #readout the bins to use
    sliceVarsQ,obsVarsQ={},{}
    fIn=ROOT.TFile.Open(url)
    for k in fIn.Get('quantiles').GetListOfKeys():
        if 'obs_' in k.GetName():
            obsVarsQ[k.GetName()] = fIn.Get('quantiles/%s'%k.GetName())
            obsVarsQ[k.GetName()].SetDirectory(0)
        else:
            sliceVarsQ[k.GetName()] = fIn.Get('quantiles/%s'%k.GetName())
            sliceVarsQ[k.GetName()].SetDirectory(0)
    fIn.Close()
    return  sliceVarsQ,obsVarsQ

"""
defines the variables to slice the UE measurement and saves their reco quantiles
"""
def defineBaseProjectionBins(opt):

    #buidl the chain
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): t.AddFile(f)
    
    #define variables and quantiles to apply
    sliceVars={ 'pt_ttbar':[], 'nch':[], 'mll':[], 'dphill':[], 'nj':[]}
    sliceVarsQDef={
        'pt_ttbar' : [100.*x/6.  for x in xrange(0,7)],
        'nch'      : [100.*x/12. for x in xrange(0,13)],
        'mll'      : [100.*x/12. for x in xrange(0,13) ],
        'dphill'   : [100.*x/12. for x in xrange(0,13) ],
        }
    obsVars={'nch':[], 'ptsum':[], 'avgpt':[]}

    #loop the available events
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()
            
        if (t.passSel&0x1)==0 : continue
        
        #count particles
        nch,ptsum=0,0
        for n in xrange(0,t.n):
            if (t.isInBFlags[n] & 0x1) : continue
            nch   +=1
            ptsum += t.pt[n]
        avgpt=ptsum/nch if nch>0 else 0.

        sliceVars['pt_ttbar'].append(t.rec_pt_ttbar[0])
        sliceVars['nch'].append(nch)
        sliceVars['mll'].append(t.mll)
        sliceVars['nj'].append(t.nj[0]-2)
        sliceVars['dphill'].append(t.dphill)
        obsVars['nch'].append(nch)
        obsVars['ptsum'].append(ptsum)
        obsVars['avgpt'].append(avgpt)

    #determine quantiles
    sliceVarsQ={'nj':[0,1,2,3,4]}
    for s in sliceVars:
        if s=='nj' : continue
        sliceVarsQ[s]=[]      
        for q in sliceVarsQDef[s]:
            sliceVarsQ[s].append( np.percentile(sliceVars[s], q) )
    obsVarsQ={}
    for s in obsVars:
        obsVarsQ[s]=[]
        for q in [100.*x/10. for x in xrange(0,11)]:
            obsVarsQ[s].append( np.percentile(obsVars[s], q) )
            
    #save to ROOT file
    fOut=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'RECREATE')
    outDir=fOut.mkdir('quantiles')
    outDir.cd()
    for s in sliceVarsQ:
        h=ROOT.TH1F('slice_%s'%s,';%s;Events/bin'%s,
                    len(sliceVarsQ[s])-1, 
                    array.array('d',sliceVarsQ[s])
                    )
        h.SetDirectory(outDir)
        h.Sumw2()
        for val in sliceVars[s] : 
            bin=h.GetXaxis().FindBin(val)
            h.Fill( val, 1./h.GetXaxis().GetBinWidth(bin) )
        h.Write()
    for s in obsVarsQ:
        h=ROOT.TH1F('obs_%s'%s,';%s;Events/bin'%s,
                    len(obsVarsQ[s])-1, 
                    array.array('d',obsVarsQ[s])
                    )
        h.SetDirectory(outDir)
        h.Sumw2()
        for val in obsVars[s] : 
            bin=h.GetXaxis().FindBin(val)
            h.Fill( val, 1./h.GetXaxis().GetBinWidth(bin) )
        h.Write()
    fOut.Close()

"""
loops over a set of files with common name to fill the migration matrices
"""        
def fillMigrationMatrices(opt,varIdx=0,wgtIdx=0):

    sliceVarsQ,obsVarsQ=readSlicesAndObservables(url='%s/UEanalysis.root'%opt.out)

    #build the migration matrices
    migMatrix={}
    for s in sliceVarsQ:
        nRecBins=sliceVarsQ[s].GetXaxis().GetNbins()
        nGenBins=nRecBins if 'nj' in s else nRecBins/2
        for i in xrange(0,nGenBins+2):
            for j in xrange(0,nRecBins+2):
                for o in obsVarsQ:
                    if o==s : continue

                    obsNRecBins=1 if j==0 else obsVarsQ[o].GetXaxis().GetNbins()
                    obsNGenBins=1 if i==0 else obsVarsQ[o].GetXaxis().GetNbins()/2

                    for recCat in ['tow','trans','away']:
                        for genCat in ['tow','trans','away']:
                            if obsNRecBins==1 : recCat=''
                            if obsNGenBins==1 : genCat=''
                            name='%s%d%d_%s%s%s'%(s,i,j,o,genCat,recCat)
                            key=(s,i,j,o,genCat,recCat)
                            if key in migMatrix : continue
                            migMatrix[key]=ROOT.TH2F(name,
                                                     '%s;gen %s;rec %s'%(name,o,o),
                                                     obsNGenBins,0,obsNGenBins,
                                                     obsNRecBins,0,obsNRecBins)

    #build the chain
    t=ROOT.TChain('tue')
    for f in glob.glob(opt.input):
        t.AddFile(f)

    #loop over the available events
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #measure at GEN level
        gen_passSel=t.gen_passSel
        gen_nch,gen_ptsum,gen_avgpt=[0,0,0],[0,0,0],[0,0,0]
        if gen_passSel:
            for n in xrange(0,t.gen_n):
                dphi_gen=ROOT.TMath.Abs( ROOT.TVector2.Phi_mpi_pi( t.gen_phi[n]-t.gen_phi_ttbar ) )
                idx_gen=getRegionFor(dphi_gen)
                gen_nch[idx_gen]+=1
                gen_ptsum[idx_gen]+=t.gen_pt[n]
            gen_avgpt = [ gen_ptsum[k]/gen_nch[k] if gen_nch[k]>0. else 0. for k in xrange(0,len(gen_nch)) ]


        #measure at RECO level
        passSel=((t.passSel>>varIdx) & 0x1)
        nch,ptsum,avgpt = [[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]
        nchPerReg,ptsumPerReg,avgPtPerReg=[0,0,0],[0,0,0],[0,0,0]
        if passSel:
            for n in xrange(0,t.n):
                isInB = ((t.isInBFlags[n]>>varIdx) & 0x1)
                if isInB : continue
                dphi_rec=ROOT.TMath.Abs( ROOT.TVector2.Phi_mpi_pi( t.phi[n]-t.rec_phi_ttbar[varIdx] ) )
                idx_rec=getRegionFor(dphi_rec)

                dphi_gen=ROOT.TMath.Abs( ROOT.TVector2.Phi_mpi_pi( t.phi[n]-t.gen_phi_ttbar ) )
                idx_gen=getRegionFor(dphi_gen)

                nch[idx_gen][idx_rec]+=1
                nchPerReg[idx_rec]+=1
                ptsum[idx_gen][idx_rec]+=t.pt[n]
                ptsumPerReg[idx_rec]+=t.pt[n]
            for l in xrange(0,3):
                for k in xrange(0,3):
                    avgpt[l][k]  = ptsum[l][k]/nch[l][k] if nch[l][k]>0. else 0.
                avgPtPerReg[l]=ptsumPerReg[l]/nchPerReg[l] if nchPerReg[l]>0. else 0.
                
        #determine bins for the SLICE variable
        pt_ttbar_recBin  = getBinForVariable(sliceVarsQ['slice_pt_ttbar'],  t.rec_pt_ttbar[varIdx])
        nj_recBin        = getBinForVariable(sliceVarsQ['slice_nj'],        ROOT.TMath.Max(0,t.nj[varIdx]-2))
        nch_recBin       = getBinForVariable(sliceVarsQ['slice_nch'],       np.sum(nch))
        mll_recBin       = getBinForVariable(sliceVarsQ['slice_mll'],       t.mll)
        dphill_recBin    = getBinForVariable(sliceVarsQ['slice_dphill'],    t.dphill)
        if not passSel:
            pt_ttbar_recBin=0
            nj_recBin=0
            nch_recBin=0
            mll_recBin=0
            dphill_recBin=0
        pt_ttbar_genBin  = getBinForVariable(sliceVarsQ['slice_pt_ttbar'], t.gen_pt_ttbar)/2
        nj_genBin        = getBinForVariable(sliceVarsQ['slice_nj'],       t.gen_nj)
        nch_genBin       = getBinForVariable(sliceVarsQ['slice_nch'],      np.sum(gen_nch))/2
        mll_genBin       = getBinForVariable(sliceVarsQ['slice_mll'],      t.gen_mll)/2
        dphill_genBin    = getBinForVariable(sliceVarsQ['slice_dphill'],   t.gen_dphill)/2
        if not gen_passSel :
            pt_ttbar_genBin=0            
            nj_genBin=0
            nch_genBin=0
            mll_genBin=0
            dphill_genBin=0
        else :
            pt_ttbar_genBin += 1
            #as the number of gen and reco bins for nj is the same, no need to do nj_genBin += 1
            nch_genBin+=1
            mll_genBin+=1
            dphill_genBin+=1

        #fill the migration matrices
        #(notice the bins returned are always subtracted by -1 to correspond to the bin lower edge)
        for l in xrange(0,3):            
            if passSel and gen_passSel:
                for k in xrange(0,3):
                                        
                    obs_nch_recBin=getBinForVariable(obsVarsQ['obs_nch'],nch[l][k])-1
                    obs_nch_genBin=getBinForVariable(obsVarsQ['obs_nch'],gen_nch[l])/2
                    key=('slice_pt_ttbar',pt_ttbar_genBin,pt_ttbar_recBin,'obs_nch',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(obs_nch_genBin,obs_nch_recBin,t.weight[wgtIdx])                    
                    key=('slice_nj',nj_genBin,nj_recBin,'obs_nch',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(obs_nch_genBin,obs_nch_recBin,t.weight[wgtIdx])
                    key=('slice_mll',mll_genBin,mll_recBin,'obs_nch',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(obs_nch_genBin,obs_nch_recBin,t.weight[wgtIdx])
                    key=('slice_dphill',dphill_genBin,dphill_recBin,'obs_nch',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(obs_nch_genBin,obs_nch_recBin,t.weight[wgtIdx])
                                        
                    ptsum_recBin=getBinForVariable(obsVarsQ['obs_ptsum'],ptsum[l][k])-1
                    ptsum_genBin=getBinForVariable(obsVarsQ['obs_ptsum'],gen_ptsum[l])/2
                    key=('slice_pt_ttbar',pt_ttbar_genBin,pt_ttbar_recBin,'obs_ptsum',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(ptsum_genBin,ptsum_recBin,t.weight[wgtIdx])                   
                    key=('slice_nj',nj_genBin,nj_recBin,'obs_ptsum',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(ptsum_genBin,ptsum_recBin,t.weight[wgtIdx])
                    key=('slice_nch',nch_genBin,nch_recBin,'obs_ptsum',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(ptsum_genBin,ptsum_recBin,t.weight[wgtIdx])
                    key=('slice_mll',mll_genBin,mll_recBin,'obs_ptsum',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(ptsum_genBin,ptsum_recBin,t.weight[wgtIdx])
                    key=('slice_dphill',dphill_genBin,dphill_recBin,'obs_ptsum',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(ptsum_genBin,ptsum_recBin,t.weight[wgtIdx])
                    
                    avgpt_recBin=getBinForVariable(obsVarsQ['obs_ptsum'],avgpt[l][k])-1
                    avgpt_genBin=getBinForVariable(obsVarsQ['obs_ptsum'],gen_avgpt[l])/2
                    key=('slice_pt_ttbar',pt_ttbar_genBin,pt_ttbar_recBin,'obs_avgpt',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(avgpt_genBin,avgpt_recBin,t.weight[wgtIdx])
                    key=('slice_nj',nj_genBin,nj_recBin,'obs_avgpt',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(avgpt_genBin,avgpt_recBin,t.weight[wgtIdx])
                    key=('slice_nch',nch_genBin,nch_recBin,'obs_avgpt',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(avgpt_genBin,avgpt_recBin,t.weight[wgtIdx])
                    key=('slice_mll',mll_genBin,mll_recBin,'obs_avgpt',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(avgpt_genBin,avgpt_recBin,t.weight[wgtIdx])
                    key=('slice_dphill',dphill_genBin,dphill_recBin,'obs_avgpt',getRegionName(l),getRegionName(k))
                    migMatrix[key].Fill(avgpt_genBin,avgpt_recBin,t.weight[wgtIdx])
            else:
                keys=[]
                if passSel and not gen_passSel:
                    binVal=(0,getBinForVariable(obsVarsQ['obs_nch'],nchPerReg[l])-1)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_nch','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_nch','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_nch','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_nch','',getRegionName(l))] )
                    binVal=(0,getBinForVariable(obsVarsQ['obs_ptsum'],ptsumPerReg[l])-1)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_ptsum','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_ptsum','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_nch',      nch_genBin,      nch_recBin,      'obs_ptsum','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_ptsum','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_ptsum','',getRegionName(l))] )
                    binVal=(0,getBinForVariable(obsVarsQ['obs_avgpt'],avgPtPerReg[l])-1)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_avgpt','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_avgpt','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_nch',      nch_genBin,      nch_recBin,      'obs_avgpt','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_avgpt','',getRegionName(l))] )
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_avgpt','',getRegionName(l))] )
                elif not passSel and gen_passSel:
                    binVal=(getBinForVariable(obsVarsQ['obs_nch'],gen_nch[l])/2,0)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_nch',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_nch',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_nch',getRegionName(l),'')] )                   
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_nch',getRegionName(l),'')] )                   
                    binVal=(getBinForVariable(obsVarsQ['obs_ptsum'],gen_ptsum[l])/2,0)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_ptsum',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_ptsum',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_nch',      nch_genBin,      nch_recBin,      'obs_ptsum',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_ptsum',getRegionName(l),'')] )                   
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_ptsum',getRegionName(l),'')] )                   
                    binVal=(getBinForVariable(obsVarsQ['obs_avgpt'],gen_avgpt[l])/2,0)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_avgpt',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_avgpt',getRegionName(l),'')] )
                    keys.append( [ binVal,('slice_nch',      nch_genBin,      nch_recBin,      'obs_avgpt',getRegionName(l),'')] )                   
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_avgpt',getRegionName(l),'')] )                   
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_avgpt',getRegionName(l),'')] )                    
                else :
                    binVal=(0,0)
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_nch','','')] )
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_nch','','')] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_nch','','')] )
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_nch','','')] )
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_ptsum','','')] )                    
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_ptsum','','')] )
                    keys.append( [ binVal,('slice_nch',      nch_genBin,      nch_recBin,      'obs_ptsum','','')] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_ptsum','','')] )
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_ptsum','','')] )
                    keys.append( [ binVal,('slice_pt_ttbar', pt_ttbar_genBin, pt_ttbar_recBin, 'obs_avgpt','','')] )                   
                    keys.append( [ binVal,('slice_nj',       nj_genBin,       nj_recBin,       'obs_avgpt','','')] )
                    keys.append( [ binVal,('slice_nch',      nch_genBin,      nch_recBin,      'obs_avgpt','','')] )
                    keys.append( [ binVal,('slice_mll',      mll_genBin,      mll_recBin,      'obs_avgpt','','')] )
                    keys.append( [ binVal,('slice_dphill',   dphill_genBin,   dphill_recBin,   'obs_avgpt','','')] )
                for binVal,key in keys:
                    migMatrix[key].Fill(binVal[0],binVal[1],t.weight[wgtIdx])

    #save to ROOT file
    fOut=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'UPDATE')
    outDirName='migmatrices_%d_%d'%(varIdx,wgtIdx)
    outDir=fOut.Get(outDirName)
    try:
        outDir.cd()
    except:
        outDir=fOut.mkdir(outDirName)
        outDir.cd()
    for k in migMatrix:
        migMatrix[k].SetDirectory(outDir)
        migMatrix[k].Write(migMatrix[k].GetName(),ROOT.TObject.kOverwrite)
    fOut.Close()


"""
"""
def showMatrices(opt,varIdx=0,wgtIdx=0):

    #prepare canvas
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(55)
    c=ROOT.TCanvas('c','c',1000,1000)
    c.SetTopMargin(0.01)
    c.SetRightMargin(0.1)
    c.SetBottomMargin(0.15)

    
    #open ROOT file with the definitions
    url='%s/UEanalysis.root'%opt.out
    sliceVarsQ,obsVarsQ=readSlicesAndObservables(url=url)
    fIn=ROOT.TFile.Open(url)
    mmDir='migmatrices_%d_%d'%(varIdx,wgtIdx)

    #build the fully combined migration matrices
    bigMatrix={}
    for s in sliceVarsQ:
        nRecBins=sliceVarsQ[s].GetXaxis().GetNbins()
        nGenBins=nRecBins if 'nj' in s else nRecBins/2
        for o in obsVarsQ:
            if o==s : continue
            obsNRecBins=obsVarsQ[o].GetXaxis().GetNbins()
            obsNGenBins=obsVarsQ[o].GetXaxis().GetNbins()/2

            key=(s,o)
            nbinsX=1+3*nGenBins*obsNGenBins
            nbinsY=1+3*nRecBins*obsNRecBins
            bigMatrix[key]=ROOT.TH2F('%s_%s'%(s,o),'',nbinsX,0,nbinsX,nbinsY,0,nbinsY)
            bigMatrix[key].SetDirectory(0)

            #fill the matrix
            genBinCtr=0
            for xbin in xrange(1,nbinsX+1):
                if xbin!=1 :
                    genBinCtr+=1
                    if genBinCtr>obsNGenBins: genBinCtr=1

                recBinCtr=0
                for ybin in xrange(1,nbinsY+1):
                    if ybin!=1:
                        recBinCtr+=1
                        if recBinCtr>obsNRecBins: recBinCtr=1                

                    sliceGenBin=(xbin-2)/(3*obsNGenBins)+1
                    obsGenRegionBin=((xbin-2)/obsNGenBins)%3
                    obsGenRegion=getRegionName(obsGenRegionBin)

                    sliceRecBin=(ybin-2)/(3*obsNRecBins)+1
                    obsRecRegionBin=((ybin-2)/obsNRecBins)%3
                    obsRecRegion=getRegionName(obsRecRegionBin)   

                    if xbin==1 and ybin==1:
                        mmH=fIn.Get('%s/%s00_%s'%(mmDir,s,o))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(1,1))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(1,1))
                    elif xbin==1:                                                
                        mmH=fIn.Get('%s/%s0%d_%s%s'%(mmDir,s,sliceRecBin,o,obsRecRegion))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(1,recBinCtr))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(1,recBinCtr))
                    elif ybin==1:
                        mmH=fIn.Get('%s/%s%d0_%s%s'%(mmDir,s,sliceGenBin,o,obsGenRegion))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(genBinCtr,1))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(genBinCtr,1))
                    else:
                        mmH=fIn.Get('%s/%s%d%d_%s%s%s'%(mmDir,s,sliceGenBin,sliceRecBin,o,obsGenRegion,obsRecRegion))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(genBinCtr,recBinCtr))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(genBinCtr,recBinCtr))
                        
                        
            #label x-axis
            sliceBin,obsBin=0,0
            for xbin in xrange(1,nbinsX+1):
                label=''
                if xbin==1: label='fail'
                else:
                    if (xbin+1-obsNGenBins/2)%(3*obsNGenBins)==0:
                        label='#color[38]{[%3.1f,%3.1f[}'%(sliceVarsQ[s].GetXaxis().GetBinLowEdge(2*sliceBin+1),sliceVarsQ[s].GetXaxis().GetBinUpEdge(2*sliceBin+2))
                        sliceBin+=1
                    elif (xbin+4)%(obsNGenBins)==0:
                        obsBin+=1
                        label='#it{%s}'%getRegionName((obsBin-1)%3)
                if len(label)==0: continue    
                bigMatrix[key].GetXaxis().SetBinLabel(xbin,label)

            #label y-axis
            sliceBin,obsBin=0,0
            for ybin in xrange(1,nbinsY+1):
                label=''
                if ybin==1: label='fail'
                else:
                    if (ybin+1-obsNRecBins/2)%(3*obsNRecBins)==0:
                        label='#color[38]{[%3.1f,%3.1f[}'%(sliceVarsQ[s].GetXaxis().GetBinLowEdge(sliceBin+1),sliceVarsQ[s].GetXaxis().GetBinUpEdge(sliceBin+1))
                        sliceBin+=1
                    elif (ybin+4)%(obsNRecBins)==0:
                        obsBin+=1
                        label='#it{%s}'%getRegionName((obsBin-1)%3)
                if len(label)==0: continue    
                bigMatrix[key].GetYaxis().SetBinLabel(ybin,label)
            
            #display the matrix
            c.Clear()
            c.SetLogz()
            bigMatrix[key].Draw('colz')
            bigMatrix[key].GetZaxis().SetLabelSize(0.03)
            bigMatrix[key].GetYaxis().SetLabelSize(0.03)
            bigMatrix[key].GetXaxis().SetLabelSize(0.03)
            bigMatrix[key].GetXaxis().SetTickLength(0)
            bigMatrix[key].GetYaxis().SetTickLength(0)
            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.025)
            tex.SetNDC()
            tex.DrawLatex(0.01,0.96,'#bf{#splitline{Reco.}{level}}')
            tex.DrawLatex(0.95,0.1,'#bf{#splitline{Gen.}{level}}')
            stit='p_{T}(t#bar{t})'
            if 'phi' in s : stit='#phi(t#bar{t})'
            if 'nj' in s : stit='N(jets)'
            if 'nch' in s : stit='N(ch)'                
            otit='N(ch)'
            if 'ptsum' in o : otit='#Sigmap_{T}(ch)'
            if 'avgpt' in o : otit='#bar{p}_{T}(ch)'
            tex.DrawLatex(0.10,0.03,'#scale[1.2]{#bf{CMS}} #it{simulation preliminary} %s vs %s'%(stit,otit))

            line=ROOT.TLine()           
            line.SetLineWidth(1)
            for xbin in xrange(1,nbinsX+1):
                if xbin==1 or (xbin-1)%(obsNGenBins)==0:
                    ls=1 if (xbin-1)%(3*obsNGenBins)==0 else 9
                    lc=1 if (xbin-1)%(3*obsNGenBins)==0 else ROOT.kGray
                    line.SetLineStyle(ls)
                    line.SetLineColor(lc)  
                    line.DrawLine(bigMatrix[key].GetXaxis().GetBinUpEdge(xbin),bigMatrix[key].GetYaxis().GetXmin(),bigMatrix[key].GetXaxis().GetBinUpEdge(xbin),bigMatrix[key].GetYaxis().GetXmax())
            for ybin in xrange(1,nbinsY+1):
                if ybin==1 or (ybin-1)%(obsNRecBins)==0:
                    ls=1 if (ybin-1)%(3*obsNRecBins)==0 else 9
                    lc=1 if (ybin-1)%(3*obsNRecBins)==0 else ROOT.kGray
                    line.SetLineStyle(ls)
                    line.SetLineColor(lc)                      
                    line.DrawLine(bigMatrix[key].GetXaxis().GetXmin(),bigMatrix[key].GetYaxis().GetBinUpEdge(ybin),bigMatrix[key].GetXaxis().GetXmax(),bigMatrix[key].GetYaxis().GetBinUpEdge(ybin))

                    
            c.Modified()
            c.Update()
            raw_input()

    #save to ROOT file
    fOut=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'UPDATE')
    outDir=fOut.Get('final_'+mmDir)
    try:
        outDir.cd()
    except:
        outDir=fOut.mkdir('final_'+mmDir)
        outDir.cd()
    for k in bigMatrix:
        bigMatrix[k].SetDirectory(outDir)
        bigMatrix[k].Write(bigMatrix[k].GetName(),ROOT.TObject.kOverwrite)
    fOut.Close()

def main():

    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',       dest='input',    help='input',               default='MC13TeV_TTJets_dilpowheg_0.root',   type='string')
    parser.add_option('-s', '--step',     dest='step',     help='step',                default=1,   type=int)
    parser.add_option('-o', '--out',      dest='out',      help='output',              default='./UEanalysis',   type='string')
    (opt, args) = parser.parse_args()

    os.system('mkdir -p %s'%opt.out)

    if opt.step==1 or opt.step<0:
        defineBaseProjectionBins(opt)
    if opt.step==2 or opt.step<0:
        fillMigrationMatrices(opt)
    if opt.step==3 or opt.step<0:
        showMatrices(opt)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
