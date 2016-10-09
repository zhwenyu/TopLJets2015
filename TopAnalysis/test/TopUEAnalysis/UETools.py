#!/usr/bin/env/python

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
    if idx==0   : return 'tow'
    elif idx==1 : return 'tra'
    return 'awa'


"""
parses the event and counts particles in each region at gen/rec levels
"""
class UEEventCounter:
    def __init__(self,t):

        self.rec_nch   = [0]*3
        self.rec_ptsum = [0]*3
        self.rec_avgpt = [0]*3
        
        self.nch   = [[0]*3,[0]*3,[0]*3]
        self.ptsum = [[0]*3,[0]*3,[0]*3]
        self.avgpt = [[0]*3,[0]*3,[0]*3]

        self.gen_nch   = [0]*3
        self.gen_ptsum = [0]*3
        self.gen_avgpt = [0]*3
                
        #reco level
        passSel=(t.passSel&0x1)
        if passSel:                    
            for n in xrange(0,t.n):
                isInB=(t.isInBFlags[n] & 0x1)
                if isInB : continue
                dphi_rec=ROOT.TMath.Abs( ROOT.TVector2.Phi_mpi_pi( t.phi[n]-t.rec_phi_ttbar[0] ) )
                idx_rec=getRegionFor(dphi_rec)
                dphi_gen=ROOT.TMath.Abs( ROOT.TVector2.Phi_mpi_pi( t.phi[n]-t.gen_phi_ttbar ) )
                idx_gen=getRegionFor(dphi_gen)
                self.rec_nch[idx_rec]        +=1
                self.nch[idx_gen][idx_rec]   +=1
                self.rec_ptsum[idx_rec]      += t.pt[n]
                self.ptsum[idx_gen][idx_rec] += t.pt[n]
            for idx_rec in xrange(0,3):
                self.rec_avgpt[idx_rec] = self.rec_ptsum[idx_rec]/self.rec_nch[idx_rec] if self.rec_nch[idx_rec]>0 else 0
                for idx_gen in xrange(0,3):
                    self.avgpt[idx_gen][idx_rec] = self.ptsum[idx_gen][idx_rec]/self.nch[idx_gen][idx_rec] if self.nch[idx_gen][idx_rec]>0 else 0

        #gen level
        passSel=(t.gen_passSel&0x1)
        if passSel:
            for n in xrange(0,t.gen_n):
                dphi_gen=ROOT.TMath.Abs( ROOT.TVector2.Phi_mpi_pi( t.gen_phi[n]-t.gen_phi_ttbar ) )
                idx_gen=getRegionFor(dphi_gen)
                self.gen_nch[idx_gen]   +=1
                self.gen_ptsum[idx_gen] += t.gen_pt[n]
            for idx_gen in xrange(0,3):
                self.gen_avgpt[idx_gen] = self.gen_ptsum[idx_gen]/self.gen_nch[idx_gen] if self.gen_nch[idx_gen]>0 else 0.


"""
return the most appropriate bin for a given value, taking into account the range available
"""
def getBinForVariable(h,val,axis):
    xmin,xmax=axis.GetXmin(),axis.GetXmax()
    if val>xmax : return axis.GetNbins()
    if val<xmin : return 0
    return axis.FindBin(val)

"""
get purity and stability for a given bin configuration
"""
def getPurityStability(h,x1,x2,y1,y2):
    NgenrecUnc=ROOT.Double(0)
    Ngenrec = h.IntegralAndError(x1, x2,            y1, y2,NgenrecUnc)    
    Nrec    = h.Integral        (0,  h.GetNbinsX(), y1, y2)
    Ngen    = h.Integral        (x1, x2,            0,  h.GetNbinsY())    
    
    pur=(Ngenrec/Ngen,NgenrecUnc/Ngen) if Ngen!=0 else (0,0)
    stab=(Ngenrec/Nrec,NgenrecUnc/Nrec) if Nrec!=0 else (0,0)
    return pur,stab


"""
project purity and stability
"""
def getPurityStabilityGraphs(h):
    purGr=ROOT.TGraphErrors()
    purGr.SetName(h.GetName()+'_pur')
    purGr.SetMarkerStyle(20)
    purGr.SetTitle('purity')
    stabGr=purGr.Clone(h.GetName()+'_stab')
    stabGr.SetMarkerStyle(24)
    stabGr.SetTitle('stability')
    for xbin in xrange(1,h.GetNbinsX()+1):
        y1 = getBinForVariable(h, h.GetXaxis().GetBinLowEdge(xbin),  h.GetYaxis())
        y2 = getBinForVariable(h, h.GetXaxis().GetBinUpEdge(xbin+1), h.GetYaxis())

        xcen     = h.GetXaxis().GetBinCenter(xbin)
        xwid     = h.GetXaxis().GetBinWidth(xbin)
        pur,stab = getPurityStability(h,xbin,xbin+1,y1,y2)

        np = purGr.GetN()
        purGr.SetPoint(np,xcen,pur[0])
        purGr.SetPointError(np,0.5*xwid,pur[1])
        stabGr.SetPoint(np,xcen,stab[0])
        stabGr.SetPointError(np,0.5*xwid,stab[1])
    return purGr,stabGr


"""
optimize migration matrix based on eff/purity criteria
"""
def optimizeMigrationMatrix(h,minStab=0.5, minPur=0.5):

    #save current binning
    xlimits,ylimits=[],[]
    for xbin in xrange(1,h.GetNbinsX()) : xlimits.append(h.GetXaxis().GetBinLowEdge(xbin))
    xlimits.append(h.GetXaxis().GetXmax())
    for ybin in xrange(1,h.GetNbinsY()) : ylimits.append(h.GetYaxis().GetBinLowEdge(ybin))
    ylimits.append(h.GetYaxis().GetXmax())

    #start by optimizing gen binning based on purity criteria
    x1,x2=1,1
    new_xlimits=[h.GetXaxis().GetXmin()]
    while x2<=h.GetNbinsX():
        pur=0
        while pur<minPur and x2<=h.GetNbinsX():
            y1=getBinForVariable(h,h.GetXaxis().GetBinLowEdge(x1),h.GetYaxis())
            y2=getBinForVariable(h,h.GetXaxis().GetBinUpEdge(x2), h.GetYaxis())
            ipur,_=getPurityStability(h,x1,x2,y1,y2)
            pur=ipur[0]
            if pur<minPur   : x2+=1
        x2val=h.GetXaxis().GetBinUpEdge(x2)
        if x2val!=new_xlimits[-1] : new_xlimits.append(x2val)
        x1,x2=x2+1,x2+1
    xmax=h.GetXaxis().GetXmax()
    if new_xlimits[-1]!=xmax: new_xlimits[-1]=xmax

    #define intermediate histo
    h_xopt=ROOT.TH2F('h_xopt','',len(new_xlimits)-1,array.array('d',new_xlimits),len(ylimits)-1,array.array('d',ylimits))
    for xbin in xrange(1,h.GetNbinsX()+1):
        for ybin in xrange(1,h.GetNbinsY()+1):
            h_xopt.Fill( h.GetXaxis().GetBinCenter(xbin), h.GetYaxis().GetBinCenter(ybin), h.GetBinContent(xbin,ybin) )
    
    #optimize the rec binning based on the stability criteria
    y1,y2=1,1
    new_ylimits=[h_xopt.GetYaxis().GetXmin()]
    while y2<=h_xopt.GetNbinsY():
        stab=0
        while stab<minStab and y2<=h_xopt.GetNbinsY():
            x1=getBinForVariable(h_xopt, h_xopt.GetYaxis().GetBinLowEdge(y1), h_xopt.GetXaxis())
            x2=getBinForVariable(h_xopt, h_xopt.GetYaxis().GetBinUpEdge(y2),  h_xopt.GetXaxis())
            _,istab=getPurityStability(h_xopt,x1,x2,y1,y2)
            stab=istab[0]
            if stab<minStab   : y2+=1            
        y2val=h_xopt.GetYaxis().GetBinUpEdge(y2)
        if y2val!=new_ylimits[-1] : new_ylimits.append(y2val)
        y1,y2=y2+1,y2+1
    ymax=h.GetYaxis().GetXmax()
    if new_ylimits[-1]!=ymax: new_ylimits[-1]=ymax

    #define final histo
    h_yopt=ROOT.TH2F('h_yopt','',len(new_xlimits)-1,array.array('d',new_xlimits),len(new_ylimits)-1,array.array('d',new_ylimits))
    for xbin in xrange(1,h.GetNbinsX()+1):
        for ybin in xrange(1,h.GetNbinsY()+1):
            h_yopt.Fill( h.GetXaxis().GetBinCenter(xbin), h.GetYaxis().GetBinCenter(ybin), h.GetBinContent(xbin,ybin) )

    c=ROOT.TCanvas('c','c',1200,800)
    c.Divide(3,2)
    c.cd(1)
    h.Draw('colz')
    c.cd(2)
    h_xopt.Draw('colz')
    c.cd(3)
    h_yopt.Draw('colz')
    c.cd(4)
    drawOpt='ap'
    hngr=getPurityStabilityGraphs(h)
    for gr in hngr:
        gr.DrawClone(drawOpt)
        drawOpt='p'
    print hngr
    c.cd(5)
    drawOpt='ap'
    hxgr=getPurityStabilityGraphs(h_xopt)
    for gr in hxgr:
        gr.Draw(drawOpt)
        drawOpt='p'
    c.cd(6)
    drawOpt='ap'
    hygr=getPurityStabilityGraphs(h_yopt)
    for gr in hygr:
        gr.Draw(drawOpt)
        drawOpt='p'
    
    raw_input()
