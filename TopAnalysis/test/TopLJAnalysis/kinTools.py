import math
import ROOT

"""
build possible dijet systems and ranks them by increasing dR or closest to MW
"""
def buildDijets(lightJets,byMass=True):
    dijets=[]
    for j in xrange(0,len(lightJets)):
        for k in xrange(j+1,len(lightJets)):
            chi=0
            if byMass:
                chi=ROOT.TMath.Abs((lightJets[j]+lightJets[k]).M()-80.4)
            else:
                chi=lightJets[j].DeltaR(lightJets[k])
            dijets.append( (lightJets[j]+lightJets[k], chi , (j,k) ) )
    dijets.sort(key=lambda x: x[1],reverse=False)
    return dijets

"""
builds possible bjj triplets ranking them by increasing dR or closest to Mtop
"""
def buildTrijets(j1,j2,bJets,byMass=True):
    trijets=[]
    wjj=j1+j2
    for j in xrange(0,len(bJets)):
        chi=0
        if byMass:
            chi=ROOT.TMath.Abs((wjj+bJets[j]).M()-172.5)
        else:
            chi=wjj.DeltaR(bJets[j])
        trijets.append( (wjj+bJets[j],chi,j) )
    trijets.sort(key=lambda x: x[1],reverse=False)
    return trijets

"""
rescale the energy of a 4-momentum vector
"""
def scaleP4(vec,scale):
    vec=ROOT.TLorentzVector(vec.Px()*scale,
                            vec.Py()*scale,
                            vec.Pz()*scale,
                            vec.E()*scale)
    return vec

"""
Impose the w mass to solve the neutrino kinematics
"""
def getNeutrinoP4(lepton,MET):
    M_W  = 80.4;
    M_mu =  0.10566;
    emu = lepton.E();
    pxmu = lepton.Px();
    pymu = lepton.Py();
    pzmu = lepton.Pz();
    pxnu = MET.Px();
    pynu = MET.Py();
    pznu = 0.;

    a = M_W*M_W - M_mu*M_mu + 2.0*pxmu*pxnu + 2.0*pymu*pynu;
    A = 4.0*(emu*emu - pzmu*pzmu);
    B = -4.0*a*pzmu;
    C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;
    tmproot = B*B - 4.0*A*C
    if tmproot<0:
        pznu = - B/(2*A);
    else:
        tmpsol1 = (-B + math.sqrt(tmproot))/(2.0*A);
        tmpsol2 = (-B - math.sqrt(tmproot))/(2.0*A);

        # two real roots, pick the one closest to pz of muon
        if abs(tmpsol2-pzmu) < abs(tmpsol1-pzmu) :
            pznu = tmpsol2
        else :
            pznu = tmpsol1
        if pznu > 300. :
            if abs(tmpsol1)<abs(tmpsol2) :
                pznu = tmpsol1
            else :
                pznu = tmpsol2

    p4nu=ROOT.TLorentzVector(pxnu,pynu,pznu,math.sqrt(pxnu**2++pynu**2+pznu**2))
    return p4nu
