import math
import ROOT

class Particle:
    def __init__(self, id, pt,eta,phi,m):
        self.id = id
        self.p4 = ROOT.TLorentzVector(0,0,0,0)
        self.p4.SetPtEtaPhiM(pt,eta,phi,m)
        self.mcTruth=None
        self.daughters=[]
        self.scaleUnc={}
    def setScaleUnc(self,name,unc):
        self.scaleUnc[name]=unc
    def setRank(self,rankVal):
        self.rankVal=rankVal
    def setMCtruth(self,part):
        self.mcTruth=part
    def DeltaR(self,part):
        return self.p4.DeltaR(part.p4)
    def addDaughter(self,part):
        self.daughters.append(part)
    def Print(self,pf=''):
        print '%s id=%d p=(%3.1f,%3.1f,%3.1f,%3.1f)'%(pf,self.id,self.p4.Pt(),self.p4.Eta(),self.p4.Phi(),self.p4.M())
        if len(self.daughters) :
            print '%s daughters list:'%pf
            for d in self.daughters: d.Print(pf=pf+'\t')
        if self.mcTruth:
            print '%s matched to '%pf
            self.mcTruth.Print(pf=pf+'\t')

"""
build possible dijet systems and ranks them by increasing dR or closest to MW
"""
def buildWjj(lightJets,orderBy='drjj'):
    WjjCandidates=[]
    for j in xrange(0,len(lightJets)):
        for k in xrange(j+1,len(lightJets)):
            chi=0
            if orderBy=='mjj':
                chi=ROOT.TMath.Abs((lightJets[j].p4+lightJets[k].p4).M()-80.4)
            elif orderBy=='sumpt':
                chi=1./(lightJets[j].p4.Pt()+lightJets[k].p4.Pt())
            else:
                chi=lightJets[j].DeltaR(lightJets[k])
            jj=lightJets[j].p4+lightJets[k].p4
            WjjCandidates.append( Particle(24,jj.Pt(),jj.Eta(),jj.Phi(),jj.M()) )
            WjjCandidates[-1].setRank(chi)
            WjjCandidates[-1].addDaughter( lightJets[j] )
            WjjCandidates[-1].addDaughter( lightJets[k] )

            for syst in ['jes','jer']:
                for var in ['up','dn']:
                    p1=ROOT.TLorentzVector(lightJets[j].p4)
                    scaleP4(p1,lightJets[j].scaleUnc[syst+var])
                    p2=ROOT.TLorentzVector(lightJets[k].p4)
                    scaleP4(p2,lightJets[k].scaleUnc[syst+var])
                    unc=(p1+p2).M()/jj.M()
                    WjjCandidates[-1].setScaleUnc(syst+var,unc)

    WjjCandidates.sort(key=lambda x: x.rankVal,reverse=False)
    return WjjCandidates

"""
Impose the w mass to solve the neutrino kinematics
"""
def buildWlnu(lepton,MET):
    M_W  = 80.4;
    M_mu =  0.10566;
    emu  = lepton.p4.E();
    pxmu = lepton.p4.Px();
    pymu = lepton.p4.Py();
    pzmu = lepton.p4.Pz();
    pxnu = MET.p4.Px();
    pynu = MET.p4.Py();
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

    #neutrino
    p4nu=ROOT.TLorentzVector(pxnu,pynu,pznu,math.sqrt(pxnu**2++pynu**2+pznu**2))
    nu=Particle(12,p4nu.Pt(),p4nu.Eta(),p4nu.Phi(),0.)

    #W->lnu
    wlp4=lepton.p4+nu.p4
    Wlnu=Particle(24,wlp4.Pt(),wlp4.Eta(),wlp4.Phi(),wlp4.M())
    Wlnu.addDaughter(lepton)
    Wlnu.addDaughter(nu)

    #all done here
    return Wlnu

"""
"""
def buildTTbar(bJets,Wjj,Wlnu):
    ttbarCandidates=[]
    for i in xrange(0,2):
        idx1,idx2=i%2,(i+1)%2
        tp4=bJets[idx1].p4+Wjj.p4
        tbarp4=bJets[idx2].p4+Wlnu.p4
        ttbarCandidates.append( (Particle(6,tp4.Pt(),   tp4.Eta(),   tp4.Phi(),    tp4.M()),
                                 Particle(6,tbarp4.Pt(),tbarp4.Eta(),tbarp4.Phi(), tbarp4.M()) ) )
        ttbarCandidates[-1][0].addDaughter(bJets[idx1])
        ttbarCandidates[-1][0].addDaughter(Wjj)
        ttbarCandidates[-1][1].addDaughter(bJets[idx2])
        ttbarCandidates[-1][1].addDaughter(Wlnu)
    return ttbarCandidates

"""
"""
def randomlyRotate(part):

    #choose randomly one daughter particle to rotate
    idx=int(ROOT.gRandom.Uniform(0,len(part.daughters)-1))
    p4=part.daughters[idx].p4
    rotP4=None
    while True:
        rotP4=ROOT.TLorentzVector(p4)
        en,pabs,phi,theta=rotP4.E(),rotP4.P(),ROOT.gRandom.Uniform(0,2*ROOT.TMath.Pi()),ROOT.TMath.ACos(ROOT.gRandom.Uniform(-1,1))
        rotP4.SetPxPyPzE(pabs*ROOT.TMath.Cos(phi)*ROOT.TMath.Sin(theta),
                         pabs*ROOT.TMath.Sin(phi)*ROOT.TMath.Sin(theta),
                         pabs*ROOT.TMath.Cos(theta),
                         en)

        if abs(rotP4.Eta())>2.4 or rotP4.Pt()<20 : continue
        break

    #add the p4 with the p4 of other daughters
    for i in xrange(0,len(part.daughters)):
        if i==idx: continue
        rotP4+=part.daughters[i].p4

    #all done
    return rotP4

def scaleP4(p4,scale):
    p4.SetPxPyPzE(p4.Px()*scale,p4.Py()*scale,p4.Pz()*scale,p4.E()*scale)
