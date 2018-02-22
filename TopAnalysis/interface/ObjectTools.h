#ifndef _object_tools_h_
#define _object_tools_h_

#include <vector>
#include <algorithm>

#include "TLorentzVector.h"

class Particle {
  
  public:
    Particle(TLorentzVector p4, int charge, int id, int qualityFlags, int origRef, double puppi = 1)
      : p4_(p4), charge_(charge), id_(id), qualityFlags_(qualityFlags), origRef_(origRef), puppi_(puppi) {}
   
    Particle( const Particle &p) 
      : p4_(p.p4_), charge_(p.charge_), id_(p.id_), qualityFlags_(p.qualityFlags_), origRef_(p.origRef_), puppi_(p.puppi_) {}

    double px()     { return p4_.Px();  }
    double py()     { return p4_.Py();  }
    double pz()     { return p4_.Pz();  }
    double e()      { return p4_.E();   }
    double pt()     { return p4_.Pt();  }
    double eta()    { return p4_.Eta(); }
    double phi()    { return p4_.Phi(); }
    double energy() { return p4_.E(); }
    double m()      { return p4_.M(); }
    double mass()   { return p4_.M(); }
    TLorentzVector p4()       { return p4_; }
    TLorentzVector momentum() { return p4_; }
    int charge()    { return charge_; }
    int id()        { return id_; }
    int qualityFlags() { return qualityFlags_; }
    bool hasQualityFlag(int bit) { return ((qualityFlags_>>bit)&0x1); }
    int originalReference() { return origRef_; }
    void setOriginalReference(int origRef) { origRef_=origRef; }
    double puppi()  { return puppi_; }

  private:
    TLorentzVector p4_;
    int charge_, id_, qualityFlags_,origRef_;
    double puppi_;
};

/**
   @short summarizes the information on a jet needed for the charmed meson analysis
 */
typedef std::pair<TLorentzVector,int> IdTrack;

class Jet {

  public:
    Jet(TLorentzVector p4, int flavor, int idx)
      : p4_(p4), flavor_(flavor), idx_(idx), overlap_(0) {}
    Jet(TLorentzVector p4, float csv, int idx)
      : p4_(p4), csv_(csv), idx_(idx) {}
    Jet(const Jet &j) 
      : p4_(j.p4_), particles_(j.particles_), trks_(j.trks_), csv_(j.csv_), flavor_(j.flavor_), idx_(j.idx_), overlap_(j.overlap_) {}
    ~Jet() {}

    double pt()     { return p4_.Pt();  }
    double eta()    { return p4_.Eta();  }
    TLorentzVector p4()       { return p4_; }
    TLorentzVector momentum() { return p4_; }
    std::vector<Particle> &particles() { return particles_; }
    int flavor()  { return flavor_; }
    int overlap() { return overlap_; }
    
    void addParticle(Particle p) { particles_.push_back(p); }
    void setFlavor(int flavor)   { flavor_ = flavor; }
    void setOverlap(int overlap) { overlap_ = overlap; }
    
    void addTrack(TLorentzVector p4, int pfid) { trks_.push_back( IdTrack(p4,pfid) ); }
    TLorentzVector &getVec() { return p4_; }
    float &getCSV() { return csv_; }
    void setCSV(float csv) { csv_=csv; }
    int &getJetIndex() { return idx_; }
    std::vector<IdTrack> &getTracks() { return trks_; }
    void sortTracksByPt() { sort(trks_.begin(),trks_.end(), sortIdTracksByPt); }
    
    static bool sortJetsByPt(Jet i, Jet j)  { return i.getVec().Pt() > j.getVec().Pt(); }
    static bool sortJetsByCSV(Jet i, Jet j) { return i.getCSV() > j.getCSV(); }
  
  private:
    static bool sortIdTracksByPt(IdTrack i, IdTrack j)  { return i.first.Pt() > j.first.Pt(); }
    
    TLorentzVector p4_;
    std::vector<Particle> particles_;
    std::vector<IdTrack> trks_;
    float csv_;
    int flavor_;
    int idx_;
    int overlap_;
};
#endif
