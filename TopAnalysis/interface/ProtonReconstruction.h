#ifndef TopLJets2015_TopAnalysis_ProtonReconstruction_h
#define TopLJets2015_TopAnalysis_ProtonReconstruction_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/Provenance/interface/RunID.h"

#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"

#include <iostream>

namespace CTPPSAlCa
{
  class RPAlignmentConstants
  {
    public:
      class Quantities
      {
        public:
          Quantities() : x( 0. ), err_x( 0. ), y( 0. ), err_y( 0. ) {}
          Quantities( float x, float y, float err_x=0., float err_y=0. ) : x( x ), err_x( err_x ), y( y ), err_y( err_y ) {}

          bool operator==( const Quantities& rhs ) const { return ( x==rhs.x && err_x==rhs.err_x && y==rhs.y && err_y==rhs.err_y ); }
          void operator=( const Quantities& rhs ) { x = rhs.x; err_x = rhs.err_x; y = rhs.y; err_y = rhs.err_y; }
          friend std::ostream& operator<<( std::ostream&, const Quantities& );

        public:
          float x, err_x, y, err_y;
      };
      typedef std::map<unsigned short, Quantities> Map;

    public:
      RPAlignmentConstants();
      friend std::ostream& operator<<( std::ostream&, const RPAlignmentConstants& );
      bool operator==( const RPAlignmentConstants& ) const;
      void operator=( const RPAlignmentConstants& );

      Map::const_iterator begin() const { return map_.begin(); }
      Map::const_iterator end() const { return map_.end(); }

      void setQuantities( unsigned short detid, const Quantities& quant );
      const Quantities quantities( unsigned short detid ) const;

    private:
      Map map_;
  };
}

class XiInterpolator
{
  public:
    XiInterpolator();
    XiInterpolator( const char* filename );
    ~XiInterpolator();

    void loadInterpolationGraphs( const char* filename, bool conversion=false );
    void setAlignmentConstants( const CTPPSAlCa::RPAlignmentConstants& ac ) { align_ = ac; }

    void computeXiSpline( unsigned int arm, unsigned int pot, float track_x, float& xi, float& err_xi );

  private:
    std::shared_ptr<TSpline3> isLF_, isLN_, isRF_, isRN_;
    CTPPSAlCa::RPAlignmentConstants align_;
};

#endif
