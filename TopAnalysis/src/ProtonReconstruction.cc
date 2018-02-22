#include "TopLJets2015/TopAnalysis/interface/ProtonReconstruction.h"

XiInterpolator::XiInterpolator()
{}

XiInterpolator::XiInterpolator( const char* filename )
{
  loadInterpolationGraphs( filename );
}

XiInterpolator::~XiInterpolator()
{}

void
XiInterpolator::loadInterpolationGraphs( const char* filename, bool conversion )
{
  TFile f( filename );
  if ( !f.IsOpen() ) {
    edm::LogError("XiInterpolator") << "Failed to load the interpolation graphs file";
    return;
  }

  // already "processed" curves
  isRN_ = std::make_shared<TSpline3>( *dynamic_cast<TSpline3*>( f.Get( "s_x_to_xi_R_1_N" )->Clone() ) );
  isRF_ = std::make_shared<TSpline3>( *dynamic_cast<TSpline3*>( f.Get( "s_x_to_xi_R_1_F" )->Clone() ) );
  isLN_ = std::make_shared<TSpline3>( *dynamic_cast<TSpline3*>( f.Get( "s_x_to_xi_L_1_N" )->Clone() ) );
  isLF_ = std::make_shared<TSpline3>( *dynamic_cast<TSpline3*>( f.Get( "s_x_to_xi_L_1_F" )->Clone() ) );

  edm::LogInfo("XiInterpolator") << "Interpolation graphs successfully loaded from file " << filename;
}

void
XiInterpolator::computeXiSpline( unsigned int arm, unsigned int pot, float track_x, float& xi, float& err_xi )
{
  xi = err_xi = 0.;

  // retrieve the alignment parameters
  const CTPPSAlCa::RPAlignmentConstants::Quantities ac = align_.quantities( arm*100+pot );
  //std::cout << "--> for this pot:\n" << ac << std::endl;

  // retrieve the proper interpolation curve
  TSpline3 *interp = 0;
  switch ( arm ) { // 0 = sector 45, 1 = sector 56
    case 0: {
      if ( pot==2 ) interp = isLN_.get();
      if ( pot==3 ) interp = isLF_.get();
    } break;
    case 1: {
      if ( pot==2 ) interp = isRN_.get();
      if ( pot==3 ) interp = isRF_.get();
    } break;
    default: return;
  }
  if ( !interp ) return;

  const float de_x = 0.4e-3, // m
              de_rel_dx = 0.1;

  // apply the alignment
  const float x_corr = track_x + ac.x*1.e-3; // convert to m

  xi = interp->Eval( x_corr );

  const float de_xi = interp->Eval( x_corr + de_x ) - xi;
  err_xi = std::sqrt( std::pow( de_xi, 2 )
                    + std::pow( de_rel_dx * xi, 2 ) );
}

namespace CTPPSAlCa
{
  RPAlignmentConstants::RPAlignmentConstants()
  {
    map_[2] = Quantities();
    map_[3] = Quantities();
    map_[102] = Quantities();
    map_[103] = Quantities();
  }

  bool
  RPAlignmentConstants::operator==( const RPAlignmentConstants& rhs ) const
  {
    return ( map_==rhs.map_ );
  }

  void
  RPAlignmentConstants::operator=( const RPAlignmentConstants& rhs )
  {
    map_ = rhs.map_;
  }

  void
  RPAlignmentConstants::setQuantities( unsigned short detid, const Quantities& quant )
  {
    map_[detid] = quant;
  }

  const RPAlignmentConstants::Quantities
  RPAlignmentConstants::quantities( unsigned short detid ) const
  {
    Map::const_iterator it = map_.find( detid );
    if ( it!=map_.end() ) return it->second;
    edm::LogWarning("RPAlignmentConstants") << "Failed to retrieve quantities for RP with DetId " << detid;
    return Quantities();
  }

  std::ostream&
  operator<<( std::ostream& os, const RPAlignmentConstants& align )
  {
    for ( RPAlignmentConstants::Map::const_iterator it=align.map_.begin(); it!=align.map_.end(); ++it ) {
      os << "DetId = " << it->first << ": " << it->second << std::endl;
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const RPAlignmentConstants::Quantities& quant )
  {
    return os << "x-shift: " << quant.x << " +/- " << quant.err_x << "\t"
              << "y-shift: " << quant.y << " +/- " << quant.err_y;
  }
}
