#include "TopLJets2015/TopAnalysis/interface/AlignmentLUTHandler.h"

namespace CTPPSAlCa
{
  AlignmentLUTHandler::AlignmentLUTHandler( const char* file ) :
    valid_( false ),
    rgx_hdr_ ( std::regex( "\\[\\w+\\/fill_(\\d+)\\/.+\\]" ) ),
    rgx_algn_( std::regex( "id=(\\d+),sh_x=([0-9.+-]+)(?:,sh_x_unc=([0-9.+-]+),sh_y=([0-9.+-]+),sh_y_unc=([0-9.+-]+))?" ) )
  {
    // load the alignment constants from external file
    loadConstants( file );
  }

  RPAlignmentConstants&
  AlignmentLUTHandler::retrieveConstants( const fill_t& fill_num )
  {
    LUT::iterator it = align_map_.find( fill_num );
    if ( it!=align_map_.end() ) return it->second;
    throw cms::Exception( "NoFillInfo" ) << "Failed to retrieve the alignment parameters for fill " << fill_num;
  }

  const RPAlignmentConstants
  AlignmentLUTHandler::getAlignmentConstants( const fill_t& fill_num ) const
  {
    if ( !valid_ ) throw cms::Exception( "InvalidLUT" ) << "Invalid look-up table";
    LUT::const_iterator it = align_map_.find( fill_num );
    if ( it!=align_map_.end() ) return it->second;

    if ( fill_num<align_map_.begin()->first ) {
      return align_map_.begin()->second; // retrieve the alignment from the first fill in the LUT
    }
    if ( fill_num>align_map_.rbegin()->first ) {
      return align_map_.rbegin()->second; // retrieve the alignment from the last fill in the LUT
    }
    return RPAlignmentConstants();
  }

  void
  AlignmentLUTHandler::loadConstants( const char* file )
  {
    std::ifstream f( file );
    if ( !f.is_open() ) return;
    std::string ss;
    fill_t fill_num = 0;

    std::smatch match;

    while ( f >> ss ) {
      if ( std::regex_match( ss, match, rgx_hdr_ ) and match.size()==2 ) { // new fill information
        fill_num = atoi( match[1].str().c_str() );
        align_map_.insert( std::make_pair( fill_num, RPAlignmentConstants() ) );
      }
      else if ( std::regex_match( ss, match, rgx_algn_ ) and match.size()>0 ) { // new alignment parameters
        RPAlignmentConstants& align = retrieveConstants( fill_num );
        const unsigned short rp_id = atoi( match[1].str().c_str() );

        const float align_x = ( match.size()>2 ) ? atof( match[2].str().c_str() ) : 0.;
        const float err_align_x = ( match.size()>3 ) ? atof( match[3].str().c_str() ) : 0.;
        const float align_y = ( match.size()>4 ) ? atof( match[4].str().c_str() ) : 0.;
        const float err_align_y = ( match.size()>5 ) ? atof( match[5].str().c_str() ) : 0.;

        align.setQuantities( rp_id, RPAlignmentConstants::Quantities( align_x, align_y, err_align_x, err_align_y ) );
      }
      else {
        throw cms::Exception( "InvalidLine" ) << "Failed to process line:\n\t" << ss;
        return;
      }
    }
    valid_ = true;
  }
}
