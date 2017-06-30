#ifndef DiphotonAnalyzer_TreeProducer_FillNumberLUTHandler_h
#define DiphotonAnalyzer_TreeProducer_FillNumberLUTHandler_h

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

namespace CTPPSAlCa
{
  class FillNumberLUTHandler
  {
    public:
      typedef unsigned short fill_t;
      typedef std::map<unsigned short,fill_t> LUT;
    public:
      inline FillNumberLUTHandler( const char* file ):
        valid_( false )
      {
        std::ifstream f( file );
        if ( !f.is_open() ) return;
        std::string ss;
        while ( f >> ss ) {
          size_t pos = ss.find( ":" );
          if ( pos==std::string::npos ) continue;
          unsigned short run  = atoi( ss.substr( 0, pos ).c_str() ),
                         fill = atoi( ss.substr( pos+1, ss.size() ).c_str() );
          fillInRun_[run] = fill;
          if ( std::find( fills_.begin(), fills_.end(), fill )==fills_.end() ) { fills_.push_back( fill ); }
        }
        valid_ = true;
      }

      inline size_t getFillNumber( unsigned short run ) const {
        if ( !valid_ or !fillInRun_.count( run ) ) return 0;
        return fillInRun_.at( run );
      }
      inline LUT::const_iterator begin() const { return fillInRun_.begin(); }
      inline LUT::const_iterator end() const { return fillInRun_.end(); }

      std::vector<fill_t> fills() const { return fills_; }

      inline bool isValid() const { return valid_; }

    private:
      bool valid_;
      LUT fillInRun_;
      std::vector<fill_t> fills_;
  };
}

#endif
