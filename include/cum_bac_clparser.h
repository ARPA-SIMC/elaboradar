#ifndef ARCHIVIATORE_CUM_BAC_CL_PARSER
#define ARCHIVIATORE_CUM_BAC_CL_PARSER

#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

struct CUM_BAC_CLOPT {
   std::string filename;
   int filetype;
   std::string sito;
   bool do_medium;
   bool do_quality;
   bool do_beamblocking;
   bool do_declut;
   bool do_bloccor;
   bool do_vpr;
   bool do_clean;
   bool do_class;
   bool do_devel;
   bool do_anaprop;

   bool data_in_odim;
   bool do_readStaticMap;  
   bool do_intermediateProd;
   bool do_SaveBothRanges; 
   bool do_SaveFullRes;
}  ;


void parseOptions(int argc, char** argv, struct CUM_BAC_CLOPT *opt);
//CUM_BAC_CLOPT parseOptions(int argc, char** argv);
void PrintOptions(struct CUM_BAC_CLOPT *opt);

#endif
