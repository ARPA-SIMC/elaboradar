/*
 * =====================================================================================
 *
 *       Filename:  cum_bac_clparser.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  26/02/2014 09:12:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "cum_bac_clparser.h"
#include <iostream>
void parseOptions(int argc, char** argv , CUM_BAC_CLOPT *opt)
{
    //    CUM_BAC_CLOPT opt;
  try {


	TCLAP::CmdLine cmd("CUM_BAC ", ' ', "0.1" );
	//
	// Define required arguments
	//
//Filename
	TCLAP::UnlabeledValueArg<std::string> filename( "FileName", "Name of data file", true, "NULL", "std::string");
        cmd.add(filename);

//file type
	TCLAP::UnlabeledValueArg<int> filetype("Filetype", "Type of data file", true,9999, "int");
        cmd.add(filetype);

//sito
	TCLAP::UnlabeledValueArg<std::string> sito("Sito", "Name of radar site", true,"NULL", "std::string");
        cmd.add(sito);

	//
	// Define optional arguments
	//
//	TCLAP::SwitchArg Short("S", "ShortPulse", "Use parameters for short pulse", false);
	TCLAP::SwitchArg Medium("M", "MediumPulse", "Use parameters for medium pulse", false);
        cmd.add(Medium);
//	std::vector<TCLAP::Arg*> xorpulse;
//	xorpulse.push_back(&Short);
//	xorpulse.push_back(&Medium);
//	cmd.xorAdd( xorpulse );

	TCLAP::SwitchArg Quality("Q", "Quality", "Calculate quality", false);
        cmd.add(Quality);
	TCLAP::SwitchArg BeamBlocking("B", "BeamBlocking", "Performe BeamBlocking correction", false);
        cmd.add(BeamBlocking);
	TCLAP::SwitchArg Declut("D", "Declut", "Performe static declutter correction", false);
        cmd.add(Declut);
	TCLAP::SwitchArg BlocNoCor("b", "BlocNoCor", "Questo non so cosa sia", false);
        cmd.add(BlocNoCor);
	TCLAP::SwitchArg VPR("V", "VPR", "Performe VPR correction", false);
        cmd.add(VPR);
	TCLAP::SwitchArg Clean("c", "Clean", "Performe Wimax-second trip cleaning", false);
        cmd.add(Clean);
	TCLAP::SwitchArg Class("C", "Class", "Evaluate stratiform-convective classification", false);
        cmd.add(Class);
	TCLAP::SwitchArg Stampe("s", "StampeExtra", "Stampe Extra per devel", false);
        cmd.add(Stampe);

	TCLAP::SwitchArg StaticMap("m", "UseStaticMap", "Use Static Map", false);
        cmd.add(Stampe);
	
	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);

//
// check variables
//
	
	opt->filename=filename.getValue();
	opt->filetype=filetype.getValue();
	opt->sito=sito.getValue();
//	if (Short.isSet()) opt->do_medium=false;
	opt->do_medium=Medium.getValue();
	opt->do_quality=Quality.getValue();
	opt->do_beamblocking=BeamBlocking.getValue();
	opt->do_declut=Declut.getValue();
	opt->do_bloccor=BlocNoCor.getValue();
	opt->do_vpr=VPR.getValue();
	opt->do_clean=Clean.getValue();
	opt->do_class=Class.getValue();
	opt->do_devel=Stampe.getValue();
	opt->do_readStaticMap=StaticMap.getValue();  

  } catch ( TCLAP::ArgException& e )
  { std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl; }

  return;
}


void PrintOptions(struct CUM_BAC_CLOPT *opt){

   std::cout <<"Filename :"<<opt->filename.c_str()<<std::endl;
   std::cout <<"Filetype :"<<opt->filetype<<std::endl;
   std::cout <<"Sito :"<<opt->sito.c_str()<<std::endl;
   std::cout <<"do_medium :"<<opt->do_medium<<std::endl;
   std::cout <<"do_quality :"<<opt->do_quality<<std::endl;
   std::cout <<"do_beamblocking :"<<opt->do_beamblocking<<std::endl;
   std::cout <<"do_declutter :"<<opt->do_declut<<std::endl;
   std::cout <<"do_bloccor"<<opt->do_bloccor<<std::endl;
   std::cout <<"do_vpr"<<opt->do_vpr<<std::endl;
   std::cout <<"do_clean"<<opt->do_clean<<std::endl;
   std::cout <<"do_class :"<<opt->do_class<<std::endl;
   std::cout <<"do_devel :"<<opt->do_devel<<std::endl;
   std::cout <<"do_readStaticMap"<<opt->do_readStaticMap<<std::endl;

   return ;
}


