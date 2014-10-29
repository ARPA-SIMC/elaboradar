
/*----------------------------------------------------------------------------*/
  /*    INCLUDE file                                      */
/*----------------------------------------------------------------------------*/

#include "cum_bac.h"
#include "logging.h"
#include "config.h"
#include "cum_bac_clparser.h"

#include <memory>
#include <cstdlib>
#include <cmath>

#include <setwork.h>

using namespace std;

/*----------------------------------------------------------------------------*/
/*      FINE SEZIONE    INCLUDE                   */
/*----------------------------------------------------------------------------*/
// DEFINIZIONI PREPROCESSORE NB: VPR E CLASS IMPLICANO QUALITY ,

/*----------------------------------------------------------------------------*/
/*  DICHIARATIVE GLOBALI                              */
/*----------------------------------------------------------------------------*/

  /*non cancellare*/
  /* int trovo_top() */
  /* { */
  /*   int ier,i,k,l; */
  /*   FILE *file0; */

  /*   for(i=0; i<NUM_AZ_X_PPI; i++){ */
  /*     for (k=0; k<volume.vol_pol[0][i].ray.size(); k++){ */
  /*       for (l=first_level_static[i][k]; l<NEL; l++)     */
  /*    {  */
  /*      if (BYTEtoDB(volume.vol_pol[l][i].ray[k]) > 10. ) */
  /*        top[i][k]=(unsigned char)(quota_f((float)(volume.vol_pol[l][i].teta_true)*CONV_RAD+0.45*DTOR,k)/100.);{ //top in ettometri */
  /*        if (l >= NEL -1 ) top[i][k]=0; */
  /*      } */
  /*    } */
  /*     } */
  /*   } */
  /*   return ier; */
  /* } */

/**
 * Write a banner with program information, compile flags, environment and so on.
 */
static void startup_banner(CUM_BAC_CLOPT *opt)
{
    LOG_CATEGORY("radar.banner");
#ifdef WRITE_DBP
//            " WRITE_DBP"   eventualmente da re-inserire
#endif
#ifdef WRITE_DBP_REORDER
//            " WRITE_DBP_REORDER" eventiualemnte da re-inserire
#endif

#ifdef Z_AVERAGE
//            " Z_AVERAGE"  obsoleta da togliere
#endif
#ifdef R_AVERAGE
//            " R_AVERAGE"  obsoleta da togliere
#endif

    /// Log initial program state

    LOG_INFO("Lancio Programma");
    LOG_INFO("-----------------------------------------------------------------");
    std::string FlagRunTime ="Flag di RunTime: ";

    FlagRunTime=FlagRunTime+" "+opt->sito;
    if(opt->do_declut)FlagRunTime=FlagRunTime+" DECLUTTER";
    if(!opt->do_medium)FlagRunTime=FlagRunTime+" SHORT";
    if(opt->do_medium)FlagRunTime=FlagRunTime+" MEDIUM";
    if(opt->do_beamblocking)FlagRunTime=FlagRunTime+" BEAMBLOCKING";
    if(opt->do_quality)FlagRunTime=FlagRunTime+" QUALITY";
    if(opt->do_readStaticMap) FlagRunTime=FlagRunTime + " STATIC";  
 
    LOG_INFO(FlagRunTime.c_str());

    LOG_INFO("-----------------------------------------------------------------");
    LOG_INFO("Variabili d'Ambiente:");
    LOG_INFO("LISTA_FILE = %s", getenv("LISTA_FILE"));
    LOG_INFO("LAST_FILE = %s", getenv("LAST_FILE"));
    LOG_INFO("ANAP_STAT_FILE = %s", getenv("ANAP_STAT_FILE"));
    LOG_INFO("BLOC_STAT_FILE = %s", getenv("BLOC_STAT_FILE"));
    LOG_INFO("ELEV_STAT_FILE = %s", getenv("ELEV_STAT_FILE"));
    LOG_INFO("DIR_OUT_PP_BLOC = %s", getenv("DIR_OUT_PP_BLOC"));
    LOG_INFO("FILE_DEM_SPC = %s", getenv("FILE_DEM_SPC"));
    LOG_INFO("FILE_DEM_GAT = %s", getenv("FILE_DEM_GAT"));
    LOG_INFO("DIR_QUALITY = %s", getenv("DIR_QUALITY"));
    LOG_INFO("FIRST_LEVEL_FILE = %s", getenv("FIRST_LEVEL_FILE"));
    LOG_INFO("OUTPUT_Z_DIR = %s", getenv("OUTPUT_Z_DIR"));
    LOG_INFO("OUTPUT_RAIN_DIR = %s", getenv("OUTPUT_RAIN_DIR"));
    LOG_INFO("OUTPUT_Z_LOWRIS_DIR = %s", getenv("OUTPUT_Z_LOWRIS_DIR"));
    LOG_INFO("-----------------------------------------------------------------");
}

/// Report a command line error and quit
void commandline_error(const char* progname, const char* msg)
{
    fprintf(stderr, "%s\n\n", msg);
    fprintf(stderr, "Usage: %s volume_file file_type site_name\n", progname);
    exit(1);
}

/* ================================ */
int main (int argc, char **argv)
    /* ================================ */
{
    char *nome_file;
    int ier_main=0;//uscite errore generico (lettura volume e succ.anap) , di test_file, del main
    char *sito;//GAT O SPC
    int file_type; // -- x definire n_elev e reolution e se è =1 esco

    CUM_BAC_CLOPT CL_opt;
    parseOptions(argc,argv,&CL_opt);

    // Initialize logging
    Logging logging;

    LOG_CATEGORY("radar.main");

    int MyMAX_BIN = 512;
    if(CL_opt.do_medium) MyMAX_BIN=1024; 
    //------- verifica n0 argomenti ------

    if (argc < 4)
        commandline_error(argv[0], "some argument is missing.");

    nome_file 		= (char *)CL_opt.filename.c_str(); 
    file_type		=	  CL_opt.filetype;
    sito		= (char *)CL_opt.sito.c_str();

    if(CL_opt.do_medium && CL_opt.filetype == 3) MyMAX_BIN = 512;   // questo dovrebbe essere il caso del medio vecchio 

    elaboradar::Config cfg;

    setwork(sito);  //-------setto ambiente lavoro (se var amb lavoro non settate le setta in automatico) ------

    startup_banner(&CL_opt);


    unique_ptr<elaboradar::CUM_BAC> cb(new elaboradar::CUM_BAC(cfg, sito, CL_opt.do_medium,MyMAX_BIN));


    // Set feature flags
    cb->do_clean 	= CL_opt.do_clean;
    cb->do_quality 	= CL_opt.do_quality;
    cb->do_beamblocking = CL_opt.do_beamblocking;
    cb->do_declutter 	= CL_opt.do_declut;
    cb->do_bloccorr 	= CL_opt.do_bloccor;
    cb->do_vpr 		= CL_opt.do_vpr;
    cb->do_class 	= CL_opt.do_class;
    cb->do_devel 	= CL_opt.do_devel;
    cb->do_readStaticMap= CL_opt.do_readStaticMap;
    cb->do_zlr_media	= true;
    cb->do_anaprop	= CL_opt.do_anaprop;
    try {
   //    cb->StampoFlag(); 
       if (cb->esegui_tutto(nome_file, file_type,  CL_opt.data_in_odim)){
            ier_main = 0;
//	    unsigned irange=60000/cb->volume.scan(0).cell_size;
//            std::cout<<"cell size "<<cb->volume.scan(0).cell_size<<"\t Beam_count"<<cb->volume.scan(0).beam_count<<std::endl;
//            for (unsigned beam=130; beam <=170; beam ++) 
//               LOG_WARN(" Raggio %3d bin %4d -- db@60Km %3f",beam,irange,cb->volume.scan(0)(beam,irange)); 
//           
        }
        else
            ier_main = 1;
    } catch (std::exception& e) {
        LOG_ERROR("Errore nella processazione: %s", e.what());
        ier_main = 1;
    }

    // è stato tolto il loop sui volumi

    LOG_INFO("End of processing, result: %d", ier_main);

    return ier_main;
}

