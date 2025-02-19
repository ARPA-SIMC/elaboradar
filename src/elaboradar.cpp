/**
 *  @file
 *  @defgroup elaboradar 
 *  @brief progetto elaborazione dati radar per ottenere campi di Z da convertire in R
 *  @details elaborazione dati radar utilizzando dati da radiosondaggio, temperatura e dati in uscita da programma beam blocking che esegue un controllo di qualita' del volume, rimuove la propagazione anomala, corregge il beam blocking,   calcola il profilo verticale, calcola la  qualità, classifica le aree convettive  
*/

/*----------------------------------------------------------------------------*/
/*    INCLUDE file                                      */
/*----------------------------------------------------------------------------*/

#include "cum_bac.h"
#include <radarelab/logging.h>
#include "config.h"
#include "site.h"
#include "cartproducts.h"
#include "cum_bac_clparser.h"

#include <memory>
#include <cstdlib>
#include <cmath>

#include <setwork.h>

using namespace std;
using namespace radarelab;
using namespace elaboradar;

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
//            " WRITE_DBP_REORDER" eventualemnte da re-inserire
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
 
    LOG_INFO("%s", FlagRunTime.c_str());

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

void check_volume(const Volume<double>& volume, int file_type)
{
//Definizioni per test_file
#define NUM_MIN_BEAM 200
#define SHORT_DEC         0
#define SHORT_FULL_VOLUME 1
#define SHORT_HAIL        2
#define MEDIUM_PULSE      3
#define SHORT_212	  4
    LOG_CATEGORY("radar.io");

    unsigned n_elev = 0;
    int expected_size_cell = 0;// != volume.resolution?

    //--- switch tra tipo di file per definire nelev = elevazioni da testare e la risoluzione
    switch (file_type)
    {
        case SHORT_DEC:
            LOG_INFO("CASO SHORT_DEC");
//            if (!volume.load_info->declutter_rsp)
//                throw runtime_error("File Senza Declutter Dinamico--cos' è???");
            expected_size_cell = 250;
            n_elev=4;
            break;
            //------------se tipo =1 esco
        case SHORT_FULL_VOLUME://-----??? DUBBIO
            LOG_INFO("CASO SHORT_FULL_VOLUME");
            if (volume.load_info->declutter_rsp)
                throw runtime_error("File con Declutter Dinamico");
            expected_size_cell = 250;
            n_elev=4;
            break;
        case SHORT_HAIL://-----??? DA BUTTARE NON ESISTE PIÙ
            LOG_INFO("CASO SHORT_HAIL");
            expected_size_cell = 250;
            n_elev=3;
            break;
        case MEDIUM_PULSE:
            LOG_INFO("CASO MEDIO OLD");
            expected_size_cell = 1000;
            n_elev=4;
            break;
        case SHORT_212://----- CORRISPONDE A VOL_NEW - da questo si ottengono il corto e il medio
            LOG_INFO("CASO SHORT_212");
//            if (!volume.load_info->declutter_rsp)
//                throw runtime_error("File senza Declutter Dinamico");
            expected_size_cell = 250;
            n_elev=4;
            break;
    }

    //----------verifico che la risoluzione di ogni elevazione sia la stessa
    if ( ! volume.is_unique_cell_size())
    {
        LOG_ERROR("File Risoluzione/size_cell non costante nel volume");
        throw runtime_error("File Risoluzione/size_cell non costante");
    }
    //----------se la risoluzione del file è diversa da quella prevista dal tipo_file dà errore ed esce (perchè poi probabilmente le matrici sballano ?)
    if (volume[0].cell_size != expected_size_cell)
    {
        LOG_ERROR("File Risoluzione/size_cell Sbagliata %f", volume[0].cell_size);
        throw runtime_error("File Risoluzione/size_cell Sbagliata");
    }
    //------eseguo test su n0 beam  sulle prime 4 elevazioni, se fallisce  esco ------------

    if (volume.size() < n_elev)
    {
        LOG_ERROR("Volume has %zd elevations, but we are expecting at least %d", volume.size(), n_elev);
        throw runtime_error("Insufficient elevation count");
    }

unsigned count_good = 0, bin_in_vol=0;
    for (unsigned k = 0; k < volume.size(); k++) /* testo solo le prime 4 elevazioni */
    {
        LOG_INFO("Numero beam presenti: %4u -- elevazione %d", volume[k].beam_count, k);

        if (volume[k].beam_count < NUM_MIN_BEAM)
            // se numero beam < numero minimo---Scrivolog ed esco !!!!!!!!!!!!!!!!!!!
        {
            LOG_ERROR("Trovati Pochi Beam Elevazione %2d - num.: %3d", k, volume[k].beam_count);
            throw runtime_error("Insufficient beam count");
        }

	count_good += (volume[k].array() > -19.0).count();
	bin_in_vol += volume.beam_count*volume[k].beam_size;

    }

// Check if volume data are bad (too many data probably fault due to AFC control
    if  ((float) (count_good)/bin_in_vol > 0.9 ){
      LOG_ERROR("Trovati troppi dati > -19 dBZ  %5.2f\%, volume con possibile problema di AFC", (float) (count_good)/bin_in_vol*100.);
      throw runtime_error ("Possible bad volume (AFC)");
    }
}

/*
void CUM_BAC::StampoFlag(){
    std::cout<<" Flag do_medium       :"<< (this->do_medium?" true":" false")<<std::endl;
    std::cout<<" Flag do_clean        :"<< (this->do_clean?" true":" false")<<std::endl;  
    std::cout<<" Flag do_quality      :"<< (this->do_quality?" true":" false")<<std::endl;
    std::cout<<" Flag do_beamblocking :"<< (this->do_beamblocking?" true":" false")<<std::endl;
    std::cout<<" Flag do_declutter    :"<< (this->do_declutter?" true":" false")<<std::endl;
    std::cout<<" Flag do_bloccorr     :"<< (this->do_bloccorr?" true":" false")<<std::endl;
    std::cout<<" Flag do_vpr          :"<< (this->do_vpr?" true":" false")<<std::endl;
    std::cout<<" Flag do_class        :"<< (this->do_class?" true":" false")<<std::endl;
    std::cout<<" Flag do_zlr_media    :"<< (this->do_zlr_media?" true":" false")<<std::endl;
    std::cout<<" Flag do_devel        :"<< (this->do_devel?" true":" false")<<std::endl;
    std::cout<<" Flag do_readStaticMap:"<< (this->do_readStaticMap?" true":" false")<<std::endl;  
}
*/

/* ================================ */
int main (int argc, char **argv)
    /* ================================ */
{
    char *nome_file;
    int ier_main=0;//uscite errore generico (lettura volume e succ.anap) , di test_file, del main
    char *sito;//GAT O SPC
    char *fuzzypath;
    int file_type; // -- x definire n_elev e reolution e se è =1 esco

    CUM_BAC_CLOPT CL_opt;
    parseOptions(argc,argv,&CL_opt);
    PrintOptions(&CL_opt);


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
    fuzzypath          = (char *)CL_opt.fuzzy_path.c_str();

    if(CL_opt.do_medium && CL_opt.filetype == 3) MyMAX_BIN = 512;   // questo dovrebbe essere il caso del medio vecchio 

    elaboradar::Config cfg;
    
    setwork(sito);  //-------setto ambiente lavoro (se var amb lavoro non settate le setta in automatico) ------

    startup_banner(&CL_opt);

    const Site& site(Site::get(sito));
    site.datipath = fuzzypath;
    Volume<double> volume(NUM_AZ_X_PPI);

    cout<<"volume size = "<<volume.size()<<endl;
    try {
      if (CL_opt.data_in_odim)
            // Legge e controlla il volume dal file ODIM
	CUM_BAC::read_odim_volume(volume, site, nome_file, fuzzypath, CL_opt.do_clean, CL_opt.do_medium, CL_opt.set_undetect);
        else
            // Legge e controlla il volume dal file SP20
            CUM_BAC::read_sp20_volume(volume, site, nome_file, CL_opt.do_clean, CL_opt.do_medium);
    } catch (std::exception& e) {
        LOG_ERROR("Errore nel caricamento del volume: %s", e.what());
        return 2;
    }
    check_volume(volume, file_type);    

    unique_ptr<elaboradar::CUM_BAC> cb(new elaboradar::CUM_BAC(volume, cfg, site, CL_opt.do_medium,MyMAX_BIN));
    // Set feature flags
    cb->do_quality 	= CL_opt.do_quality;
    cb->do_beamblocking = CL_opt.do_beamblocking;
    cb->do_declutter 	= CL_opt.do_declut;
    cb->do_bloccorr 	= CL_opt.do_bloccor;
    cb->do_class 	= CL_opt.do_class;
    cb->do_devel 	= CL_opt.do_devel;
    cb->do_readStaticMap= CL_opt.do_readStaticMap;
    cb->do_zlr_media	= true;
    cb->do_anaprop	= CL_opt.do_anaprop;

    printwork();
    std::string algos;
    if (CL_opt.do_clean) 	algos = algos+"CL";
    if (CL_opt.do_anaprop)	algos = algos+"_AP"; 
    if (CL_opt.do_beamblocking)	algos = algos+"_BB"; 
    try {
        //--------------se def anaprop : rimozione propagazione anomala e correzione beam blocking-----------------//
        LOG_INFO("inizio rimozione anaprop e beam blocking");
        for(unsigned k=0; k<volume.size(); k++) LOG_INFO(" SCAN # %2d - BeamSIZE %4d",k,volume[k].beam_size);
        cb->declutter_anaprop();
        cb->caratterizzo_volume();

        unsigned CART_DIM_ZLR  = CL_opt.do_medium ? 512 : 256;
        unsigned ZLR_N_ELEMENTARY_PIXEL = volume.at(0).cell_size == 1000. ? 1 : 4;
        if (CL_opt.do_SaveBothRanges) {
	  if (volume.max_beam_size() * volume.at(0).cell_size * 0.001 > 128 ){
	    CART_DIM_ZLR = 512;
          } else {
	    CART_DIM_ZLR = 256;
            CL_opt.do_SaveBothRanges = false;
          }
        }

        if (CL_opt.do_intermediateProd){
          CartProducts products(volume, CART_DIM_ZLR, ZLR_N_ELEMENTARY_PIXEL);
          cb->generate_maps(products);
          products.write_out(cb->assets,CART_DIM_ZLR, algos);
	  OdimProdDefs odimProd(products.z_out, products.qual_Z_1x1, products.ScaledRes);
          products.write_odim(cb->assets, CART_DIM_ZLR, algos, odimProd);
	  if (CL_opt.do_SaveFullRes){
            OdimProdDefs odimProdFullRes(products.z_fr, products.FullsizeRes);
            products.write_odim(cb->assets, 256. * ZLR_N_ELEMENTARY_PIXEL, algos, odimProdFullRes);
	  }
          if (CL_opt.do_SaveBothRanges){
            LOG_INFO("Salvo sub-image intermedie");
            products.write_out(cb->assets, 256, algos);
	    products.write_odim(cb->assets, 256, algos,odimProd);
	  }
        }

        // cb->StampoFlag();
        if (CL_opt.do_vpr) {
	  cb->want_vpr();
          if (cb->do_quality)    cb->vpr_class();
          if (CL_opt.do_vpr)		algos = algos+"_VPR"; 
          if (CL_opt.do_class)		algos = algos+"_CLASS"; 
       } 
       CartProducts products(volume, CART_DIM_ZLR, ZLR_N_ELEMENTARY_PIXEL);
       cb->generate_maps(products);
       OdimProdDefs odimProd(products.z_out, products.qual_Z_1x1, products.ScaledRes);
       if (CL_opt.do_SaveBothRanges){
            products.write_out(cb->assets,CART_DIM_ZLR,algos);
	    products.write_odim(cb->assets, CART_DIM_ZLR, algos, odimProd);
            LOG_INFO("Salvo sub-image");
            products.write_out(cb->assets, 256, algos);
	    products.write_odim(cb->assets, 256, algos, odimProd);
       }else{
            products.write_out(cb->assets,CART_DIM_ZLR,algos);
	    products.write_odim(cb->assets, CART_DIM_ZLR, algos, odimProd);
          // products.write_out(cb->assets);
       }	
       if (CL_opt.do_SaveFullRes){
         OdimProdDefs odimProdFullRes(products.z_fr, products.FullsizeRes);
         products.write_odim(cb->assets, 256. * ZLR_N_ELEMENTARY_PIXEL, algos, odimProdFullRes);
       }
    } catch (std::exception& e) {
        LOG_ERROR("Errore nella processazione: %s", e.what());
        ier_main = 1;
    }

    LOG_INFO("End of processing, result: %d", ier_main);
    return ier_main;
}

