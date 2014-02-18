
/*----------------------------------------------------------------------------*/
  /*    INCLUDE file                                      */
/*----------------------------------------------------------------------------*/

#include "cum_bac.h"
#include "logging.h"

// libreria c
#include <stdlib.h>
#include <math.h>

#include <setwork.h>

/*----------------------------------------------------------------------------*/
/*      FINE SEZIONE    INCLUDE                   */
/*----------------------------------------------------------------------------*/
// DEFINIZIONI PREPROCESSORE NB: VPR E CLASS IMPLICANO QUALITY ,

#if 0
// TODO: questi non sembrano piú servire: li cancelliamo?

//Definizioni geometriche
#define AMPLITUDE 0.9 /* esternalizzo?*/ // ampiezza fascio radar
#define DTOR  M_PI/180. /* esternalizzo?*/ //fattore conversione gradi-radianti
#define CONV_RAD 360./4096.*DTOR  // fattore conversione unità angolare radar-radianti

//Dimensioni massime range polare
#define MAX_BIN 512
#define MAX_DIM 512

// Fattori moltiplicativi per passare da unità azimut / unità elevazione a gradi
#define FATT_MOLT_EL ((double) 360./(double)4096.)
#define FATT_MOLT_AZ ((double) 360./(double)4096.)

// Soglie algoritmi
#define MISSING 0 /*valore mancante*/

// parametri ereditati da programma beam blocking:numero elevazioni da programma beam blocking ; le matrici ivi definite considerano questo
#define NSCAN 6
#endif

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
static void startup_banner()
{
    LOG_CATEGORY("radar.banner");

    /// Log initial program state
    LOG_INFO("Lancio Programma");
    LOG_INFO("-----------------------------------------------------------------");
    LOG_INFO("Flag di Compilazione: "
#ifdef BOLOGNA
            " BOLOGNA"
#endif
#ifdef SPC
            " SPC"
#endif
#ifdef WRITE_DBP
            " WRITE_DBP"
#endif
#ifdef TIME
            " TIME"
#endif
#ifdef WRITE_DBP_REORDER
            " WRITE_DBP_REORDER"
#endif
#ifdef DECLUTTER
            " DECLUTTER"
#endif
#ifdef Z_AVERAGE
            " Z_AVERAGE"
#endif
#ifdef R_AVERAGE
            " R_AVERAGE"
#endif
#ifdef Z_LOWRIS
            " Z_LOWRIS"
#endif
#ifdef SHORT
            " SHORT"
#endif
#ifdef MEDIUM
            " MEDIUM"
#endif
#ifdef STATIC
            " STATIC"
#endif
#ifdef BEAMBLOCKING
            " BEAMBLOCKING"
#endif
#ifdef QUALITY
            " QUALITY"
#endif
            ".");

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

    // Initialize logging
    Logging logging;

    LOG_CATEGORY("radar.main");

    //------- verifica n0 argomenti ------

    if (argc < 4)
        commandline_error(argv[0], "some argument is missing.");

    nome_file = argv[1];

    int file_type; // -- x definire n_elev e reolution e se è =1 esco
    if (sscanf(argv[2], "%d", &file_type) != 1)
        commandline_error(argv[0], "file_type must be a number");

    sito = argv[3];  //---- assegnazioni legate a argv[3]-- sito,  ambiente di lavoro, elevazioni
    setwork(argv[3]);  //-------setto ambiente lavoro (se var amb lavoro non settate le setta in automatico) ------

    startup_banner();

    cumbac::CUM_BAC *cb = new cumbac::CUM_BAC(sito);

    // Set feature flags
    cb->do_clean = true;
#ifdef SHORT
    cb->do_medium = false;
#else
 #ifdef MEDIUM
    cb->do_medium = true;
 #else
    #error Neither SHORT nor MEDIUM are defined
 #endif
#endif
#ifdef QUALITY
    cb->do_quality = true;
#endif
#ifdef BEAMBLOCKING
    cb->do_beamblocking = true;
#endif
#ifdef DECLUTTER
    cb->do_declutter = true;
#endif
#ifndef BLOCNOCORR
    cb->do_bloccorr = true;
#endif
#ifdef VPR
    cb->do_vpr = true;
#endif
#ifdef CLASS
    cb->do_class = true;
#endif
#ifdef STAMPE_EXTRA
    cb->do_devel = true;
#endif

    try {
        if (cb->esegui_tutto(nome_file, file_type))
            ier_main = 0;
        else
            ier_main = 1;
    } catch (std::exception& e) {
        LOG_ERROR("Errore nella processazione: %s", e.what());
        ier_main = 1;
    }
    delete cb;

    // è stato tolto il loop sui volumi

    LOG_INFO("End of processing, result: %d", ier_main);

    return ier_main;
}

