
/*----------------------------------------------------------------------------*/
  /*    INCLUDE file                                      */
/*----------------------------------------------------------------------------*/

#include "cum_bac.h"
#include "logging.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <func_Z_R.h>
#ifdef __cplusplus
}
#endif

#include <Q_components.h>


// libreria c
#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <string.h>
#include <memory.h>
#include <unistd.h>
#include <errno.h>



#ifdef QUALITY
#include <qual_par.h>
#endif
#ifdef CLASS
#include <par_class.h>
#endif

/*----------------------------------------------------------------------------*/
/*      FINE SEZIONE    INCLUDE                   */
/*----------------------------------------------------------------------------*/
// DEFINIZIONI PREPROCESSORE NB: VPR E CLASS IMPLICANO QUALITY ,

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

/*----------------------------------------------------------------------------*/
/*  DICHIARATIVE GLOBALI                              */
/*----------------------------------------------------------------------------*/

/*extern char *sys_errlist[];  ANNA 17-2-2006 sostituito con strerror che sta in string.h */
extern int errno;

namespace cumbac {

void CUM_BAC::scrivo_out_file_bin (const char *ext,const char *content,const char *dir,size_t size, const void  *matrice)
{
    char nome_file [512];
    FILE *output;
    struct tm *tempo;
    time_t time;
    /*----------------------------------------------------------------------------*/
    /*    apertura file dati di output                                          */
    /*----------------------------------------------------------------------------*/
    //definisco stringa data in modo predefinito
#ifdef SHORT
    time = NormalizzoData(volume.acq_date);
    tempo = gmtime(&time);
#endif
#ifdef MEDIUM
    tempo = gmtime(&volume.acq_date);
    time = NormalizzoData(volume.acq_date);
    tempo = gmtime(&time);
#endif
    snprintf(nome_file, 512, "%s/%04d%02d%02d%02d%02d%s",dir,
            tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
            tempo->tm_hour, tempo->tm_min, ext);

    output=controllo_apertura(nome_file,content,"w");
    printf("aperto file %s dimensione matrice %zd\n",nome_file,size);

    fwrite(matrice,size,1,output);
    fclose(output);
    return;
}

/*=======================================================================================*/

void CalcoloVPR::esegui_tutto()
{
    test_vpr=fopen(getenv("TEST_VPR"),"a+");

    LOG_INFO("processo file dati: %s", cum_bac.volume.filename.c_str());
    printf ("calcolo VPR \n") ;

    //VPR  // ------------inizializzo hvprmax ---------------

    hvprmax=INODATA;

    //VPR  // ------------chiamo combina profili con parametri sito, sito alternativo ---------------

    //  ier_comb=combina_profili(sito,argv[4]);
    ier_comb=combina_profili();
    printf ("exit status calcolo VPR istantaneo: (1--fallito 0--ok)  %i \n",ier_vpr) ; // debug
    printf ("exit status combinaprofili: (1--fallito 0--ok) %i \n",ier_comb) ; // debug


    //VPR  // ------------chiamo profile_heating che calcola riscaldamento profilo ---------------

    heating=profile_heating();
    printf ("heating %i \n", heating);
    LOG_INFO("ier_vpr %i ier_comb %i",ier_vpr,ier_comb);

    //VPR  // ------------se combina profili ok e profilo caldo correggo --------------
    if (!ier_comb && heating >= WARM){

        int ier=corr_vpr();
        printf ("exit status correggo vpr: (1--fallito 0--ok) %i \n",ier) ; // debug


        //VPR // ------------se la correzione è andata bene e il profilo è 'fresco' stampo profilo con data-------

        if ( ! ier && ! ier_vpr)
            ier_stampa_vpr=stampa_vpr();
    }
}

bool CUM_BAC::esegui_tutto(const char* nome_file, int file_type)
{
    // Legge e controlla il volume dal file SP20
    if (!read_sp20_volume(nome_file, file_type))
        return false;

    ///-------------------------ELABORAZIONE -------------------------

    //  ----- da buttare : eventuale scrittura del volume polare ad azimut e range fissi
    /* #ifdef WRITE_DBP_REORDER  */

    /*       /\*------------------------------------------------------------------  */
    /*    | eventuale scrittura del volume polare ad azimut e range fissi  |  */
    /*    -----------------------------------------------------------------*\/   */
    /*       strcat(nome_file,"_reorder");        */
    /*       ScrivoLog(6,nome_file);  */
    /*       ier=write_dbp(nome_file);  */
    /* #endif */

    //  ----- test su normalizzazione data ( no minuti strani)
    if (NormalizzoData(volume.acq_date) == -1)
        return true;

    setup_elaborazione(nome_file);

    //--------------se def anaprop : rimozione propagazione anomala e correzione beam blocking-----------------//
#ifdef ANAPROP
    LOG_INFO("inizio rimozione anaprop e beam blocking");
    int ier = elabora_dato();
#endif


    //--------------se definita la qualita procedo con il calcolo qualita e del VPR (perchè prendo solo i punti con qual > soglia?)-----------------//

    if (do_quality)
    {
        //-------------calcolo qualita' e trovo il top
        printf ("calcolo Q3D \n") ;
        caratterizzo_volume();

        /* //---------trovo il top (a X dbZ) */
        /* printf ("trovo top \n") ; */
        /* ier=trovo_top(); */


        //--------------se definita CLASS procedo con  classificazione -----------------//
        if (do_class)
            calcolo_vpr->classifica_rain();

        //--------------se definito VPR procedo con calcolo VPR -----------------//

        if (do_vpr)
            calcolo_vpr->esegui_tutto();
    }

    if (do_class)
        class_conv_fixme_find_a_name();

#ifdef TIME
    prendo_tempo();
#endif


    //--------------------da rimuovere eventuale scrittura volume ripulito

    /* #ifdef WRITE_DBP  */
    /*      /\*-------------------------------------------------  */
    /*        | eventuale scrittura del volume polare ripulito e |  */
    /*        -------------------------------------------------*\/   */
    /* #ifdef DECLUTTER  */
    /*      strcat(nome_file,"_decl");  */
    /* #else  */
    /*      strcat(nome_file,"_anap");  */
    /* #endif  */
    /*      /\*exit (1);*\/  */
    /*      ScrivoLog(8,nome_file);  */
    /*      ier=write_dbp(nome_file);  */
    /* #endif  */




    //------------------- conversione di coordinate da polare a cartesiana se ndef  WRITE_DBP -----------------------
#ifndef WRITE_DBP

    /*--------------------------------------------------
      | conversione di coordinate da polare a cartesiana |
      --------------------------------------------------*/
    LOG_INFO("%s -- Creazione Matrice Cartesiana",PrendiOra());
    creo_cart();

    //------------------- STA PRENDOTEMPO!!!!!!!!!!!!!!!
#ifdef TIME
    prendo_tempo();
#endif



    //-------------------Se definita Z_LOWRIS creo matrice 1X1  ZLR  stampo e stampo coeff MP (serve?)------------------

#ifdef Z_LOWRIS

    LOG_INFO("%s -- Estrazione Precipitazione 1X1", PrendiOra());
    creo_cart_z_lowris();

    //-------------------scritture output -----------------------
    LOG_INFO("%s -- Scrittura File Precipitazione 1X1 %s\n", PrendiOra(), nome_file);
    scrivo_out_file_bin(".ZLR","dati output 1X1",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(z_out.data),z_out.data);


    char nome_file_output[512];
    sprintf(nome_file_output,"%s/MP_coeff",getenv("OUTPUT_Z_LOWRIS_DIR"));
    FILE* output=controllo_apertura(nome_file_output,"file coeff MP","w");
    fwrite(MP_coeff,sizeof(MP_coeff),1,output);
    fclose(output);
    printf(" dopo scrivo_z_lowris\n");



#ifdef TIME
    prendo_tempo();
#endif





    //-------------------scritture output -----------------------
    if (do_quality)
    {
        //------------------ output qual  per operativo:qualità in archivio e elevazioni e anap in dir a stoccaggio a scadenza
        // temporanee
        // in archivio
        scrivo_out_file_bin(".qual_ZLR","file qualita' Z",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(qual_Z_1x1.data),qual_Z_1x1.data);

        //------------------ stampe extra per studio
#ifdef STAMPE_EXTRA
        //scrivo_out_file_bin(".corrpt","file anap",getenv("DIR_QUALITY"),sizeof(dato_corrotto),dato_corrotto);
        //scrivo_out_file_bin(".pia","file PIA",getenv("DIR_QUALITY"),sizeof(att_cart),att_cart);
        //scrivo_out_file_bin(".bloc","file bloc",getenv("DIR_QUALITY"),sizeof(beam_blocking),beam_blocking);
        //scrivo_out_file_bin(".quota","file quota",getenv("DIR_QUALITY"),sizeof(quota),quota);
        scrivo_out_file_bin(".elev","file elevazioni",getenv("DIR_QUALITY"),sizeof(elev_fin),elev_fin);
        scrivo_out_file_bin(".bloc_ZLR","file bloc",getenv("DIR_QUALITY"),sizeof(beam_blocking_1x1.data),beam_blocking_1x1.data);
        scrivo_out_file_bin(".anap_ZLR","file anap",getenv("DIR_QUALITY"),sizeof(dato_corr_1x1.data),dato_corr_1x1.data); //flag di propagazione anomala
        scrivo_out_file_bin(".quota_ZLR","file qel1uota",getenv("DIR_QUALITY"),sizeof(quota_1x1.data),quota_1x1.data);// m/100 +128
        scrivo_out_file_bin(".elev_ZLR","file elev",getenv("DIR_QUALITY"),sizeof(elev_fin_1x1.data),elev_fin_1x1.data);
        scrivo_out_file_bin(".top20_ZLR","file top20",getenv("DIR_QUALITY"),sizeof(top_1x1.data),top_1x1.data);

#endif
        //------------------ stampe correzioni da profili verticali in formato ZLR
        if (do_vpr)
        {
            scrivo_out_file_bin(".corr_ZLR","file correzione VPR",getenv("DIR_QUALITY"),sizeof(corr_1x1.data),corr_1x1.data);
            //scrivo_out_file_bin(".neve","punti di neve",getenv("DIR_QUALITY"),sizeof(neve),neve);
            // scrivo_out_file_bin(".neve_ZLR","file presunta neve ",getenv("DIR_QUALITY"),sizeof(neve_1x1),neve_1x1);
        }

        //------------------se definita CLASS  stampo punti convettivi
        if (do_class)
        {
            // in archivio?
            // scrivo_out_file_bin(".conv_ZLR","punti convettivi",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(conv_1x1),conv_1x1);
            scrivo_out_file_bin(".conv_ZLR","punti convettivi",getenv("DIR_QUALITY"),sizeof(conv_1x1.data),conv_1x1.data);
        }
    }


#endif// ifdef Z_lOWRIS
#endif // ifndef WRITE_DBP

    return true;
}

}

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
    LOG_INFO("%s -- Lancio Programma", PrendiOra());
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
#ifdef ANAPROP
            " ANAPROP"
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

#ifdef TIME
    prendo_tempo();
#endif

    cumbac::CUM_BAC *cb = new cumbac::CUM_BAC(sito);

    // Set feature flags
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

