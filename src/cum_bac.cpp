#include "cum_bac.h"
#include "logging.h"

#include <cstring>

#ifdef __cplusplus
extern "C" {
#endif
// libreria radar
//#include <radar_parameter.h>
#ifdef __cplusplus
}
#endif

//Definizioni per test_file
#define NUM_MIN_BEAM 200
#define SHORT_DEC         0
#define SHORT_FULL_VOLUME 1
#define SHORT_HAIL        2
#define MEDIUM_PULSE      3

/// This needs to be a global variable, as it is expected by libsp20
int elev_array[NEL];


CUM_BAC::CUM_BAC()
{
    logging_category = log4c_category_get("radar.cum_bac");

#ifdef VPR
    t_ground=NODATAVPR;
    chisqfin=100; //???puo' essere def in anal
    rmsefin=100;
#endif
#ifdef CLASS
    ncv=0;ncs=0;np=0;
    htbb=-9999.; hbbb=-9999.;
#endif
    memset(vol_pol,0,sizeof(vol_pol));
    memset(cart,0,sizeof(cart));
  #ifdef ZLR_MEDIA
    memset(cartm,0.,sizeof(cartm));
  #endif
    memset(z_out,0,sizeof(z_out));
    memset(nbeam_elev,0,sizeof(nbeam_elev));
    memset (first_level,0,sizeof(first_level));
    memset (first_level_static,0,sizeof(first_level_static));
  #ifdef QUALITY
    memset(dato_corrotto,0,sizeof(dato_corrotto));
    memset(dato_corr_xy,0,sizeof(dato_corr_xy));
    memset(dato_corr_1x1,0,sizeof(dato_corr_1x1));
    memset(beam_blocking,0,sizeof(beam_blocking));
    memset(beam_blocking_xy,0,sizeof(beam_blocking_xy));
    memset(beam_blocking_1x1,0,sizeof(beam_blocking_1x1));
    memset(att_cart,DBtoBYTE(0.),sizeof(att_cart));
    memset(quota_rel,0,sizeof(quota_rel));
    memset(quota,0,sizeof(quota));
    memset(quota_cart,0,sizeof(quota_cart));
    memset(quota_1x1,0,sizeof(quota_1x1));
    memset(elev_fin,0,sizeof(elev_fin));
    memset(elev_fin_xy,0,sizeof(elev_fin_xy));
    memset(elev_fin_1x1,0,sizeof(elev_fin_1x1));
    memset(qual,0,sizeof(qual));
    memset(qual_Z_cart,0,sizeof(qual_Z_cart));
    memset(qual_Z_1x1,0,sizeof(qual_Z_1x1));
    memset(top,0,sizeof(top));
    memset(topxy,0,sizeof(topxy));
    memset(top_1x1,0,sizeof(top_1x1));
  #ifdef VPR
    memset(corr_1x1,0,sizeof(corr_1x1));
    memset(neve_cart,0,sizeof(neve_cart));
    memset(neve_1x1,0,sizeof(neve_1x1));
    memset(neve,0,sizeof(neve));
    for (int i=0; i<NMAXLAYER; i++)
      vpr[i]=NODATAVPR;
    for (int l=0; l<NEL; l++)
      for (int i=0; i<NUM_AZ_X_PPI; i++)
        for(int k=0; k<MAX_BIN; k++)
      flag_vpr[l][i][k]=0;

  #ifdef CLASS
    for (int k=0; k<NUM_AZ_X_PPI*MAX_BIN;k++ ){
      lista_conv[k][0]=-999;
      lista_conv[k][1]=-999;
      lista_bckg[k][0]=-999;
      lista_bckg[k][1]=-999;
    }

    memset(stratiform,0,sizeof(stratiform));
  #endif

  #endif
  #endif
  #ifdef BEAMBLOCKING
    memset(bb_first_level,0,sizeof(bb_first_level));
  #endif

    memset(stat_anap_tot,0,sizeof(stat_anap_tot));
    memset(stat_anap,0,sizeof(stat_anap));
    memset(stat_bloc,0,sizeof(stat_bloc));
    memset(stat_elev,0,sizeof(stat_elev));


    // ----- alloco i puntatori ---------
    nome_fl=(char *)malloc(200*sizeof(char));
    nome_dem=(char *)malloc(200*sizeof(char));

    //-----  FINE INIZIALIZZAZIONI---------//
}

bool CUM_BAC::test_file(int file_type)
{
    FILE *f_aus;
    time_t last_time; //> old_data_header.norm.maq.acq_date?
    int n_elev, resolution;// != old_data_header.norm.maq.resolution?

    //--- switch tra tipo di file per definire nelev = elevazioni da testare e la risoluzione

    switch (file_type)
    {
        case SHORT_DEC:
            if (!old_data_header.norm.maq.declutter_rsp)
            {
                LOG_WARN("File Senza Declutter Dinamico--cos' è???");
                return false;
            }
            resolution=2;
            n_elev=4;
            break;
            //------------se tipo =1 esco
        case SHORT_FULL_VOLUME://-----??? DUBBIO
            if (old_data_header.norm.maq.declutter_rsp)
            {
                LOG_WARN("File con Declutter Dinamico");
                return false;
            }
            resolution=2;
            n_elev=4;
            break;
        case SHORT_HAIL://-----??? DA BUTTARE NON ESISTE PIÙ
            resolution=2;
            n_elev=3;
            LOG_INFO("CASO SHORT_HAIL");
            break;
        case MEDIUM_PULSE:
            resolution=4;
            n_elev=4;
            break;
    }

    //----------se la risoluzione del file è diversa da quella prevista dal tipo_file dà errore ed esce (perchè poi probabilmente le matrici sballano ?)
    if (old_data_header.norm.maq.resolution != resolution)
    {
        LOG_ERROR("File Risoluzione Sbagliata %1d", old_data_header.norm.maq.resolution);
        return false;
    }
    //------eseguo test su n0 beam  sulle prime 4 elevazioni, se fallisce  esco ------------

    for (int k = 0; k < n_elev; k++) /* testo solo le prime 4 elevazioni */
    {
        LOG_INFO("Numero beam presenti : %4d  -- elevazione%2d\n", nbeam_elev[k], k);

        if (nbeam_elev[k] <  NUM_MIN_BEAM)
            // se numero beam < numero minimo---Scrivolog ed esco !!!!!!!!!!!!!!!!!!!
        {
            //---Scrivolog!!!!!!!!!!!!!!!!!!!
            LOG_ERROR("Trovati Pochi Beam Elevazione %2d - num.: %3d",k,nbeam_elev[k]);
            return false;
        }
    }                                                             /*end for*/

    //--------verifico la presenza del file contenente l'ultima data processata-------

    const char* last_file = getenv("LAST_FILE");
    if (last_file != NULL)
    {
        if (access(last_file, 6) == 0)
        {
            /*--------------------------------------
              |  il file e' presente, leggo la data  |
              |  e la confronto con la data del      |
              |  volume dati in esame, se il dbp e'  |
              |  piu' giovane continuo altrimenti    |
              |  do' un avviso (una volta errore e uscivo) perchè processavo un file alla volta|
              |  adesso la logica è cambiata                           |
              -------------------------------------*/
            f_aus = fopen(last_file, "r+");
            fread(&last_time,4,1,f_aus);
            if(old_data_header.norm.maq.acq_date <= last_time)
            {
                fclose(f_aus);
                LOG_WARN("File Vecchio");
                //return false;
            } else {
                /*----------------------------
                  |  aggiorno la data nel file |
                  ----------------------------*/
                rewind(f_aus);
                fwrite(&old_data_header.norm.maq.acq_date,4,1,f_aus);
                fclose(f_aus);
            }
        }
        //-------- fin qui tutto un pezzo per dire che se il file è più vecchio dell'ultimo do' errore -----------

        // -------- se il file ultima data non è presente lo scrivo
        else
        {
            /*--------------------------------------
              |  il file non e' presente, scrivo la  |
              |  data del volume in esame            |
              --------------------------------------*/
            f_aus = fopen(last_file, "w");
            fwrite(&old_data_header.norm.maq.acq_date,4,1,f_aus);
            fclose(f_aus);
        }
    }
    // ------- se ok status di uscita:1
    return true;
}

bool CUM_BAC::read_sp20_volume(const char* nome_file, const char* sito, int file_type)
{
    // ----- definisco array delle elevazioni che è diverso per i due siti ---------
    if (!(strcmp(sito,"SPC")) ) {
        for (int i = 0; i < NEL; ++i)
            elev_array[i] = elev_array_spc[i];
    }

    if (!(strcmp(sito,"GAT")) ) {
        for (int i = 0; i < NEL; ++i)
            elev_array[i] = elev_array_spc[i];
    }

    LOG_INFO("Reading %s for site %s and file type %d", nome_file, sito, file_type);

    //--------lettura volume------
    int tipo_dati_richiesti = INDEX_Z;
    int ier = read_dbp_SP20((char*)nome_file,vol_pol,&old_data_header,
                            tipo_dati_richiesti,nbeam_elev);

    if (ier != OK)
        LOG_ERROR("Reading %s returned error code %d", nome_file, ier);

#ifdef TIME
    prendo_tempo();
#endif

    //  ----- Test sul volume test_file.......  --------
    if (!test_file(file_type))
    {
        LOG_ERROR("test_file failed");
        return false;
    }

    return ier == OK;
}

/*
 * Questa funzione al momento non va tolta. C'è del codice che va a leggere
 * fuori da dove dovrebbe, e se togliamo questa variabile statica, il layout
 * del programma in ram cambia sufficientemente per farlo segfaultare.
 *
 * Mettiamo in programma di fare qualche giro con valgrind per vedere dove sono
 * le letture sbagliate, e una volta messe a posto finalmente possiamo togliere
 * la ScrivoLog, che non è piú chiamata da nessuno.
 */
void ScrivoLog(int i, const char* stringa)
{
    static FILE *log = 0;
}

char *PrendiOra()
{
    time_t clock;
    struct tm *tempo;

    clock=time(0);

    tempo=gmtime(&clock);

    return   asctime(tempo);
}

void prendo_tempo()
{
    static time_t time_tot = 0,time1 = 0,time2 = 0;
    LOG_CATEGORY("radar.timing");

    if(time1 == 0){
        time1=time(&time1);
        time_tot = time1;
    }
    time2 = time(&time2);

    LOG_INFO("tempo parziale %ld ---- totale %ld\n", time2-time1, time2-time_tot);
    time1=time2;
    return;
}

/* time è in secondi, itime è un intero che rappresenta il numero intero di intervalli da 5 minuti*/
time_t NormalizzoData(time_t time)
{
    int itime;

    itime = time/(NMIN*60);

    /*
       printf(" esco da Normalizzo %d %d %d \n",time,itime,(time - itime*NMIN*60));
       printf("%s\n",ctime(&time));
       printf("%d\n",(NMIN-MAX_TIME_DIFF)*60);
       */

    if(time - itime*NMIN*60 <MAX_TIME_DIFF*60) return (itime*NMIN*60); /* se la differenza è meno di tre minuti vado al 5° min. prec*/
    if(time - itime*NMIN*60 >(NMIN-MAX_TIME_DIFF)*60) return ((itime+1)*NMIN*60); /* se la differenza è più di tre minuti vado al 5° min. successivo*/
    //altrimenti ritorno -1
    return -1;
}


