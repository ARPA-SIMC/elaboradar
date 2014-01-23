#include "cum_bac.h"

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

int CUM_BAC::test_file(int file_type)
{
    FILE *f_aus;
    time_t last_time; //> old_data_header.norm.maq.acq_date?
    int n_elev, resolution;// != old_data_header.norm.maq.resolution?
    int k;

    // ---- attribuisce la stringa di input (che è un numero 0-3) a un intero: file_type-----
    printf("tipo file %1d\n",file_type);

    //--- switch tra tipo di file per definire nelev = elevazioni da testare e la risoluzione

    switch (file_type)
    {
        case SHORT_DEC:
            if(!old_data_header.norm.maq.declutter_rsp )
            {
                //---- eddaje sta scrivolog!!!!------
                sprintf(errori,"File Senza Declutter Dinamico--cos' è???");
                ScrivoLog(16,errori);
            }
            resolution=2;
            n_elev=4;
            break;
            //------------se tipo =1 esco
        case SHORT_FULL_VOLUME://-----??? DUBBIO
            if(old_data_header.norm.maq.declutter_rsp )
            {
                sprintf(errori,"File con Declutter Dinamico");
                ScrivoLog(16,errori);
                return 0;
            }
            resolution=2;
            n_elev=4;
            break;
        case SHORT_HAIL://-----??? DA BUTTARE NON ESISTE PIÙ
            resolution=2;
            n_elev=3;
            printf("CASO SHORT_HAIL\n");
            break;
        case MEDIUM_PULSE:
            resolution=4;
            n_elev=4;
            break;
    }

    //----------se la risoluzione del file è diversa da quella prevista dal tipo_file dà errore ed esce (perchè poi probabilmente le matrici sballano ?)
    if(old_data_header.norm.maq.resolution != resolution)
    {
        sprintf(errori,"File Risoluzione Sbagliata %1d", old_data_header.norm.maq.resolution);
        ScrivoLog(16,errori);
        return 0;
    }
    //------eseguo test su n0 beam  sulle prime 4 elevazioni, se fallisce  esco ------------

    for(k=0; k<n_elev; k++) /* testo solo le prime 4 elevazioni */
    {

        printf(" numero beam presenti : %4d  -- elevazione%2d\n",
                nbeam_elev[k],k);

        if(nbeam_elev[k] <  NUM_MIN_BEAM)
            // se numero beam < numero minimo---Scrivolog ed esco !!!!!!!!!!!!!!!!!!!
        {
            //---Scrivolog!!!!!!!!!!!!!!!!!!!
            sprintf(errori,"Trovati Pochi Beam Elevazione %2d - num.: %3d",k,nbeam_elev[k]);
            printf("Trovati Pochi Beam Elevazione %2d - num.: %3d",k,nbeam_elev[k]);
            ScrivoLog(16,errori);
            exit (1);
        }
    }                                                             /*end for*/
    sprintf(errori,"primi test passati");
    ScrivoLog(16,errori);


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
                sprintf(errori,"File Vecchio");
                ScrivoLog(16,errori);
                //return 0;
            }
            /*----------------------------
              |  aggiorno la data nel file |
              ----------------------------*/
            rewind(f_aus);
            fwrite(&old_data_header.norm.maq.acq_date,4,1,f_aus);
            fclose(f_aus);
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
    return 1;
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


    int tipo_dati_richiesti = INDEX_Z;
    //--------lettura volume------
    int ier = read_dbp_SP20((char*)nome_file,vol_pol,&old_data_header,
                            tipo_dati_richiesti,nbeam_elev);

    // ----- TEMPO E LOG------
    //  -----Scrivolog --------
    printf("fatta lettura\n");
    ScrivoLog(5,nome_file);

#ifdef TIME
    prendo_tempo();
#endif
    // ----- FINE TEMPO E LOG------


    //--test legato a argv[2]

    //  ----- Test sul volume test_file.......  --------
    int ier_test=test_file(file_type);

    printf("ier -- test  %d  %d\n",ier,ier_test);

    return ier == OK && ier_test;
}


void ScrivoLog(int i, const char* stringa)
    /*======================scrivo log=================================================================*/
{
    static FILE *log = 0;
    if (!log)
    {
        const char* log_fname = getenv("LOG_FILE");
        // If LOG_FILE is not set, just do not log
        if (log_fname == NULL)
            return;
        log=fopen(log_fname, "a");
        if(log==NULL)
        {
            printf(" impossibile aprire il file di log \n");
            exit(1);
        }
    }
    switch (i)
    {
        case  0:
            fprintf(log,"%s -- Lancio Programma\n",PrendiOra());
            break;
        case  1:
            fprintf(log,"-----------------------------------------------------------------\n");
            fprintf(log,"Flag di Compilazione\n");
#ifdef BOLOGNA
            fprintf(log," BOLOGNA ");
#endif
#ifdef SPC
            fprintf(log," SPC ");
#endif
#ifdef WRITE_DBP
            fprintf(log," WRITE_DBP ");
#endif
#ifdef TIME
            fprintf(log," TIME ");
#endif
#ifdef WRITE_DBP_REORDER
            fprintf(log," WRITE_DBP_REORDER ");
#endif
#ifdef DECLUTTER
            fprintf(log," DECLUTTER ");
#endif
#ifdef Z_AVERAGE
            fprintf(log," Z_AVERAGE ");
#endif
#ifdef R_AVERAGE
            fprintf(log," R_AVERAGE ");
#endif
#ifdef Z_LOWRIS
            fprintf(log," Z_LOWRIS ");
#endif
#ifdef ANAPROP
            fprintf(log," ANAPROP");
#endif
#ifdef SHORT
            fprintf(log," SHORT");
#endif
#ifdef MEDIUM
            fprintf(log," MEDIUM");
#endif
#ifdef STATIC
            fprintf(log," STATIC");
#endif
#ifdef BEAMBLOCKING
            fprintf(log," BEAMBLOCKING");
#endif
#ifdef QUALITY
            fprintf(log," QUALITY");
#endif

            fprintf(log,"\n");
            fprintf(log,"-----------------------------------------------------------------\n");
            fprintf(log,"Variabili d'Ambiente\n");
            fprintf(log,"LISTA_FILE = %s\n",getenv("LISTA_FILE"));
            fprintf(log,"LAST_FILE = %s\n",getenv("LAST_FILE"));
            fprintf(log,"ANAP_STAT_FILE = %s\n",getenv("ANAP_STAT_FILE"));
            fprintf(log,"BLOC_STAT_FILE = %s\n",getenv("BLOC_STAT_FILE"));
            fprintf(log,"ELEV_STAT_FILE = %s\n",getenv("ELEV_STAT_FILE"));
            fprintf(log,"DIR_OUT_PP_BLOC = %s\n",getenv("DIR_OUT_PP_BLOC"));
            fprintf(log,"FILE_DEM_SPC = %s\n",getenv("FILE_DEM_SPC"));
            fprintf(log,"FILE_DEM_GAT = %s\n",getenv("FILE_DEM_GAT"));
            fprintf(log,"DIR_QUALITY = %s\n",getenv("DIR_QUALITY"));
            fprintf(log,"FIRST_LEVEL_FILE = %s\n",getenv("FIRST_LEVEL_FILE"));
            fprintf(log,"OUTPUT_Z_DIR = %s\n",getenv("OUTPUT_Z_DIR"));
            fprintf(log,"OUTPUT_RAIN_DIR = %s\n",getenv("OUTPUT_RAIN_DIR"));
            fprintf(log,"OUTPUT_Z_LOWRIS_DIR = %s\n",getenv("OUTPUT_Z_LOWRIS_DIR"));
            fprintf(log,"-----------------------------------------------------------------\n");
            break;
        case  2:
            fprintf(log,"%s -- Apertura File Lista %s\n",PrendiOra(),getenv("LISTA_FILE"));
            break;
        case  3:
            fprintf(log,"%s -- Errore Apertura File Lista%s\n",
                    PrendiOra(),getenv("LISTA_FILE"));
            break;
        case  4:
            fprintf(log,"%s -- Apertura File Dati %s\n",PrendiOra(),stringa);
            break;
        case  5:
            fprintf(log,"%s -- Lettura File Dati\n",PrendiOra());
            break;
        case  6:
            fprintf(log,"%s -- Scrittura File Ordinato %s\n",PrendiOra(),stringa);
            break;
        case  7:
            fprintf(log,"%s -- Cancellazione Clutter e Propagazione Anomala\n",PrendiOra());
            break;
        case  8:
            fprintf(log,"%s -- Scrittura File Polare Ripulito %s\n",PrendiOra(),stringa);
            break;
        case  9:
            fprintf(log,"%s -- Creazione Matrice Cartesiana \n",PrendiOra());
            break;
        case 13:
            fprintf(log,"%s -- Estrazione Precipitazione 1X1\n",PrendiOra());
            break;
        case 14:
            fprintf(log,"%s -- Scrittura File Precipitazione 1X1 %s\n",PrendiOra(),stringa);
            break;
        case 15:
            fprintf(log,"%s -- Fine Programma\n",PrendiOra());
            fclose(log);
            break;
        case 16:
            fprintf(log,"%s -- %s\n",PrendiOra(),stringa);
            break;
        case 17:
            fprintf(log,"%s -- %s scrittura file qualita'\n",PrendiOra(),stringa);
            break;
        case 18:
            if (stringa==NULL) fprintf(log,"%s -- %s manca argomento necessario al programma \n",PrendiOra(),stringa);
            else
                fprintf(log,"%s -- argomento passato al programma %s \n",PrendiOra(),stringa);
            break;

    }

}                        /*end funzione ScrivoLog(i,stringa)*/

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

    if(time1 == 0){
        time1=time(&time1);
        time_tot = time1;
    }
    time2 = time(&time2);
    printf(" tempo parziale %ld ---- totale %ld\n",time2-time1, time2-time_tot);
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


