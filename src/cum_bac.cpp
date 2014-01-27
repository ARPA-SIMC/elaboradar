#include "cum_bac.h"
#include "logging.h"
#include "utils.h"

#include <cstring>

#ifdef __cplusplus
extern "C" {
#endif
// libreria radar
//#include <radar_parameter.h>
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include <func_Z_R.h>
#ifdef __cplusplus
}
#endif

#include <qual_par.h>

//Definizioni per test_file
#define NUM_MIN_BEAM 200
#define SHORT_DEC         0
#define SHORT_FULL_VOLUME 1
#define SHORT_HAIL        2
#define MEDIUM_PULSE      3

//Definizioni per statistica anap
#define STEP_STAT_ANAP_RANGE  40 /*dim range griglia per stat anap*/
#define STEP_STAT_ANAP_AZ     25 /*dim azim griglia per stat anap*/
#define N_MIN_BIN           500 /*--- numero minimo di celle presenti in un
                  settore per la statistica            ---*/

// Soglie algoritmi
#define MAX_DIF_OR 30            /* differenzio limiti controllo anap      */
#define MIN_VALUE_OR -10         /* a seconda che sia alla prima o success.*/
#define MAX_DIF_NEXT_OR 15       /*/ elevazione                             */
#define MIN_VALUE_NEXT_OR 0
#define THR_CONT_ANAP 1 /* limite in numero occorrenze anaprop sul raggio dopo i 30 km per non togliere =*/
#define OVERBLOCKING 51 /* minimo BB non accettato*/

// anaprop
#define LIMITE_ANAP 240/* LIMITE in numero bins per cambiare controllo anaprop*/

#define DTOR  M_PI/180. /* esternalizzo?*/ //fattore conversione gradi-radianti
#define CONV_RAD 360./4096.*DTOR  // fattore conversione unità angolare radar-radianti

// parametri ereditati da programma beam blocking:numero elevazioni da programma beam blocking ; le matrici ivi definite considerano questo
#define NSCAN 6

/// This needs to be a global variable, as it is expected by libsp20
int elev_array[NEL];


CUM_BAC::CUM_BAC()
    : do_quality(false), do_beamblocking(false), do_declutter(false), do_bloccorr(false)
{
    logging_category = log4c_category_get("radar.cum_bac");

    t_ground=NODATAVPR;
    chisqfin=100; //???puo' essere def in anal
    rmsefin=100;
    ncv=0;ncs=0;np=0;
    htbb=-9999.; hbbb=-9999.;
    memset(vol_pol,0,sizeof(vol_pol));
    memset(cart,0,sizeof(cart));
    memset(cartm,0.,sizeof(cartm));
    memset(z_out,0,sizeof(z_out));
    memset(nbeam_elev,0,sizeof(nbeam_elev));
    memset (first_level,0,sizeof(first_level));
    memset (first_level_static,0,sizeof(first_level_static));
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

    for (int k=0; k<NUM_AZ_X_PPI*MAX_BIN;k++ ){
      lista_conv[k][0]=-999;
      lista_conv[k][1]=-999;
      lista_bckg[k][0]=-999;
      lista_bckg[k][1]=-999;
    }

    memset(stratiform,0,sizeof(stratiform));

    memset(bb_first_level,0,sizeof(bb_first_level));

    memset(stat_anap_tot,0,sizeof(stat_anap_tot));
    memset(stat_anap,0,sizeof(stat_anap));
    memset(stat_bloc,0,sizeof(stat_bloc));
    memset(stat_elev,0,sizeof(stat_elev));


    //-----  FINE INIZIALIZZAZIONI---------//
}

void CUM_BAC::setup_elaborazione(const char* nome_file, const char* sito)
{
    /*------------------------------------------
      | rimozione propagazione anomala e clutter |
      ------------------------------------------*/
    LOG_INFO("%s -- Cancellazione Clutter e Propagazione Anomala", nome_file);

    assets.configure(sito, old_data_header.norm.maq.acq_date);

    // --- ricavo il mese x definizione first_level e  aMP bMP ---------
    //definisco stringa data in modo predefinito
    Time = NormalizzoData(old_data_header.norm.maq.acq_date);
    tempo = gmtime(&Time);
    month=tempo->tm_mon+1;

    // scrivo la variabile char date con la data in formato aaaammgghhmm
    sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
            tempo->tm_mday,tempo->tm_hour, tempo->tm_min);

    // ------definisco i coeff MP in base alla stagione( mese) che servono per calcolo VPR e attenuazione--------------
    if ( month > 4 && month < 10 )  {
        aMP=aMP_conv;
        bMP=bMP_conv;
    }
    else {
        aMP=aMP_strat;
        bMP=bMP_strat;

    }
    MP_coeff[0]=(unsigned char)(aMP/10);
    MP_coeff[1]=(unsigned char)(bMP*10);
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
        LOG_INFO("Numero beam presenti: %4d -- elevazione %d", nbeam_elev[k], k);

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

    //  ----- Test sul volume test_file.......  --------
    if (!test_file(file_type))
    {
        LOG_ERROR("test_file failed");
        return false;
    }

    return ier == OK;
}

int CUM_BAC::elabora_dato()
{
    int i,l,k;
    int test_an;
    int el_inf,el_up;
    float bin_low,bin_high,bin_low_low,cont_anap;
    float fondo_scala,elevaz;
    unsigned char flag_anap;

    //-------------calcolo fondo scala   ------------
    flag_anap = 1;
    fondo_scala = BYTEtoDB(flag_anap);/*-19.7 dBZ*/


    //-------------leggo mappa statica ovvero first_level (funzione leggo_first_level)------------
    leggo_first_level();

    //-------------se definita qualita' leggo dem e altezza fascio (funzioni legg_dem e leggo_hray)(mi servono per calcolare qualità)
    if (do_quality)
    {
        leggo_dem();
        leggo_hray();
    }

    //------------se definito DECLUTTER , non rimuovo anap e riscrivo  volume polare facedndo declutter solo con mappa statica.... ancora valido?

    if (do_declutter)
    {
        for(i=0; i<NUM_AZ_X_PPI; i++)
        {
            for(k=0; k<vol_pol[0][i].b_header.max_bin; k++)
            {
                //---assegno el_inf a mappa statica
                el_inf = first_level_static[i][k];
                //---ricopio valori a mappa statica sotto
                for(l=0; l<=el_inf; l++)
                {
                    vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
                    //------------se definito BEAM BLOCKING e non definito BLOCNOCORR (OPZIONE PER non correggere il beam blocking a livello di mappa statica PUR SAPENDO QUANT'È)
                    if (do_beamblocking && do_bloccorr)
                    {
                        vol_pol[l][i].ray[k]=DBtoBYTE(BYTEtoDB(vol_pol[l][i].ray[k])-10*log10(1.-(float)beam_blocking[i][k]/100.));
                        //vol_pol[l][i].ray[k]=vol_pol[l][i].ray[k]+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5);
                    }

                }
                if (do_quality)
                    elev_fin[i][k]=el_inf;
            }
        }
        return 0;
    }

    //------------se non definito DECLUTTER inizio rimozione propagazione anomala al livello mappa dinamica e elaborazioni accessorie


    /* 26-5-2004 : se sono alla 1 o successive elevazioni
       e range > 60 km cambio le soglie, in modo
       da evitare di riconoscere come anaprop una pioggia shallow
       Il criterio diventa: - se la differenza tra Z all'elevazione più bassa della
       corrente e la Z corrente è <10 dbZ allora
       rendo inefficaci i limiti di riconoscimento anaprop. */

    //--------ciclo sugli azimut e bins per trovare punti con propagazione anomala----------------

    for(i=0; i<NUM_AZ_X_PPI; i++)
    {
        flag_anap = 0;
        cont_anap=0;// aggiunto per risolvere problema di uso con preci shallow
        for(k=0; k<vol_pol[0][i].b_header.max_bin; k++)
            //------------- incremento statistica tot ------------------
        {
            stat_anap_tot[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
            // ------------assegno l'elevazione el_inf a first_level e elev_fin a el_inf---------
            el_inf = first_level[i][k];

            if (do_quality)
                elev_fin[i][k]=el_inf;

            // ------------assegno a el_up il successivo di el_inf e se >=NEL metto bin_high=fondo_scala
            el_up = el_inf +1;
            if ( el_up >= NEL )
                bin_high = fondo_scala ;

            // ------------assegno  bin_low_low (cioè il valore sotto il bin base)
            if (el_inf>0) bin_low_low=BYTEtoDB(vol_pol[el_inf-1][i].ray[k]);
            else  bin_low_low=fondo_scala+1;

            // ------------assegno bin_low bin_high anche
            bin_low  = BYTEtoDB(vol_pol[el_inf][i].ray[k]);
            bin_high = BYTEtoDB(vol_pol[el_up][i].ray[k]);

            //------------assegno le soglie per anaprop : se sono oltre 60 km e se la differenza tra il bin sotto il base e quello sopra <10 non applico test (cambio i limiti per renderli inefficaci)
            MAX_DIF=MAX_DIF_OR;
            MAX_DIF_NEXT=MAX_DIF_NEXT_OR;
            MIN_VALUE=MIN_VALUE_OR;
            MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;
            //----------questo serviva per evitare di tagliare la precipitazione shallow ma si dovrebbe trovare un metodo migliore p.es. v. prove su soglia
            if((el_inf>=1)&&(k>LIMITE_ANAP)&&(bin_low_low-bin_low<10)) //-----------ANNULLO EFFETTO TEST ANAP
            {
                MAX_DIF_NEXT=BYTEtoDB(255);
                MAX_DIF=BYTEtoDB(255);
                MIN_VALUE=BYTEtoDB(0);
                MIN_VALUE_NEXT= BYTEtoDB(0);  }

            // ------------separo i diversi casi x analisi anaprop: ho dati sia al livello base che sopra o no  e ho trovato anaprop in precedenza sul raggio o no
            if (cont_anap> THR_CONT_ANAP || k < 80  )
                test_an=(bin_low > fondo_scala && bin_high >= fondo_scala );
            else
                test_an=(bin_low > fondo_scala && bin_high > fondo_scala );

            //------------------ se ho qualcosa sia al livello base che sopra allora effettuo il confronto-----------------
            //if(bin_low > fondo_scala && bin_high >= fondo_scala ) // tolto = per bin_high
            if(test_an )
            {
                //------------------ se ho trovato anap prima nel raggio cambio le soglie le abbasso)-----------------
                if(flag_anap)
                {
                    //-----------caso di propagazione anomala ---------
                    if(bin_low-bin_high >= MAX_DIF_NEXT || bin_high <= MIN_VALUE_NEXT )
                    {

                        //---------assegno l'indicatore di presenza anap nel raggio e incremento statistica anaprop, assegno matrici che memorizzano anaprop  e elevazione_finale e azzero beam blocking perchè ho cambiato elevazione
                        flag_anap = 1;
                        cont_anap=cont_anap+1;

                        //--------ricopio valore a el_up su tutte elev inferiori--------------
                        for(l=0; l<el_up; l++) {
                            vol_pol[l][i].ray[k]=vol_pol[el_up][i].ray[k];
                        } //

                        //--------azzero beam_blocking ( ho cambiato elevazione, non ho disponible il bbeam blocking all' elev superiore)--------------
                        if (do_beamblocking)
                            beam_blocking[i][k]=0;/* beam blocking azzerato */

                        //--------------------incremento la statitica anaprop e di cambio elevazione-------------
                        stat_anap[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
                        if (el_up > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++; //incremento la statistica cambio elevazione

                        //-------------------memorizzo dati di qualita '-------------
                        if (do_quality)
                        {
                            dato_corrotto[i][k]=ANAP_YES;/*  risultato test: propagazione anomala*/
                            elev_fin[i][k]=el_up;
                        }
                    }
                    else
                        //-----non c'è propagazione anomala:ricopio su tutte e elevazioni il valore di el_inf e correggo il beam blocking e incremento la statistica beam_blocking, assegno matrice anaprop a 0 nel punto e assegno a 0 indicatore anap nel raggio-----------
                    {
                        flag_anap = 0;
                        if (do_beamblocking && do_bloccorr)
                        {
                            vol_pol[el_inf][i].ray[k]=DBtoBYTE(BYTEtoDB(vol_pol[l][i].ray[k])-10*log10(1.-(float)beam_blocking[i][k]/100.));
                            //    vol_pol[l][i].ray[k]=vol_pol[l][i].ray[k]+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5); //correggo beam blocking
                            stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]= stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]+ beam_blocking[i][k]; // incremento statistica beam blocking
                        }
                        for(l=0; l<=el_up; l++){
                            vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k]; // assegno a tutti i bin sotto el_inf il valore a el_inf (preci/Z a el_inf nella ZLR finale)

                        }
                        if (do_quality)
                        {
                            dato_corrotto[i][k]=ANAP_OK;/* matrice risultato test: no propagazione anomala*/
                            elev_fin[i][k]=el_inf;
                        }
                        if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;//incremento la statistica cambio elevazione

                    }
                }
                //------------------se invece flag_anap == 0 cioè non ho ancora trovato anaprop nel raggio ripeto le operazioni ma con i limiti più restrittivi (controlla elev_fin)------------------

                else
                {
                    if(bin_low-bin_high >= MAX_DIF || bin_high <= MIN_VALUE  )

                    {
                        //--------ricopio valore a el_up su tutte elev inferiori--------------
                        for(l=0; l<el_up; l++){
                            vol_pol[l][i].ray[k]=1;
                            vol_pol[l][i].ray[k]=vol_pol[el_up][i].ray[k];  //ALTERN
                        }
                        //---------assegno l'indicatore di presenza anap nel raggio e incremento statistica anaprop, assegno matrici che memorizzano anaprop e elevazione_finale e azzero beam blocking perchè ho cambiato elevazione
                        flag_anap = 1;
                        cont_anap=cont_anap+1;
                        stat_anap[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
                        if (do_quality)
                        {
                            dato_corrotto[i][k]=ANAP_YES;/*matrice risultato test: propagazione anomala*/
                            elev_fin[i][k]=el_up;
                        }
                        if (el_up > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;//incremento la statistica cambio elevazione
                        if (do_beamblocking)
                            beam_blocking[i][k]=0;
                    }

                    //-----non c'è propagazione anomala:ricopio su tutte e elevazioni il valore di el_inf e correggo il beam blocking,  incremento la statistica beam_blocking, assegno matrice anaprop a 0 nel punto , assegno a 0 indicatore anap nel raggio, assegno elevazione finale e incremento statisica cambio elevazione se el_inf > first_level_static[i][k]-----------
                    else
                    {
                        for(l=0; l<=el_inf; l++)
                        {
                            vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
                            if (do_beamblocking && do_bloccorr)
                            {
                                vol_pol[l][i].ray[k]=DBtoBYTE(BYTEtoDB(vol_pol[l][i].ray[k])-10*log10(1.-(float)beam_blocking[i][k]/100.));
                                //vol_pol[l][i].ray[k]=vol_pol[l][i].ray[k]+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5);
                                stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]= stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]+ beam_blocking[i][k];
                            }
                        }

                        if (do_quality)
                        {
                            dato_corrotto[i][k]=ANAP_OK;
                            elev_fin[i][k]=el_inf;
                        }

                        if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;//incremento la statistica cambio elevazione
                        flag_anap = 0;

                    } /*endif test anaprop*/
                }/*endif flaganap*/
            }/*endif bin_low > fondo_scala && bin_high >= fondo_scala*/

            //----------------se al livello base non ho dato riempio con i valori di el_up tutte le elevazioni sotto (ricostruisco il volume) e assegno beam_blocking 0
            else if (bin_low < fondo_scala)
            {
                for(l=0; l<el_up; l++)
                {
                    vol_pol[l][i].ray[k]=vol_pol[el_up][i].ray[k];
                    if(!vol_pol[l][i].flag)
                    {
                        vol_pol[l][i].flag = 1;
                        vol_pol[l][i].b_header.alfa =(short)(i*.9/FATT_MOLT_AZ);
                        vol_pol[l][i].b_header.teta = elev_array[l];
                        vol_pol[l][i].b_header.tipo_gran=vol_pol[el_up][i].b_header.tipo_gran;
                        vol_pol[l][i].b_header.max_bin=vol_pol[el_up][i].b_header.max_bin;
                    }
                }
                //----------------controlli su bin_high nel caso in cui bin_low sia un no data per assegnare matrice anap  (dato_corrotto[i][k])
                if (do_quality)
                {
                    if (bin_high<fondo_scala)   dato_corrotto[i][k]=ANAP_NODAT;/*manca dato sotto e sopra*/
                    if (cont_anap< THR_CONT_ANAP )
                        test_an=(bin_high>=fondo_scala); //modificato per contemplare > o >=
                    else
                        test_an=(bin_high>fondo_scala);

                    // if (bin_high>=fondo_scala)  dato_corrotto[i][k]=ANAP_NOCONTROL;/*manca controllo (sotto non ho nulla)*/ //messo l=
                    if (test_an) dato_corrotto[i][k]=ANAP_NOCONTROL;
                    if (bin_high==fondo_scala) dato_corrotto[i][k]=ANAP_OK;/*non piove (oppure sono sopra livello preci...)*/
                }

                if (do_beamblocking)
                    beam_blocking[i][k]=0;
            }

            //----------------se bin_low == fondo_scala riempio matrice vol_pol con dato a el_inf (mi resta il dubbio di quest'if se seve o basti un else ) azzero matrice anap (dato ok)
            else if (bin_low == fondo_scala || bin_high <= fondo_scala)/* quel che resta da (bin_low > fondo_scala && bin_high >= fondo_scala) e (bin_low < fondo_scala) ; messo =per bin_high*/

            {

                for(l=0; l<el_inf; l++)//riempio con i valori di el_inf tutte le elevazioni sotto (ricostruisco il volume)
                {
                    vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
                    if(!vol_pol[l][i].flag)
                    {
                        vol_pol[l][i].flag = 1;
                        vol_pol[l][i].b_header.alfa =(short)(i*.9/FATT_MOLT_AZ);
                        vol_pol[l][i].b_header.teta = elev_array[l];  //perchè ridefinisce ??
                        vol_pol[l][i].b_header.tipo_gran=vol_pol[el_inf][i].b_header.tipo_gran;
                        vol_pol[l][i].b_header.max_bin=vol_pol[el_inf][i].b_header.max_bin;
                    }
                }

                if (do_quality)
                {
                    dato_corrotto[i][k]=ANAP_OK; // dubbio
                    elev_fin[i][k]=el_inf;
                }

                if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
            }
            /*-----------------------------------------------------------fine di tutti gli if-----------*/
            //-----finiti tutti i controlli assegno le varibili di qualita definitive: elevazione, quota calcolata sull'elevazione reale con propagazione standard , e quota relativa al suolo calcolata con elevazione nominale e propagazione da radiosondaggio.

            if (do_quality)
            {
                elevaz=(float)(vol_pol[elev_fin[i][k]][i].teta_true)*CONV_RAD;
                // elev_fin[i][k]=first_level_static[i][k];//da togliere
                quota[i][k]=(unsigned short)(quota_f(elevaz,k));
                quota_rel[i][k]=(unsigned short)(hray[k][elev_fin[i][k]]-dem[i][k]);/*quota sul suolo in m con elev nominale e prop da radiosondaggio (v. programma bloc_grad.f90)*/
            }
        }
    }

    LOG_INFO("elabora_dato completed");
    ScrivoStatistica();

    return 0;
}                         /*end funzione elabora_dato()*/

void CUM_BAC::leggo_first_level()
{
    FILE *file;

#ifdef STATIC
    /*-------------------
      Leggo mappa  statica
      -------------------*/
    // lettura dimensioni matrice mappa statica da file esterno
    int dim = assets.read_file_first_level_dim();
    //leggo mappa statica con dimensioni appena lette
    file = assets.open_file_first_level();
    for(int i=0; i<NUM_AZ_X_PPI; i++)
        fread(&first_level_static[i][0],dim,1,file);
    // copio mappa statica su matrice first_level
    memcpy(first_level,first_level_static,sizeof(first_level));
    fclose(file);
#endif

    if (do_beamblocking)
    {
        /*----------------------------
          Leggo file elevazioni per BB
          ----------------------------*/
        file = assets.open_file_first_level_bb_el();
        for(int i=0; i<NUM_AZ_X_PPI; i++)
            fread(&bb_first_level[i][0],MAX_BIN,1,file);
        fclose(file);
        /*------------------------
          Leggo file valore di BB
          ------------------------*/
        file = assets.open_file_first_level_bb_bloc();
        /* Se elevazione clutter statico < elevazione BB, prendi elevazione BB,
           altrimeti prendi elevazione clutter statico e metti a 0 il valore di BB*/
        for(int i=0; i<NUM_AZ_X_PPI; i++){   /*ciclo sugli azimut*/
            fread(&beam_blocking[i][0],MAX_BIN,1,file);

            for (int j=0; j<MAX_BIN; j++) /*ciclo sul range  */
            {
                if (do_bloccorr)
                {
                    if (first_level_static[i][j]<=bb_first_level[i][j])
                        first_level[i][j]=bb_first_level[i][j];
                    else
                    {  beam_blocking[i][j]=0;
                        first_level[i][j]=first_level_static[i][j]; }
                } else {
                    if (first_level_static[i][j]>bb_first_level[i][j])
                        beam_blocking[i][j]=0;
                    if (first_level_static[i][j]<bb_first_level[i][j])
                        beam_blocking[i][j]=OVERBLOCKING;
                }
            }
        }
        fclose(file);
    }

    /*-------------------------------
      patch per espandere il clutter
      -------------------------------*/
#ifdef MEDIUM
    {
        unsigned char first_level_tmp[NUM_AZ_X_PPI][MAX_BIN];
        int k;
        memcpy(first_level_tmp,first_level,sizeof(first_level));
        for (int i=NUM_AZ_X_PPI; i<800; i++)
        {
            for (int j=0; j<MAX_BIN; j++)
            {
                for (k=i-1; k<i+2; k++)
                    if(first_level[i%NUM_AZ_X_PPI][j] < first_level_tmp[k%NUM_AZ_X_PPI][j])
                        first_level[i%NUM_AZ_X_PPI][j]=first_level_tmp[k%NUM_AZ_X_PPI][j];
            }
        }
    }
#endif
}

void CUM_BAC::leggo_hray( )
{
    FILE *file;

    /*--------------------------
      Leggo quota centro fascio
      --------------------------*/
    file = assets.open_file_hray();
    fscanf(file,"%f ",&dtrs);
    for(int i=0; i<MAX_BIN; i++){
        for(int j=0; j<NSCAN;j++)
            fscanf(file,"%f ",&hray[i][j]);
    }
    fclose(file);

    file = assets.open_file_hray_inf();
    fscanf(file,"%f ",&dtrs);
    for(int i=0; i<MAX_BIN; i++){
        for(int j=0; j<NSCAN;j++)
            fscanf(file,"%f ",&hray_inf[i][j]);
    }
    fclose(file);

    return  ;
}

void CUM_BAC::leggo_dem()
{
    /*---------------------
      Leggo dem
      ---------------------*/
    FILE *file = assets.open_file_dem();
    for (int i=0; i<MAX_BIN; i++){
        for (int j=0; j<NUM_AZ_X_PPI;j++)
            fscanf(file,"%f ",&dem[j][i]);
    }
    fclose(file);
    return ;
}

//------------funzione quota_f-----------------------
//---------funzione che calcola la quota in metri del centro del fascio-----------------------
//--------distanza=k*dimensionecella +semidimensionecella in metri ----------------------
//--------quota=f(distinkm, rstinkm, elevazinrad) in metri
float CUM_BAC::quota_f(float elevaz, int k) // quota funzione di elev(radianti) e range
{
    float dist;
    float q_st;
    dist=k*(size_cell[old_data_header.norm.maq.resolution])+(size_cell[old_data_header.norm.maq.resolution])/2.;
    q_st=(sqrt(pow(dist/1000.,2)+(rst*rst)+2.0*dist/1000.*rst*sin(elevaz))-rst)*1000.; /*quota in prop. standard da elevazione reale  */;

    return q_st;
}

void CUM_BAC::ScrivoStatistica()
{
    int az,ran;
    unsigned char statistica[DIM1_ST][DIM2_ST];
    unsigned char statistica_bl[DIM1_ST][DIM2_ST];
    unsigned char statistica_el[DIM1_ST][DIM2_ST];

    LOG_INFO("scrivo statistica");
    memset(statistica,255,DIM1_ST*DIM2_ST);
    memset(statistica_bl,255,DIM1_ST*DIM2_ST);
    memset(statistica_el,255,DIM1_ST*DIM2_ST);

    for(az=0; az<DIM1_ST; az++)
        for(ran=0; ran<DIM2_ST; ran++)

            if(stat_anap_tot[az][ran] >= N_MIN_BIN)
            {
                statistica[az][ran] =
                    (unsigned char)((stat_anap[az][ran]*100)/stat_anap_tot[az][ran]);
                statistica_bl[az][ran] =
                    (unsigned char)((stat_bloc[az][ran]*100)/stat_anap_tot[az][ran]);
                statistica_el[az][ran] =
                    (unsigned char)((stat_elev[az][ran]*100)/stat_anap_tot[az][ran]);
            }

    FILEFromEnv f_stat;

    if (f_stat.open_from_env("ANAP_STAT_FILE", "a"))
    {
        fwrite(date,12,1,f_stat);
        fwrite(statistica,DIM1_ST*DIM2_ST,1,f_stat);
    }

    if (f_stat.open_from_env("BLOC_STAT_FILE", "a"))
    {
        fwrite(date,12,1,f_stat);
        fwrite(statistica_bl,DIM1_ST*DIM2_ST,1,f_stat);
    }

    if (f_stat.open_from_env("ELEV_STAT_FILE", "a"))
    {
        fwrite(date,12,1,f_stat);
        fwrite(statistica_el,DIM1_ST*DIM2_ST,1,f_stat);
    }

    return ;
}

FILE *CUM_BAC::controllo_apertura (const char *nome_file, const char *content, const char *mode)
{
    FILE *file;
    int ier_ap=0;

    //printf("controllo apertura %s\n",nome_file);

    if (strcmp(mode,"r") == 0)
        ier_ap=access(nome_file,R_OK);
    else ier_ap=0;

    if (!ier_ap) {
        if (strcmp(mode,"r") == 0) file = fopen(nome_file,"r");
        else file=fopen(nome_file,"w");
    }
    else
    {
        LOG_ERROR("Errore Apertura %s %s", content, nome_file);
        exit(1);
    }
    if (file == NULL) {
        LOG_ERROR("Errore Apertura %s %s", content, nome_file);
        exit(1);
    }
    return(file);
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


