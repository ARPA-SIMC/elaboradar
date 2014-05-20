#include "cum_bac.h"
#include "logging.h"
#include "utils.h"
#include "site.h"
#include "volume/sp20.h"
#include "volume/odim.h"
#include "cylindrical.h"
#include "steiner.h"
#include "viz.h"
#include "interpola_vpr.h"
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <math.h>
#include <iostream>
#include <unistd.h>
#include "setwork.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <func_Z_R.h>
#ifdef __cplusplus
}
#endif

#include <func_Q3d.h>

#include <qual_par.h>
#include <par_class.h>

#ifdef NEL
#undef NEL
#endif

//Definizioni per test_file
#define NUM_MIN_BEAM 200
#define SHORT_DEC         0
#define SHORT_FULL_VOLUME 1
#define SHORT_HAIL        2
#define MEDIUM_PULSE      3
#define SHORT_212	  4

// Soglie algoritmi
#define MAX_DIF_OR 30            /* differenzio limiti controllo anap      */
#define MIN_VALUE_OR -10         /* a seconda che sia alla prima o success.*/
#define MAX_DIF_NEXT_OR 15       /*/ elevazione                             */
#define MIN_VALUE_NEXT_OR 0
#define THR_CONT_ANAP 1 /* limite in numero occorrenze anaprop sul raggio dopo i 30 km per non togliere =*/
#define OVERBLOCKING 51 /* minimo BB non accettato*/
#define SOGLIA_TOP 20 // soglia per trovare top
#define THRES_ATT 0 /* minimo valore di Z in dBZ per calcolare att rate */
#define MISSING 0 /*valore mancante*/

//Definizioni geometriche
#define AMPLITUDE 0.9 /* esternalizzo?*/ // ampiezza fascio radar

// anaprop
#define LIMITE_ANAP 240/* LIMITE in numero bins per cambiare controllo anaprop*/

#define DTOR  M_PI/180. /* esternalizzo?*/ //fattore conversione gradi-radianti
#define CONV_RAD 360./4096.*DTOR  // fattore conversione unità angolare radar-radianti

using namespace std;

namespace cumbac {

namespace {
/**
 *  @brief funzione che legge la quota del centro fascio e del limite inferiore del fascio da file
 *  @details legge la quota del centro fascio e del limite inferiore del fascio da file e li memorizza nei vettori hray_inf e hray
 *  @return 
 */
struct HRay
{
    static const int NSCAN = 6;

    float* hray;
    // distanza temporale radiosondaggio
    float dtrs;

    HRay() : hray(0) { }
    ~HRay() { if (hray) delete[] hray; }

    float* operator[](unsigned idx) { return hray + idx * NSCAN; }
    const float* operator[](unsigned idx) const { return hray + idx * NSCAN; }

    void load_hray(Assets& assets)
    {
        // quota centro fascio in funzione della distanza e elevazione
        load_file(assets.open_file_hray());
    }
    void load_hray_inf(Assets& assets)
    {
        // quota limite inferiore fascio in funzione della distanza e elevazione
        load_file(assets.open_file_hray_inf());
    }


private:
    void load_file(FILE* file)
    {
        if (hray)
        {
            delete[] hray;
            hray = 0;
        }
        hray = new float[MAX_BIN * NSCAN];

        /*--------------------------
          Leggo quota centro fascio
          --------------------------*/
        fscanf(file,"%f ",&dtrs);
        for(int i=0; i<MAX_BIN; i++){
            for(int j=0; j<NSCAN;j++)
                fscanf(file,"%f ", hray + i * NSCAN + j);
        }
        fclose(file);
    }
};
}

GridStats::GridStats()
{
}

GridStats::~GridStats()
{
    if (stat_anap) delete[] stat_anap;
    if (stat_tot) delete[] stat_tot;
    if (stat_bloc) delete[] stat_bloc;
    if (stat_elev) delete[] stat_elev;
}

void GridStats::init(const Volume<double>& volume)
{
    size_az = volume.scan(0).beam_count / step_stat_az + 1;
    size_beam = volume.scan(0).beam_size / step_stat_range + 1;

    stat_anap = new unsigned[size_az * size_beam];
    stat_tot = new unsigned[size_az * size_beam];
    stat_bloc = new unsigned[size_az * size_beam];
    stat_elev = new unsigned[size_az * size_beam];

    for (unsigned i = 0; i < size_az * size_beam; ++i)
        stat_anap[i] = stat_tot[i] = stat_bloc[i] = stat_elev[i] = 0;
}

CUM_BAC::CUM_BAC(const char* site_name, bool medium, unsigned max_bin)
    : MyMAX_BIN(max_bin), site(Site::get(site_name)),
      do_medium(medium),
      elev_fin(volume, load_info),
      calcolo_vpr(0),
      first_level(MyMAX_BIN), first_level_static(MyMAX_BIN),
      bb_first_level(MyMAX_BIN), beam_blocking(MyMAX_BIN),dem(MyMAX_BIN),
      qual(0), cappi(MyMAX_BIN)
{
    logging_category = log4c_category_get("radar.cum_bac");
}

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

CUM_BAC::~CUM_BAC()
{
    if (qual) delete qual;
    if (calcolo_vpr) delete calcolo_vpr;
}

void CUM_BAC::setup_elaborazione(const char* nome_file)
{
    /*------------------------------------------
      | rimozione propagazione anomala e clutter |
      ------------------------------------------*/
    LOG_INFO("%s -- Cancellazione Clutter e Propagazione Anomala", nome_file);

    assets.configure(site, load_info.acq_date);

    grid_stats.init(volume);

    // --- ricavo il mese x definizione first_level e  aMP bMP ---------
    //definisco stringa data in modo predefinito
    time_t Time = NormalizzoData(load_info.acq_date);
    struct tm* tempo = gmtime(&Time);
    int month=tempo->tm_mon+1;

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

    dato_corrotto.resize(400, volume.max_beam_size());
    dato_corrotto.fill(0);
    quota.resize(400, volume.max_beam_size());
    quota.fill(0);

    compute_top();

    //--------------se definito VPR procedo con ricerca t_ground che mi serve per classificazione per cui la metto prima-----------------//
    if (do_vpr) calcolo_vpr = new CalcoloVPR(*this);
    LOG_INFO(" End setup_elaborazione");
}

void CUM_BAC::compute_top()
{
    top.resize(400, volume.max_beam_size());
    top.fill(0);
    for (unsigned l=0; l<volume.size(); l++)
    {
        const auto& scan = volume.scan(l);
        const auto& scan_info = load_info.scan(l);
        for (int i=0; i<NUM_AZ_X_PPI; i++)
        {
            const double elevaz = scan_info.get_elevation_rad(i); //--- elev reale
            for (unsigned k = 0; k < scan.beam_size; ++k)
                if (scan.get(i, k) > SOGLIA_TOP)
                    top(i, k) = (unsigned char)((quota_f(elevaz, k))/100.); //top in ettometri
        }
    }
}

bool CUM_BAC::test_file(int file_type)
{
    int n_elev = -1;
    int expected_size_cell = 0;// != volume.resolution?

    //--- switch tra tipo di file per definire nelev = elevazioni da testare e la risoluzione

    switch (file_type)
    {
        case SHORT_DEC:
            if (!load_info.declutter_rsp)
            {
                LOG_WARN("File Senza Declutter Dinamico--cos' è???");
                return false;
            }
            expected_size_cell = 250;
            n_elev=4;
            break;
            //------------se tipo =1 esco
        case SHORT_FULL_VOLUME://-----??? DUBBIO
            if (load_info.declutter_rsp)
            {
                LOG_WARN("File con Declutter Dinamico");
                return false;
            }
            expected_size_cell = 250;
            n_elev=4;
            break;
        case SHORT_HAIL://-----??? DA BUTTARE NON ESISTE PIÙ
            expected_size_cell = 250;
            n_elev=3;
            LOG_INFO("CASO SHORT_HAIL");
            break;
        case MEDIUM_PULSE:
            expected_size_cell = 1000;
            n_elev=4;
            LOG_INFO("CASO MEDIO OLD");
            break;
        case SHORT_212://----- CORRISPONDE A VOL_NEW - da questo si ottengono il corto e il medio
            if (!load_info.declutter_rsp)
            {
                LOG_WARN("File senza Declutter Dinamico");
                return false;
            }
            expected_size_cell = 250;
            n_elev=4;
            LOG_INFO("CASO SHORT_212");
            break;
    }

    //----------se la risoluzione del file è diversa da quella prevista dal tipo_file dà errore ed esce (perchè poi probabilmente le matrici sballano ?)
    if (load_info.size_cell != expected_size_cell)
    {
        LOG_ERROR("File Risoluzione/size_cell Sbagliata %f", (double)load_info.size_cell);
        return false;
    }
    //------eseguo test su n0 beam  sulle prime 4 elevazioni, se fallisce  esco ------------

    for (int k = 0; k < n_elev; k++) /* testo solo le prime 4 elevazioni */
    {
        volume::PolarScanLoadInfo info = load_info.scan(k);
        LOG_INFO("Numero beam presenti: %4d -- elevazione %d", info.count_rays_filled(), k);

        if (info.count_rays_filled() <  NUM_MIN_BEAM)
            // se numero beam < numero minimo---Scrivolog ed esco !!!!!!!!!!!!!!!!!!!
        {
            //---Scrivolog!!!!!!!!!!!!!!!!!!!
            LOG_ERROR("Trovati Pochi Beam Elevazione %2d - num.: %3d",k,info.count_rays_filled());
            return false;
        }
    }                                                             /*end for*/

    //--------verifico la presenza del file contenente l'ultima data processata-------
    bool is_new = assets.save_acq_time(load_info.acq_date);
    if (!is_new)
        LOG_WARN("File Vecchio");

    // ------- se ok status di uscita:1
    return true;
}

bool CUM_BAC::read_sp20_volume(const char* nome_file, int file_type)
{
    LOG_INFO("Reading %s for site %s and file type %d", nome_file, site.name.c_str(), file_type);

    volume::SP20Loader loader(site, do_medium, do_clean, MyMAX_BIN);
    loader.load_info = &load_info;
    loader.vol_z = &volume;
    loader.load(nome_file);

    elev_fin.init();

    /*
    printf("fbeam ϑ%f α%f", volume.scan(0)[0].teta, volume.scan(0)[0].alfa);
    for (unsigned i = 0; i < 20; ++i)
        printf(" %d", (int)volume.scan(0).get_raw(0, i));
    printf("\n");
    */

    //  ----- Test sul volume test_file.......  --------
    if (!test_file(file_type))
    {
        LOG_ERROR("test_file failed");
        return false;
    }

    return true;
}

bool CUM_BAC::read_odim_volume(const char* nome_file, int file_type)
{
    LOG_INFO("Reading %s for site %s and file type %d", nome_file, site.name.c_str(), file_type);

    volume::ODIMLoader loader(site, do_medium, do_clean, MyMAX_BIN);
    loader.load_info = &load_info;
    loader.vol_z = &volume;
    loader.load(nome_file);

    elev_fin.init();

    /*
    printf("fbeam ϑ%f α%f", this->volume.scan(0)[0].teta, this->volume.scan(0)[0].alfa);
    for (unsigned i = 0; i < 20; ++i)
        printf(" %d", (int)this->volume.scan(0).get_raw(0, i));
    printf("\n");
    */

    /*
    int numRaggi»···»···= scan->getNumRays();
    NUM_AZ_X_PPI

    NEL

        se due scan per stessa elecvazione, prendo il primo

        guardare se il passo di griglia è 0.9 o dare errore
        sennò prendere il beam che ha l'angolo piú vicino

        fill_bin in sp20lib

        leggere DBZH o TH (fare poi DBtoBYTE)
        */

    /*
    struct VOL_POL volume.scan(NEL)[NUM_AZ_X_PPI];
    T_MDB_data_header   old_data_header;

    //--------lettura volume------
    int tipo_dati_richiesti = INDEX_Z;
    int ier = read_dbp_SP20((char*)nome_file,volume.vol_pol,&old_data_header,
                            tipo_dati_richiesti,volume.nbeam_elev);

    if (ier != OK)
        LOG_ERROR("Reading %s returned error code %d", nome_file, ier);

    //  ----- Test sul volume test_file.......  --------
    if (!test_file(file_type))
    {
        LOG_ERROR("test_file failed");
        return false;
    }
    */

    // TODO: look for the equivalent of declutter_rsp and check its consistency
    // like in test_file

    //return ier == OK;
    return true;
}

void CUM_BAC::elabora_dato()
{
    const float fondo_scala = BYTEtoDB(1); // -19.7 dBZ
    HRay hray;

    //-------------leggo mappa statica ovvero first_level (funzione leggo_first_level)------------
    leggo_first_level();

    //-------------se definita qualita' leggo dem e altezza fascio (funzioni legg_dem e leggo_hray)(mi servono per calcolare qualità)
    if (do_quality)
    {
        assets.load_dem(dem);
        hray.load_hray(assets);
    }

    //------------se definito DECLUTTER , non rimuovo anap e riscrivo  volume polare facedndo declutter solo con mappa statica.... ancora valido?

    if (do_declutter)
    {
        for(unsigned i=0; i<NUM_AZ_X_PPI; i++)
        {
            for(unsigned k=0; k<volume.scan(0).beam_size; k++)
            {
                //---assegno el_inf a mappa statica
                unsigned el_inf = first_level_static(i, k);
                //---ricopio valori a mappa statica sotto
                for(unsigned l=0; l<=el_inf; l++)
                {
                    // Enrico: cerca di non leggere/scrivere fuori dal volume effettivo
                    if (k >= volume.scan(l).beam_size) continue;
                    if (k < volume.scan(el_inf).beam_size)
                        volume.scan(l).set(i, k, volume.scan(el_inf).get(i, k));
                    else
                        volume.scan(l).set(i, k, MISSING_DB);
                    //------------se definito BEAM BLOCKING e non definito BLOCNOCORR (OPZIONE PER non correggere il beam blocking a livello di mappa statica PUR SAPENDO QUANT'È)
                    if (do_beamblocking && do_bloccorr)
                    {
                        volume.scan(l).set(i, k, BeamBlockingCorrection(volume.scan(l).get(i, k), beam_blocking(i, k)));
                        //volume.scan(l).set_raw(i, k, DBtoBYTE(BYTEtoDB(volume.scan(l).get_raw(i, k))-10*log10(1.-(float)beam_blocking[i][k]/100:tab.)));
                        //volume.scan(l).set_raw(i, k, volume.scan(l).get_raw(i, k)+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5));
                    }

                }
                if (do_quality)
                    elev_fin[i][k]=el_inf;
            }
        }
        return;
    }

    //------------se non definito DECLUTTER inizio rimozione propagazione anomala al livello mappa dinamica e elaborazioni accessorie


    /* 26-5-2004 : se sono alla 1 o successive elevazioni
       e range > 60 km cambio le soglie, in modo
       da evitare di riconoscere come anaprop una pioggia shallow
       Il criterio diventa: - se la differenza tra Z all'elevazione più bassa della
       corrente e la Z corrente è <10 dbZ allora
       rendo inefficaci i limiti di riconoscimento anaprop. */

    //--------ciclo sugli azimut e bins per trovare punti con propagazione anomala----------------

    for(unsigned i=0; i<NUM_AZ_X_PPI; i++)
    {
        bool flag_anap = false;
        unsigned cont_anap=0;// aggiunto per risolvere problema di uso con preci shallow
        for(unsigned k=0; k<volume.scan(0).beam_size; k++)
        {
            //------------- incremento statistica tot ------------------
            grid_stats.incr_tot(i, k);
            // ------------assegno l'elevazione el_inf a first_level e elev_fin a el_inf---------
            int loc_el_inf = first_level(i, k);
            while ( k >= volume.scan(loc_el_inf).beam_size)
            {
                LOG_INFO("Decremento el_inf per k fuori range (i,k,beam_size,el_inf_dec) (%d,%d,%d,%d)",i,k,volume.scan(loc_el_inf).beam_size,loc_el_inf-1);
                loc_el_inf--;
            }
            if (loc_el_inf < 0) throw std::runtime_error("loc_el_inf < 0");
            //const int el_inf = first_level(i, k);
            const unsigned el_inf = loc_el_inf;

            if (do_quality)
                elev_fin[i][k]=el_inf;

            // ------------assegno a el_up il successivo di el_inf e se >=NEL metto bin_high=fondo_scala
            const unsigned el_up = el_inf +1;

            // ------------assegno bin_low e bin_high 
   
            const float bin_low  = volume.scan(el_inf).get(i, k);
	    float bin_high;
	    if (el_up >= volume.size() ) {
	        bin_high=BYTEtoDB(0);
	    } else {               
	        if ( k >= volume.scan(el_up).beam_size){
	            bin_high=BYTEtoDB(0);
                } else{
                    bin_high = volume.scan(el_up).get(i, k);
	        }
	    }
            // ------------assegno  bin_low_low (cioè il valore sotto il bin base)
           const float bin_low_low = (el_inf > 0)
                ? volume.scan(el_inf-1).get(i, k)
                : fondo_scala+1;

            //------------assegno le soglie per anaprop : se sono oltre 60 km e se la differenza tra il bin sotto il base e quello sopra <10 non applico test (cambio i limiti per renderli inefficaci)
            /* differenza massima tra le due elevazioni successive perchè non sia clutter e valore minimo a quella superiore pe il primo e per i successivi (NEXT) bins*/

            int MAX_DIF=MAX_DIF_OR;
            int MAX_DIF_NEXT=MAX_DIF_NEXT_OR;
            int MIN_VALUE=MIN_VALUE_OR;
            int MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;
            //----------questo serviva per evitare di tagliare la precipitazione shallow ma si dovrebbe trovare un metodo migliore p.es. v. prove su soglia
            if((el_inf>=1)&&(k>LIMITE_ANAP)&&(bin_low_low-bin_low<10)) //-----------ANNULLO EFFETTO TEST ANAP
            {
                // FIXME: perché BYTEtoDB se assegnamo a degli interi? [Enrico]
                MAX_DIF_NEXT=BYTEtoDB(255);
                MAX_DIF=BYTEtoDB(255);
                MIN_VALUE=BYTEtoDB(0);
                MIN_VALUE_NEXT= BYTEtoDB(0);
            }

            // ------------separo i diversi casi x analisi anaprop: ho dati sia al livello base che sopra o no  e ho trovato anaprop in precedenza sul raggio o no
            bool test_an;
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
                        flag_anap = true;
                        cont_anap=cont_anap+1;

                        //--------ricopio valore a el_up su tutte elev inferiori--------------
                        for(unsigned l=0; l<el_up; l++) {
                            volume.scan(l).set(i, k, volume.scan(el_up).get(i, k));
                        } //

                        //--------azzero beam_blocking ( ho cambiato elevazione, non ho disponible il bbeam blocking all' elev superiore)--------------
                        if (do_beamblocking)
                            beam_blocking(i, k)=0;/* beam blocking azzerato */

                        //--------------------incremento la statitica anaprop e di cambio elevazione-------------
                        grid_stats.incr_anap(i, k);
                        if (el_up > first_level_static(i, k)) grid_stats.incr_elev(i, k); //incremento la statistica cambio elevazione

                        //-------------------memorizzo dati di qualita '-------------
                        if (do_quality)
                        {
                            dato_corrotto(i, k)=ANAP_YES;/*  risultato test: propagazione anomala*/
                            elev_fin[i][k]=el_up;
                        }
                    }
                    else
                        //-----non c'è propagazione anomala:ricopio su tutte e elevazioni il valore di el_inf e correggo il beam blocking e incremento la statistica beam_blocking, assegno matrice anaprop a 0 nel punto e assegno a 0 indicatore anap nel raggio-----------
                    {
                        flag_anap = false;
                        if (do_beamblocking && do_bloccorr)
                        {
                            // FIXME: cosa dovrebbe essere l qui? Non siamo
                            // dentro a un ciclo for che itera su l [Enrico]
                            // Nella get, abbiamo messo el_inf
                            volume.scan(el_inf).set(i, k, BeamBlockingCorrection(volume.scan(el_inf).get(i, k), beam_blocking(i, k)));
                            //volume.scan(el_inf).get_raw(i, k)=DBtoBYTE(BYTEtoDB(volume.scan(l).get_raw(i, k))-10*log10(1.-(float)beam_blocking[i][k]/100.));
                            //    volume.scan(l).get_raw(i, k)=volume.scan(l).get_raw(i, k)+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5); //correggo beam blocking
                            grid_stats.incr_bloc(i, k, beam_blocking(i, k)); // incremento statistica beam blocking
                        }
// 20140128 - errore nel limite superiore ciclo
// for(l=0; l<=el_up; l++){
                        for(unsigned l=0; l<el_inf; l++){
                            volume.scan(l).set(i, k, volume.scan(el_inf).get(i, k)); // assegno a tutti i bin sotto el_inf il valore a el_inf (preci/Z a el_inf nella ZLR finale)

                        }
                        if (do_quality)
                        {
                            dato_corrotto(i, k)=ANAP_OK;/* matrice risultato test: no propagazione anomala*/
                            elev_fin[i][k]=el_inf;
                        }
                        if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k); //incremento la statistica cambio elevazione

                    }
                }
                //------------------se invece flag_anap == 0 cioè non ho ancora trovato anaprop nel raggio ripeto le operazioni ma con i limiti più restrittivi (controlla elev_fin)------------------

                else
                {
                    if(bin_low-bin_high >= MAX_DIF || bin_high <= MIN_VALUE  )

                    {
                        //--------ricopio valore a el_up su tutte elev inferiori--------------
                        for(unsigned l=0; l<el_up; l++){
           //                 volume.scan(l).set_raw(i, k, 1);
                            volume.scan(l).set(i, k, volume.scan(el_up).get(i, k));  //ALTERN
                        }
                        //---------assegno l'indicatore di presenza anap nel raggio e incremento statistica anaprop, assegno matrici che memorizzano anaprop e elevazione_finale e azzero beam blocking perchè ho cambiato elevazione
                        flag_anap = true;
                        cont_anap=cont_anap+1;
                        grid_stats.incr_anap(i, k);
                        if (do_quality)
                        {
                            dato_corrotto(i, k)=ANAP_YES;/*matrice risultato test: propagazione anomala*/
                            elev_fin[i][k]=el_up;
                        }
                        if (el_up > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                        if (do_beamblocking)
                            beam_blocking(i, k)=0;
                    }

                    //-----non c'è propagazione anomala:ricopio su tutte e elevazioni il valore di el_inf e correggo il beam blocking,  incremento la statistica beam_blocking, assegno matrice anaprop a 0 nel punto , assegno a 0 indicatore anap nel raggio, assegno elevazione finale e incremento statisica cambio elevazione se el_inf > first_level_static(i, k)-----------
                    else
                    {
                        for(unsigned l=0; l<=el_inf; l++)
                        {
                            volume.scan(l).set(i, k, volume.scan(el_inf).get(i, k));
                            if (do_beamblocking && do_bloccorr)
                            {
                                volume.scan(l).set(i, k, BeamBlockingCorrection(volume.scan(l).get(i, k), beam_blocking(i, k)));
                                //volume.scan(l).set_raw(i, k, DBtoBYTE(BYTEtoDB(volume.scan(l).get_raw(i, k))-10*log10(1.-(float)beam_blocking(i, k)/100.)));
                                //volume.scan(l).set_raw(i, k, volume.scan(l).get_raw(i, k)+ceil(-3.1875*10.*log10(1.-(float)beam_blocking(i, k)/100.)-0.5));
                                grid_stats.incr_bloc(i, k, beam_blocking(i, k));
                            }
                        }

                        if (do_quality)
                        {
                            dato_corrotto(i, k)=ANAP_OK;
                            elev_fin[i][k]=el_inf;
                        }

                        if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                        flag_anap = false;

                    } /*endif test anaprop*/
                }/*endif flaganap*/
            }/*endif bin_low > fondo_scala && bin_high >= fondo_scala*/

            //----------------se al livello base non ho dato riempio con i valori di el_up tutte le elevazioni sotto (ricostruisco il volume) e assegno beam_blocking 0
            else if (bin_low < fondo_scala)
            {
                for(unsigned l=0; l<el_up; l++)
                {
                    if (volume.scan(l).beam_size > k)
                        volume.scan(l).set(i, k, volume.scan(el_up).get(i, k));
                }
                //----------------controlli su bin_high nel caso in cui bin_low sia un no data per assegnare matrice anap  (dato_corrotto(i, k))
                if (do_quality)
                {
                    if (bin_high<fondo_scala)   dato_corrotto(i, k)=ANAP_NODAT;/*manca dato sotto e sopra*/
                    bool test_an1;
                    if (cont_anap< THR_CONT_ANAP )
                        test_an1=(bin_high>=fondo_scala); //modificato per contemplare > o >=
                    else
                        test_an1=(bin_high>fondo_scala);

                    // if (bin_high>=fondo_scala)  dato_corrotto(i, k)=ANAP_NOCONTROL;/*manca controllo (sotto non ho nulla)*/ //messo l=
                    if (test_an1) dato_corrotto(i, k)=ANAP_NOCONTROL;
                    if (bin_high==fondo_scala) dato_corrotto(i, k)=ANAP_OK;/*non piove (oppure sono sopra livello preci...)*/
                }

                if (do_beamblocking)
                    beam_blocking(i, k)=0;
            }

            //----------------se bin_low == fondo_scala riempio matrice volume.vol_pol con dato a el_inf (mi resta il dubbio di quest'if se seve o basti un else ) azzero matrice anap (dato ok)
            else if (bin_low == fondo_scala || bin_high <= fondo_scala)/* quel che resta da (bin_low > fondo_scala && bin_high >= fondo_scala) e (bin_low < fondo_scala) ; messo =per bin_high*/

            {

                for(unsigned l=0; l<el_inf; l++)//riempio con i valori di el_inf tutte le elevazioni sotto (ricostruisco il volume)
                {
                    if (volume.scan(l).beam_size > k)
                        volume.scan(l).set(i, k, volume.scan(el_inf).get(i, k));
                }

                if (do_quality)
                {
                    dato_corrotto(i, k)=ANAP_OK; // dubbio
                    elev_fin[i][k]=el_inf;
                }

                if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k);
            }
            /*-----------------------------------------------------------fine di tutti gli if-----------*/
            //-----finiti tutti i controlli assegno le varibili di qualita definitive: elevazione, quota calcolata sull'elevazione reale con propagazione standard , e quota relativa al suolo calcolata con elevazione nominale e propagazione da radiosondaggio.

            if (do_quality)
            {
                // FIXME: this reproduces the truncation we had by storing angles as short ints between 0 and 4096
                //float elevaz=(float)(volume.ray_at_elev_preci(i, k).teta_true)*CONV_RAD;
                const float elevaz = elev_fin.elevation_rad_at_elev_preci(i, k);
                // elev_fin[i][k]=first_level_static(i, k);//da togliere
                quota(i, k)=(unsigned short)(quota_f(elevaz,k));
            }
        }
    }

    LOG_INFO("elabora_dato completed");
    ScrivoStatistica();
}                         /*end funzione elabora_dato()*/

void CUM_BAC::leggo_first_level()
{
    if (do_readStaticMap)
    {
        // Leggo mappa statica
        assets.load_first_level(first_level_static);
        // Allargo per coprire la dimensione del volume
        if (first_level_static.cols() < volume.max_beam_size())
            first_level_static.resize_beams_and_propagate_last_bin(volume.max_beam_size());

        // copio mappa statica su matrice first_level
        first_level = first_level_static;
        LOG_INFO("Letta mappa statica");
    }

    if (do_beamblocking)
    {
        // Leggo file elevazioni per BB
        assets.load_first_level_bb_el(bb_first_level);
        // Allargo per coprire la dimensione del volume
        if (bb_first_level.cols() < volume.max_beam_size())
            bb_first_level.resize_beams_and_propagate_last_bin(volume.max_beam_size());

        // Leggo file valore di BB
        assets.load_first_level_bb_bloc(beam_blocking);
        // Allargo per coprire la dimensione del volume
        if (beam_blocking.cols() < volume.max_beam_size())
            beam_blocking.resize_beams_and_propagate_last_bin(volume.max_beam_size());

        /* Se elevazione clutter statico < elevazione BB, prendi elevazione BB,
           altrimeti prendi elevazione clutter statico e metti a 0 il valore di BB*/
        for(unsigned i=0; i < first_level.rows(); ++i)
            for (unsigned j=0; j < first_level.cols(); ++j)
            {
                if (do_bloccorr)
                {
                    if (first_level_static(i, j)<=bb_first_level(i, j))
                        first_level(i, j)=bb_first_level(i, j);
                    else
                    {
                        beam_blocking(i, j)=0;
                        first_level(i, j)=first_level_static(i, j);
                    }
                } else {
                    if (first_level_static(i, j)>bb_first_level(i, j))
                        beam_blocking(i, j)=0;
                    if (first_level_static(i, j)<bb_first_level(i, j))
                        beam_blocking(i, j)=OVERBLOCKING;
                }
            }
    }

    /*-------------------------------
      patch per espandere il clutter
      -------------------------------*/
    if(do_medium){
        PolarScan<unsigned char> first_level_tmp(first_level);
        LOG_INFO(" Dentro patch %u ",MyMAX_BIN);
        for (unsigned i=NUM_AZ_X_PPI; i<800; i++)
            for (unsigned j=0; j<MyMAX_BIN; j++)
                for (unsigned k=i-1; k<i+2; k++)
                    if(first_level(i%NUM_AZ_X_PPI, j) < first_level_tmp(k%NUM_AZ_X_PPI, j))
                        first_level(i%NUM_AZ_X_PPI, j)=first_level_tmp(k%NUM_AZ_X_PPI, j);
        LOG_INFO(" fine patch %u ",MyMAX_BIN);
    }
}

//------------funzione quota_f-----------------------
//---------funzione che calcola la quota in metri del centro del fascio-----------------------
//--------distanza=k*dimensionecella +semidimensionecella in metri ----------------------
//--------quota=f(distinkm, rstinkm, elevazinrad) in metri
double CUM_BAC::quota_f(double elevaz, int k) // quota funzione di elev(radianti) e range
{
    double dist = k * load_info.size_cell + load_info.size_cell / 2.;
    // quota in prop. standard da elevazione reale
    return (sqrt(pow(dist / 1000., 2) + (rst * rst) + 2.0 * dist / 1000. * rst * sin(elevaz)) - rst) * 1000.;
}

void CUM_BAC::ScrivoStatistica()
{
    //Definizioni per statistica anap
    static const int DIM1_ST = 16;
    static const int DIM2_ST = 13;
    /*--- numero minimo di celle presenti in un
      settore per la statistica            ---*/
    static const int N_MIN_BIN = 500;

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
            if (grid_stats.count(az, ran) >= N_MIN_BIN)
            {
                statistica[az][ran] = grid_stats.perc_anap(az, ran);
                statistica_bl[az][ran] = grid_stats.perc_bloc(az, ran);
                statistica_el[az][ran] = grid_stats.perc_elev(az, ran);
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

FILE *controllo_apertura (const char *nome_file, const char *content, const char *mode)
{
    LOG_CATEGORY("radar.io");
    FILE *file;
    int ier_ap=0;

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
        throw std::runtime_error("errore apertura");
    }
    if (file == NULL) {
        LOG_ERROR("Errore Apertura %s %s", content, nome_file);
        throw std::runtime_error("errore apertura");
    }
    return(file);
}

/*
   comstart caratterizzo_volume
   idx calcola qualita' volume polare
   calcola qualita' volume polare
   NB il calcolo è fatto considerando q=0 al di sotto della mappa dinamica.
   per ora drrs=dist nche nel caso di Gattatico, mentre dtrs è letto da file
   si puo' scegliere tra qualita' rispetto a Z e rispetto a R, in realtà per ora sono uguali.

   float quo:  quota bin
   float el:   elevazione
   float rst:  raggio equivalente in condizioni standard
   float dZ:   correzione vpr
   float sdevZ:  stand. dev. correzione vpr

//--- nb. non ho il valore di bb sotto_bb_first_level
comend
*/
void CUM_BAC::caratterizzo_volume()
{
    LOG_DEBUG("start caratterizzo_volume");

    HRay hray_inf; /*quota limite inferiore fascio in funzione della distanza e elevazione*/
    hray_inf.load_hray_inf(assets);

    qual = new Volume<unsigned char>(volume, 0);

    // path integrated attenuation
    double PIA;
    // dimensione verticale bin calcolata tramite approcio geo-ottico
    float dh=1.;
    // distanza radiosondaggio,
    float dhst=1.;
    // tempo dal radiosondaggio
    float drrs=1.;
    // distanza dal radar
    float dist=1.;
    // beam blocking
    unsigned char bb=0;
    // indice clutter da anaprop
    unsigned char cl=0;

    //----------ciclo su NSCAN(=6), cioè sul numero di elevazioni (nominali) per le quali ho calcolato il beam blocking
    /* a questo punto servono: bb, cl,  PIA, dtrs e drrs radiosond, quota, hsup e hinf beam-----------------*/

    //for (l=0; l<NSCAN; l++)/*ciclo elevazioni*/// NSCAN(=6) questo lascia molti dubbi sul fatto che il profilo verticale alle acquisizioni 48, 19 etc..  sia realmente con tutti i dati! DEVO SOSTITUIRE CON nel E FARE CHECK.

    for (unsigned l=0; l<volume.size(); l++)/*ciclo elevazioni*/// VERIFICARE CHE VADA TUTTO OK
    {
        const auto& scan = volume.scan(l);
        const auto& scan_info = load_info.scan(l);
        for (int i=0; i<NUM_AZ_X_PPI; i++)/*ciclo azimuth*/
        {
            //-----elevazione reale letta da file* fattore di conversione 360/4096
            // FIXME: this reproduces the truncation we had by storing angles as short ints between 0 and 4096
            //elevaz=(float)(volume.scan(l)[i].teta_true)*CONV_RAD;//--- elev reale
            //elevaz=(float)(volume.scan(l)[i].elevation*DTOR);//--- elev reale
            const double elevaz = scan_info.get_elevation_rad(i);//--- elev reale

            //--assegno PIA=0 lungo il raggio NB: il ciclo nn va cambiato in ordine di indici!
            PIA=0.;

            for (unsigned k=0; k<scan.beam_size; k++)/*ciclo range*/
            {
                double sample = volume.scan(l).get(i, k);

                //---------distanza in m dal radar (250*k+125 x il corto..)
                dist= k*load_info.size_cell+load_info.size_cell/2.;/*distanza radar */

                //-----distanza dal radiosondaggio (per GAT si finge che sia colocato ..), perchè? (verificare che serva )
                drrs=dist;
                /* if (!(strcmp(sito,"GAT")) ) {  */
                /*     drrs=dist; */
                /* } */
                /* if (!(strcmp(sito,"SPC")) ) {  */
                /*     drrs=dist; */
                /* } */


                //assegno la PIA (path integrated attenuation) nel punto e POI la incremento  (è funzione dell'attenuazione precedente e del valore nel punto)
                PIA=attenuation(DBtoBYTE(sample),PIA);

                //------calcolo il dhst ciè l'altezza dal bin in condizioni standard utilizzando la funzione quota_f e le elevazioni reali
                dhst =quota_f(elevaz+0.45*DTOR,k)-quota_f(elevaz-0.45*DTOR,k);


                //----qui si fa un po' di mischione: finchè ho il dato dal programma di beam blocking uso il dh con propagazione da radiosondaggio, alle elevazioni superiori assegno dh=dhst  e calcolo quota come se fosse prop. standard, però uso le elevazioni nominali

                if (l<hray_inf.NSCAN-1   ) {
                    dh=hray_inf[k][l+1]-hray_inf[k][l]; /* differenza tra limite sup e inf lobo centrale secondo appoccio geo-ott*/
                }
                else {
                    dh=dhst; /* non ho le altezze oltre nscan-1 pero' suppongo che a tali elevazioni la prop. si possa considerare standard*/
                    // Enrico: commentato: non viene piú letto
                    // hray[k][l]=quota_f(elevaz,k);//non lo assegno
                }

                if (l < elev_fin[i][k]) {
                    cl=ANAP_YES;
                    bb=BBMAX;
                } else if (l == elev_fin[i][k]) {
                    cl=dato_corrotto(i, k);  /*cl al livello della mappa dinamica*/
                    bb=beam_blocking(i, k);  /*bb al livello della mappa dinamica *///sarebbe da ricontrollare perchè con la copia sopra non è più così
                } else if (l > elev_fin[i][k]) {
                    cl=0;       /*per come viene scelta la mappa dinamica si suppone che al livello superiore cl=0 e bb=0*/
                    bb=0;   // sarebbe if (l-bb_first_level(i, k) >0  bb=0;  sopra all'elevazione per cui bb<soglia il bb sia =0 dato che sono contigue o più però condiz. inclusa
                }

                //------dato che non ho il valore di beam blocking sotto i livelli che ricevo in ingresso ada progrmma beam blocking e
                //--------dato che sotto elev_fin rimuovo i dati come fosse anaprop ( in realtà c'è da considerare che qui ho pure bb>50%)
                //--------------assegno qualità zero sotto il livello di elev_fin (si può discutere...), potrei usare first_level_static confrontare e in caso sia sotto porre cl=1
                if (l < elev_fin[i][k]) {
                    qual->scan(l).set(i, k, 0);
                    cl=2;
                } else {
                    //--bisogna ragionare di nuovo su definizione di qualità con clutter se si copia il dato sopra.--

                    //--calcolo la qualità--
                    // FIXME: qui tronca: meglio un round?
                    qual->scan(l).set(i, k, (unsigned char)(func_q_Z(cl,bb,dist,drrs,hray_inf.dtrs,dh,dhst,PIA)*100));
                }

                if (qual->scan(l).get(i, k) ==0) qual->scan(l).set(i, k, 1);//????a che serve???
                if (do_vpr)
                {
                    /* sezione PREPARAZIONE DATI VPR*/
                    if(cl==0 && bb<BBMAX_VPR )   /*pongo le condizioni per individuare l'area visibile per calcolo VPR, riduco il bb ammesso (BBMAX_VPR=20)*/ //riveder.....?????
                        calcolo_vpr->flag_vpr->scan(l).set(i, k, 1);
                }
            }
        }
    }

    LOG_DEBUG("End caratterizzo_volume");
    return;
}

double CUM_BAC::attenuation(unsigned char DBZbyte, double  PIA)  /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/
{
    double Zhh,att_rate,R;/* PIA diventa att_tot devo decidere infatti se PIA sarà 3d percio' temp. uso  nomi diversi*/
    double att_tot;

    //---ricevo in ingresso il dato e l'attenuazione fino  quel punto
    //---la formula recita che l'attenuazione è pari una funzione di Z reale (quindi corretta dell'attenuazione precedente). ovviamente devo avere un segnale per correggere.
    //--------- CALCOL
    att_tot=PIA;
    Zhh=(double)(BYTEtoZ(DBZbyte));
    if (10*log10(Zhh) > THRES_ATT )
    {
        Zhh=pow(10., (log10(Zhh)+ 0.1*att_tot));
        R=pow((Zhh/aMP),(1.0/bMP));
        att_rate=0.0018*pow(R,1.05);
        att_tot=att_tot+2.*att_rate*0.001*load_info.size_cell;
        if (att_tot>BYTEtoDB(254)) att_tot=BYTEtoDB(254);
    }
    return att_tot;
}

void CalcoloVPR::classifica_rain()
{
    LOG_CATEGORY("radar.class");
    int hmax=-9999;

    /* ;---------------------------------- */
    /* ;          FASE 0 :                  */
    /* ;---------------------------------- */
    // DEFINISCO QUOTE DELLA BASE E DEL TOP DELLA BRIGHT BAND USANDO IL DATO quota del picco  DEL PRECEDENTE RUN O, SE NON PRESENTE LA QUOTA DELLO ZERO DA MODELLO

    // Lettura quota massimo da VPR  calcolo base e top bright band
    LOG_INFO("data= %s",cum_bac.date);
    // calcolo il gap
    gap = cum_bac.assets.read_profile_gap();
    //-- se gap < memory leggo hmax da VPR
    if (gap<=MEMORY){
        hmax = cum_bac.assets.read_vpr_hmax();
            //---suppongo una semiampiezza massima della bright band di 600 m e definisco htopbb e hbasebb come hmassimo +600 m (che da clima ci sta) e hmassimo -600 m
    }

    if (hmax >= 0)
    {
        hbbb=(hmax-600.)/1000.;
        htbb=(hmax+600.)/1000.;
    } else {
        //-- se gap > memory o se non ho trovato il file
        // Lettura 0 termico da modello, e calcolo base e top bright band
        LOG_INFO("leggo 0termico per class da file %s",getenv("FILE_ZERO_TERMICO"));
        //  leggo informazioni di temperatura da modello*/
        float zeroterm;//zerotermico
        if (cum_bac.assets.read_0term(zeroterm))
        {
            //-- considerato che lo shift medio tra il picco e lo zero è tra 200 e 300 m, che il modello può avere un errore, definisco cautelativamente htbb come quota zero + 400 m e hbbb come  quota zero -700 m  .
            htbb=zeroterm/1000. + 0.4; // se non ho trovato il vpr allora uso un range più ristretto, potrebbe essere caso convettivo
            hbbb=zeroterm/1000. - 1.0;
        } else {
            LOG_ERROR("non ho trovato il file dello zero termico");
            LOG_INFO("attenzione, non ho trovat zero termico ne da vpr ne da radiosondaggio");
            htbb=0; // discutibile così faccio tutto con VIZ
            hbbb=0;
        }
    }

    // se hbasebb è <0 metto 0
    if (hbbb<0) hbbb=0;

    LOG_INFO("calcolati livelli sopra e sotto bright band hbbb=%f  htbb=%f",hbbb,htbb);

    // TODO: remove duplication with CylindricalVolume::resample
    const Volume<double>& volume = cum_bac.volume;
    const double size_cell = cum_bac.load_info.size_cell;
    double range_min=0.5 * size_cell/1000.;
    double range_max=(MyMAX_BIN-0.5) * size_cell/1000.;
    double xmin=floor(range_min*cos(volume.elevation_max()*DTOR)); // distanza orizzontale minima dal radar
    double zmin=pow(pow(range_min,2.)+pow(4./3*REARTH,2.)+2.*range_min*4./3.*REARTH*sin(volume.elevation_min() * DTOR),.5) -4./3.*REARTH+h_radar; // quota  minima in prop standard
    double xmax=floor(range_max*cos(volume.elevation_min()*DTOR)); // distanza orizzontale massima dal radar
    double zmax=pow(pow(range_max,2.)+pow(4./3*REARTH,2.)+2.*range_max*4./3.*REARTH*sin(volume.elevation_max() * DTOR),.5) -4./3.*REARTH+h_radar;//quota massima


    x_size=(xmax-xmin)/RES_HOR_CIL; //dimensione orizzontale
    // FIXME: usiamo volume.max_beam_size invece di MyMAX_BIN?
    if (x_size > MyMAX_BIN) x_size=MyMAX_BIN;
    z_size=(zmax-zmin)/RES_VERT_CIL; //dimensione verticale

    // ricampionamento del volume in coordinate cilindriche
    CylindricalVolume cil(NUM_AZ_X_PPI, x_size, z_size, 0);
    cil.resample(cum_bac.volume, MyMAX_BIN, cum_bac.load_info.size_cell);

    //-------------------------------------------------------------------------------------------------------------------------
    // faccio la classificazione col metodo Vertical Integrated Reflectivity
    CalcoloVIZ viz(cil, htbb, hbbb, t_ground);
    viz.classifico_VIZ();

    //classificazione con STEINER
    //  if (hmax > 2000.) {// per evitare contaminazioni della bright band, si puo' tunare
    // if (hbbb > 500.) {// per evitare contaminazioni della bright band, si puo' tunare
    CalcoloSteiner steiner(cum_bac.volume, cum_bac.elev_fin, cum_bac.MyMAX_BIN, cum_bac.load_info.size_cell);
    steiner.calcolo_background();
    steiner.classifico_STEINER();
    //  }
    merge_metodi(steiner, viz);
    return ;
}

void CalcoloVPR::merge_metodi(const CalcoloSteiner& steiner, const CalcoloVIZ& viz)
{
    //inizializzazione vettori e matrici
    for (unsigned i=0; i<NUM_AZ_X_PPI; i++){
        conv[i]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;

        for (unsigned j=0; j<x_size; j++){ // cambiato da x_size
            conv[i][j]=MISSING;
        }
    }

    for (unsigned j=0; j<NUM_AZ_X_PPI; j++)
        for (unsigned k=0; k<x_size; k++)
            if (steiner.conv_STEINER(j, k) == viz.conv_VIZ(j, k) && steiner.conv_STEINER(j, k) > 0 && viz.stratiform(j, k) < 1)
                conv[j][k] = viz.conv_VIZ(j, k);
}

//----------ALGORITMO
/*  combina il profilo verticale corrente con quello precedente tramite il metodo di Germann (2003)
    a) calcolo gap tra ultimo profilo e istante corrente
    b) se gap < MEMORY leggo profilo
    c) se gap> MEMORY peso 0 il profilo storico e cerco oltre la data (per casi vecchi)
    d) faccio func_vpr
    f) cerco il profilo con cui combinare (->proprio, se gap<MEMORY ->dell'altro radar se gap_res<MEMORY e profile_heating_res=WARM)
    g) Combino livelli con peso sottostante
    Dati cv e ct, volume totale e volume precipitante il peso del vpr istantaneo è calcolato come segue:
    c0=2*cv;
    peso=(float)ct/(c0+ct)
    long int c0,cv,ct; costanti di combinazione (v. ref.)
    h) trovo livello minimo, se livello minimo profilo combinato più alto del precedente calcolo la diff media e sommo al vecchio
    e) ricalcolo livello minimo
    float vpr0[NMAXLAYER],vpr1[NMAXLAYER],vpr[NMAXLAYER]; profilo precedente, ultimo e combinato
    float alfat,noval; peso, nodata
    FILE *file;
    int mode,ilay;  modalità calcolo profilo (0=combinazione, 1=istantaneo),indice di strato
*/
int CalcoloVPR::combina_profili()
{
    LOG_CATEGORY("radar.vpr");
    long int c0, cv = 0, ct = 0;
    vector<float> vpr0(NMAXLAYER, NODATAVPR);
    vector<float> vpr1(NMAXLAYER, NODATAVPR);
    float vpr_dbz;
    float alfat,noval;
    int mode,ilay,i,foundlivmin=0,il,ier_ap,combinante=0; // combinante: variabile che contiene presenza vpr alternativo
    vector<long int> area(NMAXLAYER, NODATAVPR);
    long int ar=0;
    int n=0,diff=0;
    char nomefile[150],stringa[100];
    struct tm *T_tempo;
    time_t Time,T_Time;
    FILE *file;

    mode=MOD_VPR;
    noval=NODATAVPR;


    /* questo per fare ciclo sul vpr vecchio*/
    Time = cum_bac.NormalizzoData(cum_bac.load_info.acq_date);

    //--------inizializzo cv e ct-------------//
    //-----calcolo del profilo istantaneo:faccio func_vpr-----//

    ier_vpr=func_vpr(&cv,&ct,vpr1,area_vpr); // ho fatto func_vpr, il profilo istantaneo
    LOG_INFO("fatta func vpr %d", ier_vpr);


    /*modalità VPR combinato*/

    if(mode == 0) {

        /*----calcolo il peso c0 per la combinazione dei profili*/

        c0=2*(cv);

        /*------calcolo la distanza temporale che separa l'ultimo profilo calcolato dall'istante attuale--*/
        /* (dentro il file LAST_VPR c'è una data che contiene la data cui si riferisce il vpr in n0 di secondi dall'istante di riferimento)*/

        gap = cum_bac.assets.read_profile_gap();


        /*------leggo il profilo vecchio più recente di MEMORY ----*/
        /*------nota bene: è in R ovvero  pioggia!! ----*/

        if (!cum_bac.assets.read_vpr0(vpr0, area))
        {
            LOG_WARN("non esiste file vpr vecchio: %s",getenv("VPR0_FILE"));

            //----se file non esiste assegno gap=100----
            gap=100;
        }

        //------------se gap < MEMORY leggo vpr e area per ogni strato-----------
        //--------qui dentro c'è la funzione controllo_apertura, per la quale rimandiamo a dopo qualsiasi commento--------

        if (gap <= MEMORY)
        {
            combinante=1;
        } else {
            //-----Se gap > MEMORY

            //a)----- tento .. sono in POST-ELABORAZIONE:----

            //-----devo andare a ricercare tra i profili 'buoni' in archivio quello con cui combinare il dato----
            //---- trattandosi di profili con data nel nome del file, costruisco il nome a partire dall'istante corrente ciclando su un numero di quarti d'ora
            //---- pari a memory finchè non trovo un profilo. se non lo trovo gap resta=100

            for (i=0;i<MEMORY;i++){

                //---calcolo della data---//

                T_Time=Time+i*900;
                T_tempo=gmtime(&T_Time);

                sprintf(nomefile,"%s/%04d%02d%02d%02d%02d_vpr_%s",getenv("DIR_STORE_VPR"), //--questa non sarebbe una dir_arch più che store?... contesto il nome...
                        T_tempo->tm_year+1900, T_tempo->tm_mon+1, T_tempo->tm_mday,
                        T_tempo->tm_hour, T_tempo->tm_min,getenv("SITO"));
                ier_ap=access(nomefile,R_OK);

                //---- se non ho errore apertura metto gap=0 e metto profilo 'caldo', leggo il profilo e lo converto in R e interrompo il ciclo di ricerca!
                if  (!ier_ap){
                    file=fopen(nomefile,"r");
                    gap=0;
                    heating=WARM;
                    fscanf(file," %s %s %s %s" ,stringa ,stringa,stringa,stringa);
                    for (ilay=0;  ilay<NMAXLAYER; ilay++){
                        fscanf(file," %i %f %li", &il, &vpr0[ilay], &ar);  //---NB il file in archivio è in dBZ e contiene anche la quota----

                        //---- converto in R il profilo vecchio--
                        if (vpr0[ilay]>0){
                            vpr_dbz=vpr0[ilay];
                            vpr0[ilay] = DBZtoR(vpr_dbz,cum_bac.aMP,cum_bac.bMP);
                            area[ilay]=ar;
                        }
                        else
                            vpr0[ilay] = NODATAVPR;
                    }
                    combinante=1;
                    break;
                }
            }
        }

        //----a fine calcolo sul sito in esame stampo il valore del gap
        LOG_INFO("gap %li",gap);

        //TOLTA: combinazione dell'istantaneo col vecchio dell'altro radar purchè sia 'caldo' (non prevista la post-combinazione)


        //-----se è andata male la ricerca dell'altro e anche il calcolo dell'istantaneo esco

        if ( !combinante && ier_vpr)
            return (1);


        //----------------se invece l'istantaneo c'è o ho trovato un file con cui combinare

        //-----se ho i due profili riempio parte bassa con differenza media  allineandoli e combino poi

        if (!ier_vpr && combinante) {
            // calcolo la diff media
            diff=0;
            for (ilay=0;  ilay<NMAXLAYER; ilay++){
                if ( vpr0[ilay]> NODATAVPR && vpr1[ilay]>NODATAVPR ){
                    diff=diff + vpr0[ilay]-vpr1[ilay];
                    n=n+1;
                }
            }
            if (n>0){
                diff=diff/n;
                for (ilay=0; ilay<livmin/TCK_VPR; ilay++){
                    if (vpr0[ilay]<= NODATAVPR && vpr1[ilay] > NODATAVPR)
                        vpr0[ilay]=vpr1[ilay]-diff;
                    if (vpr1[ilay]<= NODATAVPR && vpr0[ilay] > NODATAVPR)
                        vpr1[ilay]=vpr0[ilay]+diff;

                }
            }
            // peso vpr corrente per combinazione
            alfat=(float)ct/(c0+ct);
            for (ilay=0;  ilay<NMAXLAYER; ilay++){
                if (vpr0[ilay] > NODATAVPR && vpr1[ilay] > NODATAVPR)
                    vpr[ilay]=comp_levels(vpr0[ilay],vpr1[ilay],noval,alfat);// combino livelli
            }
        }



        else { // se il calcolo dell'istantaneo non è andato bene , ricopio l'altro vpr e la sua area
            if (combinante){
                for (ilay=0;  ilay<NMAXLAYER; ilay++){
                    area_vpr[ilay]=area[ilay];
                    vpr[ilay]=vpr0[ilay];
                }
            }
            else{
                // se il calcolo dell'istantaneo  è andato bene ricopio il profilo
                for (ilay=0; ilay<NMAXLAYER; ilay++) vpr[ilay]=vpr1[ilay];
            }
        }
    }


    /*fine mode=0 VPR combinato, mode=1 VPR istantaneo controllo se l'istantaneo è andato ok e in caso affermativo continuo*/

    else  {
        if (ier_vpr)
            return (1);
        for (ilay=0; ilay<NMAXLAYER; ilay++) vpr[ilay]=vpr1[ilay];
    }


    //------------- trovo livello minimo -------

    livmin=0;
    foundlivmin=0;
    for (ilay=0; ilay<NMAXLAYER; ilay++){
        if (vpr[ilay]> NODATAVPR && !foundlivmin) {
            livmin=ilay*TCK_VPR+TCK_VPR/2;
            foundlivmin=1;
        }
    }
    LOG_INFO(" livmin %i", livmin);

    if (livmin>=(NMAXLAYER-1)*TCK_VPR+TCK_VPR/2  || !foundlivmin) return (1);


    //-----scrivo il profilo e la sua area-----
    cum_bac.assets.write_vpr0(vpr, area_vpr);

    return(0);
}

int CalcoloVPR::profile_heating()
#include <vpr_par.h>
{
    LOG_CATEGORY("radar.vpr");
    //---leggo ultimo file contenente riscaldamento , se non esiste impongo heating=0 (verificare comando)
    int heating = cum_bac.assets.read_vpr_heating();

    //--una volta letto il file, se il calcolo del vpr è andato bene incremento di uno heating sottraendo però la differenza di date (in quarti d'ora)-1 tra gli ultimi due profili
    //--lo faccio perchè potrei avere heating più alto del dovuto se ho avuto un interruzione del flusso dei dati
    //-- heating ha un valore massimo pari  WARM dopodichè diventa heating = MEMORY e così resta finchè non sono passati MEMORY istanti non aggiornati ( va bene?)
    //---se il profil non  è stato aggiornato invece decremento la variabile riscaldamento di gap con un minimo pari a 0
    if (ier_vpr){
        heating=heating-gap; /*se il profilo non è stato aggiornato, ho raffreddamento, in caso arrivi sotto WARM riparto da 0, cioè serve riscaldamento  */
    }
    else  {
        heating=heating-gap+2; /*se il profilo è stato aggiornato, ho riscaldamento , in caso arrivi sopra WARM riparto da MEMORY  */
        if (heating>=WARM) heating=MEMORY;  /* se heating raggiunge WARM allora lo pongo uguale a MEMORY     */

        cum_bac.assets.write_last_vpr();
    }
    if (heating<0) heating=0;

    //----stampo heating su file
    cum_bac.assets.write_vpr_heating(heating);

    //----stampo log vpr
    LOG_INFO("gap %li heating %i",gap,heating);

    return(heating);
}


int CalcoloVPR::stampa_vpr()
{
    float vpr_dbz;
    int ilay;
    FILE *file;

    file=controllo_apertura(getenv("VPR_ARCH")," ultimo vpr in dBZ per il plot","w");
    fprintf(file," QUOTA   DBZ    AREA PRECI(KM^2/1000)\n" );
    for (ilay=0;  ilay<NMAXLAYER; ilay++){
        if (vpr[ilay]> 0.001 ) {
            vpr_dbz=cum_bac.RtoDBZ(vpr[ilay]);
            fprintf(file," %i %10.3f %li\n", ilay*TCK_VPR+TCK_VPR/2, vpr_dbz, area_vpr[ilay]);
        }
        else
            fprintf(file," %i %10.3f %li\n", ilay*TCK_VPR+TCK_VPR/2, NODATAVPR, area_vpr[ilay]);
    }
    fclose(file);
    return 0;
}
/*=======================================================================================*/
/*
   comstart corr_vpr
   idx corregge i dati tramite profilo verticale
   corregge i dati tramite profilo verticale a partire da quelli con valore maggiore di THR_CORR(v vpr_par.h): riporto il dato alla quota del livello 'liquido' tramite "traslazione" cioè aggiungo al valore in dBZ la differenza tra il valore del VPR alla quota 'liquida' e quello alla quota della misura, anche in caso di neve,purchè esista il livello liquido nel profilo. In questo caso pero' flaggo il bin tramte la variabile neve[][]. In caso il profilo intero sia di neve allora riporto al valore di Z al suolo (o al livello rappresentativo) perchè non ho un valore di riferimento 'liquido'.
   La correzione avviene solo se heating>warm

   int ilref : livello del suolo o della quota rappresentativa di esso nel VPR (dove considero buoni i dati a partire dati da REPR_LEV)
   int ilray : livello a cui ho il dato
   int ilay2 : livello sopra o sotto ilray a seconda che il fascio stia sopra o sotto il centro di ilray, serve per interpolare il dato vpr su retta ed evitare correzioni a gradino
   int heating,warm: indicano quanto è caldo il profilo, e la soglia di riscaldamento
   int snow : indica se il profilo è di neve (snow=1)
   int neve[NUM_AZ_X_PPI][MAX_BIN]: 1 se c'è neve, identificata se quota dem> hvprmax+300m
   float corr: correzione in dB
   float vpr_liq: valore di R nel VPR alla quota 'liquida' ricavato dalla funzione analyse_VPR

   ilref= (dem[i][k]>REPR_LEV)?(floor(dem[i][k]/TCK_VPR)):( floor(REPR_LEV/TCK_VPR)); in pratica livello dem se > REPR_LEV, livello di REPR_LEV altrimenti.
   ilray=floor((hbin[i][k])/TCK_VPR)
   corr=RtoDBZ(vpr_liq)-RtoDBZ(vpr_hray)
   comend
*/
int CalcoloVPR::corr_vpr()
    //* ====correzione profilo====================================*/

#include <vpr_par.h>

{
    LOG_CATEGORY("radar.vpr");

    int ilray,ilref,ilay2,ier_ana,snow,strat;
    float corr,vpr_liq,vpr_hray,hbin,hliq;

    /*inizializzazione variabili */
    snow=0;
    //vpr al livello liquido liv liquido e liv max
    vpr_liq=NODATAVPR;
    hliq=NODATAVPR;
    hvprmax=INODATA;

    // analisi vpr

    ier_max=trovo_hvprmax(&hvprmax);
    ier_ana=analyse_VPR(&vpr_liq,&snow,&hliq);
    LOG_INFO("ier_analisi %i",ier_ana) ;

    /* se analisi dice che non è il caso di correggere non correggo (NB in questo caso non riempio la matrice di neve)*/
    if (ier_ana) return 1;

    LOG_INFO("altezza bright band %i",hvprmax);
    LOG_INFO("CORREGGO VPR");


    //correzione vpr
    for (unsigned i=0; i<NUM_AZ_X_PPI; i++){
        for (unsigned k=0; k<cum_bac.volume.scan(0).beam_size; k++){
            corr=0.;
            /* trovo elevazione reale e quota bin*/
            //elevaz=(float)(volume_at_elev_preci(i, k).teta_true)*CONV_RAD;
            hbin=(float)cum_bac.quota(i, k);

            /* se dall'analisi risulta che nevica assegno neve ovunque*/
            if (snow) neve(i, k)=1;
            strat=1;
            if (cum_bac.do_class)
            {
                if (conv[i][k] >= CONV_VAL){
                    strat=0;
                }
            }
            //--- impongo una soglia per la correzione pari a 0 dBZ

            if (cum_bac.volume.scan(0).get(i, k) > THR_CORR && hbin > hliq  && strat){

                //---trovo lo strato del pixel, se maggiore o uguale a NMAXLAYER lo retrocedo di 2, se minore di livmn lo pongo uguale a livmin
                ilray=(hbin>=livmin)?(floor(hbin/TCK_VPR)):(floor(livmin/TCK_VPR));//discutibile :livello del fascio se minore di livmin posto=livmin
                if (ilray>= NMAXLAYER ) ilray=NMAXLAYER-2;//livello del fascio se >= NMAXLAYER posto =NMAXLAYER-2

                //---trovo ilay2 strato con cui mediare per calcolare il vpr a una quota intermedia tra 2 livelli, se l'altezza del bin è sopra metà strato prendo quello sopra altrimenti quello sotto
                if ((int)hbin%TCK_VPR > TCK_VPR/2) ilay2=ilray+1;
                else ilay2=(int)(fabs(ilray-1));
                if (ilay2< livmin) ilay2=livmin;

                //trovo ilref: livello di riferimento per ricostruire il valore vpr al suolo nel caso di neve.
                // in caso di profilo di pioggia mi riporto sempre al valore del livello liquido e questo può essere un punto critico.. vedere come modificarlo.

                ilref=(cum_bac.dem(i, k)>livmin)?(floor(cum_bac.dem(i, k)/TCK_VPR)):(floor(livmin/TCK_VPR));//livello di riferimento; se livello dem>livmin = livello dem altrimenti livmin


                if (vpr[ilref] > 0 && vpr[ilray] > 0 ){ /*devo avere dati validi nel VPR alle quote considerate!*/
                    //-- calcolo il valore del profilo alla quota di interesse
                    vpr_hray=vpr[ilray]+((vpr[ilray]-vpr[ilay2])/(ilray*TCK_VPR-TCK_VPR/2-ilay2*TCK_VPR))*(hbin-ilray*TCK_VPR-TCK_VPR/2); /*per rendere la correzione continua non a gradini */
                    //--identifico le aree dove nevica stando alla quota teorica dello zero termico

                    if (cum_bac.dem(i, k)> hvprmax+HALF_BB-TCK_VPR || snow){ /*classifico neve*/
                        neve(i, k)=1;

                    }

                    //--se nevica la correzione consiste solo nel riportare il valore del vpr al suolo: PROPOSTA: qui si potrebbe generare una mappa di intensità di neve ma deve essere rivisto tutto


                    //if(snow) //A rimosso, faccio una cosa diversa
                    if(neve(i, k)){

                        //faccio la regressione lineare dei punti del profilo sopra il punto del dem
                        //calcolo il valore al livello del dem e lo sostituisco a vpr[ilref] nella correzione
                        // faccio linearizzazione in maniera becera:
                        //vpr[ilref]=(vpr[ilref+7]-vpr[ilref+2])/(5)*(ilref-(ilref+2))+vpr[ilref+2];

                        //passaggio=BYTEtoR(volume.vol_pol,aMP_SNOW,bMP_SNOW)

                        //volpol[0][i][k]=RtoBYTE(passaggio)

                        corr=cum_bac.RtoDBZ(vpr[ilref])-cum_bac.RtoDBZ(vpr_hray);

                        cum_bac.volume.scan(0).set(i, k, RtoDBZ(
                                    BYTE_to_mp_func(
                                        DBtoBYTE(cum_bac.volume.scan(0).get(i, k)),
                                        aMP_SNOW,
                                        bMP_SNOW),
                                    aMP_class,
                                    bMP_class));

                    }
                    else
                        // -- altrimenti correggo comunque a livello liquido :
                        corr=RtoDBZ(vpr_liq,aMP_class,bMP_class)-RtoDBZ(vpr_hray,aMP_class,bMP_class);/*riporto comunque al valore liquido anche se sono sopra la bright band*/

                    // --  controllo qualità su valore correzione
                    if (corr>MAX_CORR) corr=MAX_CORR; /*soglia sulla massima correzione*/
                    if (hbin<hvprmax && corr>0.) corr=0; /*evito effetti incrementi non giustificati*/

                    //controllo qualità su valore corretto e correzione
                    double corrected = cum_bac.volume.scan(0).get(i, k) + corr;
                    if (corrected > MAXVAL_DB) // se dato corretto va fuori scala assegno valore massimo
                        cum_bac.volume.scan(0).set(i, k, MAXVAL_DB);
                    else if ( corrected < MINVAL_DB) // se dato corretto va a fodoscala assegno valore di fondo scala
                        cum_bac.volume.scan(0).set(i, k, MINVAL_DB);
                    else
                        cum_bac.volume.scan(0).set(i, k, corrected);  // correggo

                    corr_polar(i, k)=(unsigned char)(corr)+128;


                    //inserisco un ponghino per rifare la neve con aMP e bMP modificati // DA SCOMMENTARE SE DECIDO DI FARLO

                    //if (neve[i][k]) volume.scan(0).get_raw(i, k)=DBtoBYTE(RtoDBZ( BYTE_to_mp_func(volume.scan(0).get_raw(i, k),aMP_SNOW,bMP_SNOW),aMP_class,bMP_class )) ;


                }
            }
        }
    }
    return(0);
}

int CalcoloVPR::trovo_hvprmax(int *hmax)
{
    int i,imax,istart,foundlivmax;
    float h0start,peak,soglia;


    if (t_ground != NODATAVPR)
    {
        LOG_DEBUG("trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca");
        h0start=t_ground/9.8*1000 ;
        istart=h0start/TCK_VPR -2;
        if (istart< livmin/TCK_VPR) istart=livmin/TCK_VPR;
        LOG_DEBUG("t_ground h0start istart %f %f %i",t_ground,h0start,istart);
    }
    else {
        LOG_DEBUG("trovo hvprmax  a partire da livmin");
        istart=livmin/TCK_VPR+1;
    }


    /* trovo hvprmax e il suo livello a partire dal livello istart */

    //--inizializzazione
    foundlivmax=0;
    peak=NODATAVPR;
    *hmax=INODATA;
    // Enrico vprmax=NODATAVPR;
    imax=INODATA;
    soglia=DBZtoR(THR_VPR,200,1.6); // CAMBIATO, ERRORE, PRIMA ERA RtoDBZ!!!!VERIFICARE CHE IL NUMERO PARAMETRI FUNZIONE SIA CORRETTO

    //--se vpr al livello corrente e 4 layer sopra> soglia, calcolo picco
    if (vpr[istart] >soglia && vpr[istart+4] > soglia){
        peak=10*log10(vpr[istart]/vpr[istart+4]);//inizializzo il picco
        LOG_DEBUG("peak1 = %f",peak);
    }
    //----se picco > MINIMO il punto è ok
    if(peak> MIN_PEAK_VPR){
        imax=istart;
        // Enrico vprmax=vpr[imax];
        LOG_DEBUG("il primo punto soddisfa le condizioni di picco");
    }
    for (i=istart+1;i<NMAXLAYER-4;i++) //la ricerca è un po' diversa dall'originale.. trovo il picco + alto con valore  rispetto a 4 sopra > soglia
    {
        if (vpr[i] <soglia || vpr[i+4] < soglia) break;
        peak=10*log10(vpr[i]/vpr[i+4]);
        if (vpr[i]>vpr[i-1]  && peak> MIN_PEAK_VPR ) // se vpr(i) maggiore del massimo e picco sufficientemente alto
        {
            imax=i;
            // Enrico vprmax=vpr[imax];
        }

    }

    if ( imax  > INODATA ){
        foundlivmax=1;
        peak=10*log10(vpr[imax]/vpr[imax+4]);
        *hmax=imax*TCK_VPR+TCK_VPR/2;
        LOG_DEBUG("trovato ilaymax %i %i",*hmax,imax);
        LOG_DEBUG(" picco in dbR %f",peak);
    }

    LOG_DEBUG("exit status trovo_hvprmax %i",foundlivmax);
    return (foundlivmax);
}

/*
   1)hvprmax - semiampiezza BB > liv 100 m => la bright band sta sopra il suolo => interpolo il profilo per trovare il livello liquido
   se interpolazione fallisce NON CORREGGO (scelta cautelativa, potrei avere caso convettivo
   o pochi punti o molto disomogeneo)
   2)hvprmax - semiampiezza BB < liv 100 m  => la bright contiene o è sotto il livello 100 m  oppure ho un massimo 'spurio'
   quindi calcolo la Tground, temperatura media nelle stazioni più vicine al radar, se non trovo T torno al caso 1.
   2A) Tground>T_MAX_ML ----> prob. caso convettivo o max preci passa sopra radar, interpolo il profilo e calcolo livello liquido
   se interpolazione fallisce NON CORREGGO
   2B) T_MIN_ML<T0<T_MAX_ML : prob. bright band al suolo, interpolo il profilo per trovare il livello liquido
   se interpolazione fallisce NON CORREGGO
   2C) Tground<T_MIN_ML ----> prob. profilo neve NON interpolo e non calcolo vpr al livello liquido perchè non ho livello liquido

   comend
*/
int CalcoloVPR::analyse_VPR(float *vpr_liq,int *snow,float *hliq)
    /*=======analisi profilo============ */
{
    int ier=1,ier_ana=0,liv0;
    char date[]="000000000000";
    struct tm *tempo;
    time_t Time;

    // ------------inizializzazione delle variabili ----------

    //strcpy(date,"000000000000");

    int tipo_profilo=-1;
    float v600sottobb=NODATAVPR;
    float v1000=NODATAVPR;
    float v1500=NODATAVPR;
    float vliq=NODATAVPR;
    float vhliquid=NODATAVPR;
    float vprmax=NODATAVPR;
    //*togliere gli ultimi tre*/;

    //ier_max=trovo_hvprmax(&hvprmax);


    if (t_ground == NODATAVPR) //1  se non ho nè T nè il massimo esco altrimenti tipo_profilo=0
    {
        LOG_WARN("non ho T,...");

        if ( ! ier_max ) {
            LOG_ERROR(" non ho trovato hvprmax, nè T, esco");
            return 1;
        }
        tipo_profilo=0;
    }
    else
    {

        if (t_ground >= T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100.){ //1  se T > T_MAX_ML e non ho il massimo esco
            if ( ! ier_max ) {
                LOG_ERROR(" temperatura alta e non ho trovato hvprmax, esco");
                return 1;
            }
            tipo_profilo=0;
        }


        // if (t_ground >= T_MAX_SN+0.65*(float)(livmin+TCK_VPR/2)/100  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100. )
        if (t_ground >= T_MAX_SN  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100. )
        {

            if (  ier_max ) {
                LOG_INFO(" temperatura da scioglimento e massimo in quota");
                tipo_profilo=2;
            }
            else{
                LOG_ERROR(" temperatura da scioglimento ma superiore a temperatura max neve e non ho trovato hvprmax, esco");
                return 1;
            }
            // solo una scritta per descrivere cos'è accaduto
            liv0=livmin+HALF_BB;
            if (hvprmax > liv0) LOG_INFO(" il livello %i è sotto la Bright band, ma T bassa  interpolo",livmin);
            else LOG_INFO(" il livello %i potrebbe essere dentro la Bright Band, interpolo",livmin);

        }

        //if (t_ground >= T_MIN_ML  && t_ground < T_MAX_SN+0.65*(float)(livmin+TCK_VPR/2)/100.)
        if (t_ground < T_MAX_SN)
        {
            if ( ier_max ){
                LOG_INFO(" temperatura da neve o scioglimento e massimo in quota");
                tipo_profilo=2;
            }
            else {
                LOG_INFO(" temperatura da neve o scioglimento e massimo non trovato,neve , non interpolo");
                tipo_profilo=3;
                hvprmax=0;
            }
        }

    }

    // InterpolaVPR_NR iv;
    InterpolaVPR_GSL iv;

    switch
        (tipo_profilo)
        {
            case 0:
            case 1:
            case 2:
                ier=iv.interpola_VPR(vpr.data(), hvprmax, livmin);
                if (ier){
                    LOG_INFO(" interpolazione fallita");
                    switch (tipo_profilo)
                    {

                        /*Questo fallisce se la bright band non è al suolo (per evitare correzioni dannose in casi poco omogenei)*/
                        case 0:
                        case 1:
                            ier_ana=1;
                            *vpr_liq=NODATAVPR;
                            *hliq=NODATAVPR;
                            break;
                        case 2:
                            *vpr_liq=vpr[(hvprmax+1000)/TCK_VPR]*2.15;/*21 aprile 2008*/
                            *hliq=0;
                            break;
                    }
                }
                else{
                    LOG_INFO(" interpolazione eseguita con successo");
                    //
                    // stampa del profilo interpolato
                    char file_vprint[200];
                    sprintf(file_vprint,"%s_int",getenv("VPR_ARCH"));
                    FILE* file=controllo_apertura(file_vprint," vpr interpolato ","w");
                    for (unsigned i = 0; i < NMAXLAYER; ++i)
                        fprintf(file," %f \n", cum_bac.RtoDBZ(iv.vpr_int[i]));
                    fclose(file);

                    /*calcolo valore di riferimento di vpr_liq per l'acqua liquida nell'ipotesi che a[2]=quota_bright_band e a[2]-1.5*a[3]=quota acqua liquida*/
                    if (tipo_profilo == 2 ) {
                        *hliq=(iv.E-2.1*iv.G)*1000.;
                        //lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);
                        if (*hliq<0)
                            *hliq=0;  /*con casi di bright band bassa.. cerco di correggere il più possibile*/
                        *vpr_liq=vpr[(hvprmax+1000)/TCK_VPR]*2.15;
                    }
                    else {
                        *hliq=(iv.E-2.1*iv.G)*1000.;
                        //lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);
                        if ( *hliq > livmin) {
                            *vpr_liq=vpr[(int)(*hliq/TCK_VPR)]; // ... SE HO IL VALORE VPR USO QUELLO.
                        }
                        else // altrimenti tengo il valore vpr neve + 6 dB* e metto tipo_profilo=2
                        {
                            if (*hliq<0) *hliq=0;
                            tipo_profilo=2;
                            *vpr_liq=vpr[(hvprmax+1000)/TCK_VPR]*2.15;
                        }
                    }
                }
                break;
            case 3:
                *snow=1;
                *vpr_liq=NODATAVPR;
                *hliq=NODATAVPR;
                break;
        }
    LOG_INFO("TIPO_PROFILO= %i vpr_liq %f hliq %f", tipo_profilo, *vpr_liq,*hliq );


    /* parte di stampa test vpr*/

    /* nome data */
    //definisco stringa data in modo predefinito
    Time = cum_bac.NormalizzoData(cum_bac.load_info.acq_date);
    tempo = gmtime(&Time);
    sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
            tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
    if (! ier ) {
        if(*hliq > livmin +200 )
            vhliquid=cum_bac.RtoDBZ(vpr[(int)(*hliq)/TCK_VPR]);
        vliq=cum_bac.RtoDBZ(*vpr_liq);
    }
    if (ier_max) {
        if ( hvprmax-600 >= livmin )
            v600sottobb=cum_bac.RtoDBZ(vpr[(hvprmax-600)/TCK_VPR]);
        if ((hvprmax+1000)/TCK_VPR < NMAXLAYER )
            v1000=cum_bac.RtoDBZ(vpr[(hvprmax+1000)/TCK_VPR]);
        if ((hvprmax+1500)/TCK_VPR < NMAXLAYER )
            v1500=cum_bac.RtoDBZ(vpr[(hvprmax+1500)/TCK_VPR]);
        vprmax=cum_bac.RtoDBZ(vpr[(hvprmax/TCK_VPR)]);
    }

    fprintf(test_vpr,"%s %i %i %f %f %f  %f %f %f %f %f %f %f %f  %f %f %f  %f \n",date,hvprmax,tipo_profilo,stdev,iv.chisqfin,*hliq,vliq,vhliquid,v600sottobb,v1000+6,v1500+6,vprmax,iv.rmsefin,iv.B,iv.E,iv.G,iv.C,iv.F);
    fclose(test_vpr);

    // fine parte di stampa test vpr

    //---SCRIVO ALTEZZA MASSIMO PER CLASSIFICAZIONE AL RUN SUCCESSIVO

    cum_bac.assets.write_vpr_hmax(hvprmax);
    LOG_INFO("fatta scrittura hmax vpr = %d",hvprmax);

    return (ier_ana);
}

/*
idx funzione calcolo VPR istantaneo
calcola il VPR istantaneo secondo il metodo di Germann e Joss (2003)
Per il calcolo si considerano i punti con Z>THR_VPR, qualità>QMIN_VPR, BeamBlocking<20% e clutter free all'interno del volume scelto.
Il punto del vpr corrispondente ad un livello è calcolato come la somma dei volumi di pioggia di ciascuna cella polare presente nell'intervallo di quota, diviso l'area totale. Per ogni livello devo avere un'area minima coperta di pioggia,se non si verifica, ricopio il valore del primo livello che raggiunge questa copertura.
La parte del profilo alta, che sta cioè al di sopra dell'ultimo livello calcolato, viene riempite tirando una retta con coeff. amgolre negativo e costante, pari a -0.03 (valore passibile di discussione)
Il profilo è poi soggetto a quality check e viene rigettato (return(1)) se:
- estensione verticale troppo bassa (devo avere almeno 2 km di profilo da un punto che sta non più alto di 700 m dal suolo)
- cv<CT_MIN*ct, cioè frazione di volume precipitante piccola
- la deviazione standard del volume di R a 1100 m normalizzata rispetto al volume totale precipitante supera soglia prefissata


COSTANTI
THR_VPR:soglia valore Z per calcolo vpr
QMIN_VPR: qualità minima necessaria per entrare nel calcolo vpr
MIN_AREA: area minima a un livello per avere un punto del vpr
IAZ_MIN_SPC,IAZ_MAX_SPC,IAZ_MIN_GAT,I_AZ_MAX_GAT,R_MIN_VPR,R_MAX_VPR:limiti dell'area considerata per il calcolo VPR
TCK_VPR: spessore dei livelli
NMAXLAYER: n max. livelli
VEXTMIN_VPR: estensione verticale minima del profilo
THR_TDEVR: soglia di massima dev. standard di R per considerare buono il profilo
CT_MIN: frazione di volume minimo precipitante per considerare buono il profilo

VARIABILI
int flag_vpr[][][]: flag che indica l'usabilità della cella polare in termini di posizione del bin, beam blocking e clutter
long int *cv,*ct:  volume del settore libero,volume precipitante
float vpr1[]: vpr istantaneo
long int vert_ext,vol_rain: estensione verticale profilo, volume pioggia del singolo bin
long int area_vpr[NMAXLAYER]; area totale usata per calcolo vpr

*/
int CalcoloVPR::func_vpr(long int *cv, long int *ct, vector<float>& vpr1, vector<long int>& area_vpr)
{
    LOG_CATEGORY("radar.vpr");
    int i,iA,ilay,il,ilast,iaz_min,iaz_max;
    long int dist,dist_plain,vert_ext,vol_rain;
    float area,quota_true_st;
    float grad;

    //----------------inizializzazioni----------------
    stdev = -1.;
    vert_ext=0;

    *cv=0;
    *ct=0;


    /*SECTION A: cv and ct retrieval and calculation of current area- weighted vpr*/


    //------------riconoscimento sito per definizione limiti azimut---------
    iaz_min=cum_bac.site.vpr_iaz_min;
    iaz_max=cum_bac.site.vpr_iaz_max;

    for (unsigned l=0; l<cum_bac.volume.size(); l++)//ciclo elevazioni
    {
        const PolarScan<double>& scan = cum_bac.volume.scan(l);
        const volume::PolarScanLoadInfo& scan_info = cum_bac.load_info.scan(l);

        for (unsigned k=0; k < scan.beam_size; k++)/*ciclo range*/
        {
            //-------------calcolo distanza-----------
            dist=k*(long int)(cum_bac.load_info.size_cell)+(int)(cum_bac.load_info.size_cell)/2.;

            //-----ciclo settore azimut???
            for (iA=iaz_min; iA<iaz_max; iA++)//ciclo sulle unità di azimut  (0,9°)    ---------
            {
                //----------ottengo il valore dell'indice in unità di azimut (0,9°) ---------
                i=(iA+NUM_AZ_X_PPI)%NUM_AZ_X_PPI;

                //--------calcolo elevazione e quota---------
                // FIXME: this reproduces the truncation we had by storing angles as short ints between 0 and 4096
                //elevaz=(float)(cum_bac.volume.scan(l)[i].teta_true)*CONV_RAD;
                //elevaz=(float)(cum_bac.volume.scan(l)[i].elevation*DTOR);
                const float elevaz = scan_info.get_elevation_rad(i);
                quota_true_st=cum_bac.quota_f(elevaz,k);

                //--------trovo ilay---------
                ilay=floor(quota_true_st/TCK_VPR);//  in teoria quota indipendente da azimuth  , in realtà no (se c'è vento)

                if (ilay <0 || ilay >= NMAXLAYER) {
                    //fprintf(log_vpr,"ilay %d errore\n",ilay);
                    break;
                }


                vol_rain=0;
                // dist=(long int)(dist*cos((float)(cum_bac.volume_at_elev_preci(i, k).teta_true)*CONV_RAD));

                /* //---------calcolo la distanza proiettata sul piano-------------  */

                dist_plain=(long int)(dist*cos(elevaz));
                if (dist_plain <RMIN_VPR || dist_plain > RMAX_VPR )
                    flag_vpr->scan(l).set(i, k, 0);

                if (cum_bac.qual->scan(l).get(i, k) < QMIN_VPR) flag_vpr->scan(l).set(i, k, 0);

                //AGGIUNTA PER CLASS
                if(cum_bac.do_class){
                    if(conv[i][k]>= CONV_VAL){
                        flag_vpr->scan(l).set(i, k, 0);
                    }
                }

                // ------per calcolare l'area del pixel lo considero un rettangolo dim bin x ampiezzamediafascio x flag vpr/1000 per evitare problemi di memoria?
                area=cum_bac.load_info.size_cell*dist_plain*AMPLITUDE*DTOR*flag_vpr->scan(l).get(i, k)/1000.; // divido per  mille per evitare nr troppo esagerato

                // ------incremento il volume totale di area

                *cv=*cv+(long int)(area);


                //---------------------condizione per incrementare VPR contributo: valore sopra 13dbz, qualità sopra 20 flag>0 (no clutter e dentro settore)------------------
                double sample = scan.get(i, k);
                if (sample > THR_VPR &&  flag_vpr->scan(l).get(i, k) > 0 )
                {
                    //-------incremento il volume di pioggia = pioggia x area
                    vol_rain=(long int)(BYTE_to_mp_func(DBtoBYTE(sample),cum_bac.aMP,cum_bac.bMP)*area);//peso ogni cella con la sua area

                    //-------incremento l'area precipitante totale ct,aggiungendo però,cosa che avevo messo male una THR solo per ct, cioè per il peso
                    if (sample > THR_PDF)
                        *ct=*ct+(long int)(area);

                    //------se l'area in quello strato è già maggiore di 0 allora incremento il volume dello strato altrimenti lo scrivo ex novo. poi vpr1 andrà diviso per l'area
                    if (area_vpr[ilay]> 0) vpr1[ilay]=vpr1[ilay]+(float)(vol_rain);
                    else vpr1[ilay]=(float)(vol_rain);

                    //------incremento l'area dello strato----------
                    area_vpr[ilay]=area_vpr[ilay]+area;
                }
            }
        }
    }

    LOG_INFO("calcolati ct e cv ct= %li cv= %li",*ct,*cv);

    /*SECTION B: vpr quality checks and re-normalisation of vpr*/

    //--------------CONTROLLO DI QUALITA' E NORMALIZZAZIONE DEL PROFILO ISTANTANEO CALCOLATO
    //-------- se il volume supera quello minimo------
    if ((*ct) > CT_MIN*(*cv)) {

        ilast=0;
        vert_ext=0;


        //----- calcolo 'estensione verticale del profilo , negli strati dove l'area è troppo piccola assegno NODATAVPR,  NORMALIZZO il profilo, e se l'estensione è minore di VEXTMIN_VPR esco-


        for (ilay=0; ilay<NMAXLAYER; ilay++){


            LOG_INFO("  ilay %d area_vpr= %ld  ct= %ld  cv= %ld", ilay, area_vpr[ilay],*ct,*cv );

            if (area_vpr[ilay]>=MIN_AREA) {
                vert_ext=vert_ext+TCK_VPR;
                vpr1[ilay]=vpr1[ilay]/(float)(area_vpr[ilay]);

            }
            else
            {
                vpr1[ilay]=NODATAVPR;


                //----  se incontro un punto vuoto oltre 700 m ( o se sono arrivata alla fine) assegno ilast ed esco dal ciclo

                /*   if (ilast > 0 && vert_ext>VEXTMIN_VPR){ */

                if (ilay*TCK_VPR+TCK_VPR/2 > MINPOINTBASE || ilay== NMAXLAYER -1){  //cambio il criterio per calcolare estensione minima. devo avere almeno 2 km consecutivi a partire da un punto che stia almeno a 700 m (MINPOINTBASE),
                    LOG_INFO("raggiunta cima profilo");
                    ilast=ilay-1;// c'era errore!!!

                    //---------- raggiunta la cima profilo faccio check immediato sull'estensione verticale
                    if (vert_ext<VEXTMIN_VPR ){
                        LOG_INFO("estensione profilo verticale troppo bassa");
                        *ct=0;
                        ilast=0;
                        for  (il=0; il<ilast; il++) vpr1[il]=NODATAVPR;
                        return(1);
                    }

                    break; // esco dal ciclo..modifica
                }
            }
        }
    }// fine se volumeprecipitante sufficiente

    // ---------se il volume non supera quello minimo esco---------
    else {
        LOG_INFO("volume precipitante troppo piccolo");
        *ct=0;
        ilast=0;
        for  (il=0; il<NMAXLAYER; il++) vpr1[il]=NODATAVPR; //--devo riassegnare o mi rimane 'sporco' forse si potrebbe usare una ver diversa
        return(1);
    }


    //------calcolo il gradiente del profilo come media del gradiente negli ultimi 3 strati per assegnare la parte 'alta' (novità)

    grad=((vpr1[ilast]-vpr1[ilast-1]) + (vpr1[ilast-1]-vpr1[ilast-2])/(2.)+ (vpr1[ilast-2]-vpr1[ilast-3])/(3.) ) /3.;
    if (grad > 0.002)
        grad=-0.03 ;
    LOG_INFO(" %f", grad);

    //--riempio la parte alta del profilo decrementando di grad il profilo in ogni strato fino a raggiunere 0, SI PUÒ TOGLIERE E METTERE NODATA

    for (ilay=ilast+1; ilay<NMAXLAYER; ilay++) {
        if (vpr1[ilay-1] + grad > 0.002)
            vpr1[ilay]= vpr1[ilay-1]+grad;
        else
            vpr1[ilay]=0;
    }

    // HO CAMBIATO DA GRADIENTE FISSO PARI A (V. VECCHIO) A GRADIENTE RICAVATO DAL PROFILO PER LA PARTE ALTA
    //---HO TOLTO TUTTA LA PARTE CHE FA IL CHECK SULLA STDEV A 1100 M SI PUO' RIVEDERE SE METTERLA SE SERVE.


    return(0);
}

float comp_levels(float v0, float v1, float nodata, float peso)
{
    float result;
    /* if ((v0<nodata+1)&&(v1<nodata+1)) result=nodata; */
    /* if (v0<nodata+1) result=v1; */
    /* if (v1<nodata+1) result=v0;         */
    if ((v0>nodata) && (v1>nodata)  ) result=((1.-peso)*v0+peso*v1); /* in questa configurazione il vpr è di altezza costante  nel tempo ma un po' 'sconnesso' in alto*/

    else result=nodata;
    return(result);
}

void CUM_BAC::conversione_convettiva()
{
    for (unsigned i=0; i<NUM_AZ_X_PPI; i++){
        for (unsigned k=0; k<volume.scan(0).beam_size; k++){
            if (calcolo_vpr->conv[i][k] > 0){
                volume.scan(0).set(i, k,
                        ::RtoDBZ(
                            BYTE_to_mp_func(
                                DBtoBYTE(volume.scan(0).get(i, k)),
                                aMP_conv,
                                bMP_conv),
                            aMP_class,
                            bMP_class));

            }
        }
    }
}

bool CUM_BAC::esegui_tutto(const char* nome_file, int file_type)
{
    // Legge e controlla il volume dal file SP20
//    if (!read_sp20_volume(nome_file, file_type))
//        return false;
    // Legge e controlla il volume dal file ODIM
    if (!read_odim_volume(nome_file, file_type))
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
    if (NormalizzoData(load_info.acq_date) == -1)
        return true;

    setup_elaborazione(nome_file);
	printwork();

    //--------------se def anaprop : rimozione propagazione anomala e correzione beam blocking-----------------//
    LOG_INFO("inizio rimozione anaprop e beam blocking");
    elabora_dato();

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
        conversione_convettiva();

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

    /*--------------------------------------------------
      | conversione di coordinate da polare a cartesiana |
      --------------------------------------------------*/
    LOG_INFO("Creazione Matrice Cartesiana");
    Cart cart_maker(volume.max_beam_size());
    cart_maker.creo_cart(*this);


    //-------------------Se definita Z_LOWRIS creo matrice 1X1  ZLR  stampo e stampo coeff MP (serve?)------------------

    LOG_INFO("Estrazione Precipitazione 1X1");
    CartLowris cart_low(do_medium ? 512: 256);
    cart_low.creo_cart_z_lowris(*this, cart_maker);

    unsigned char MP_coeff[2]; /* a/10 e b*10 per scrivere come 2 byte */
    MP_coeff[0]=(unsigned char)(aMP/10);
    MP_coeff[1]=(unsigned char)(bMP*10);

    char nome_file_output[512];
    sprintf(nome_file_output,"%s/MP_coeff",getenv("OUTPUT_Z_LOWRIS_DIR"));
    FILE* output=controllo_apertura(nome_file_output,"file coeff MP","w");
    fwrite(MP_coeff,sizeof(MP_coeff),1,output);
    fclose(output);
    printf(" dopo scrivo_z_lowris\n");

    if (do_quality && do_devel)
    {
        H5::H5File outfile = assets.get_devel_data_output();
        //scrivo_out_file_bin(".corrpt","file anap",getenv("DIR_QUALITY"),sizeof(dato_corrotto),dato_corrotto);
        //scrivo_out_file_bin(".bloc","file bloc",getenv("DIR_QUALITY"),sizeof(beam_blocking),beam_blocking);
        //scrivo_out_file_bin(".quota","file quota",getenv("DIR_QUALITY"),sizeof(quota),quota);
        //scrivo_out_file_bin(".elev","file elevazioni",getenv("DIR_QUALITY"),sizeof(elev_fin),elev_fin);
        elev_fin.write_info_to_debug_file(outfile);
    }

    LOG_INFO("Scrittura File Precipitazione 1X1 %s\n", nome_file);
    cart_low.write_out(*this, assets);

    return true;
}


CalcoloVPR::CalcoloVPR(CUM_BAC& cum_bac)
    : cum_bac(cum_bac),
      area_vpr(NMAXLAYER, 0),
      vpr(NMAXLAYER, NODATAVPR),
      corr_polar(cum_bac.volume.max_beam_size()),
      neve(cum_bac.volume.max_beam_size()),
      flag_vpr(0)
{
    logging_category = log4c_category_get("radar.vpr");
    MyMAX_BIN=cum_bac.MyMAX_BIN;
    htbb=-9999.; hbbb=-9999.;
    t_ground=NODATAVPR;

    /*
    for (int k=0; k<NUM_AZ_X_PPI*MyMAX_BIN;k++ ){
      lista_conv[k][0]=-999;
      lista_conv[k][1]=-999;
    }
    */

    flag_vpr = new Volume<unsigned char>(cum_bac.volume, 0);

    if (cum_bac.do_vpr)
        t_ground = cum_bac.assets.read_t_ground();
}

CalcoloVPR::~CalcoloVPR()
{
    if (flag_vpr) delete flag_vpr;
}

void CalcoloVPR::esegui_tutto()
{
    test_vpr=fopen(getenv("TEST_VPR"),"a+");

    LOG_INFO("processo file dati: %s", cum_bac.load_info.filename.c_str());
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

double CUM_BAC::BeamBlockingCorrection(double val_db, unsigned char beamblocking)
{
   return val_db - 10 * log10(1. - (double)beamblocking / 100.);
}

float CUM_BAC::RtoDBZ(float rain) const
{
    return ::RtoDBZ(rain, aMP, bMP);
}

/* time è in secondi, itime è un intero che rappresenta il numero intero di intervalli da 5 minuti*/
time_t CUM_BAC::NormalizzoData(time_t time)
{
    //Parametri passare da minuti del file a minuti standard arrotondando per difetto o eccesso ( prima si arrotondava al 5° ora si arrotonda al minuto )
    // massima differenza in minuti tra data acquisizione e standard per arrotondare per difetto
    const unsigned MAX_TIME_DIFF = do_medium ? 1 : 3;
    const unsigned NMIN = 1;
    int itime;

    itime = time/(NMIN*60);

    if(time - itime*NMIN*60 <MAX_TIME_DIFF*60) return (itime*NMIN*60); /* se la differenza è meno di tre minuti vado al 5° min. prec*/
    if(time - itime*NMIN*60 >(NMIN-MAX_TIME_DIFF)*60) return ((itime+1)*NMIN*60); /* se la differenza è più di tre minuti vado al 5° min. successivo*/
    //altrimenti ritorno -1
    return -1;
}

Cart::Cart(unsigned max_bin)
    : max_bin(max_bin),
      cart(max_bin*2), 
      beam_blocking_xy(max_bin*2),
      dato_corr_xy(max_bin*2),
      elev_fin_xy(max_bin*2),
      topxy(max_bin*2),
      qual_Z_cart(max_bin*2),
      corr_cart(max_bin*2),
      neve_cart(max_bin*2),
      conv_cart(max_bin*2),
      cappi_cart(max_bin*2),
      quota_cart(max_bin*2),
      cartm(max_bin*2)
{
}

namespace {
struct CartData
{
    Image<double> azimut;
    Image<double> range;

    CartData(int max_bin=512)
        : azimut(max_bin), range(max_bin)
    {
        for(int i=0; i<max_bin; i++)
            for(int j=0; j<max_bin; j++)
            {
                range(i, j) = hypot(i+.5,j+.5);
                azimut(i, j) = 90. - atan((j+.5)/(i+.5)) * M_1_PI*180.;
            }
    }
};
}

void Cart::creo_cart(const CUM_BAC& cb)
{
    LOG_CATEGORY("radar.cart");

    //matrici per ricampionamento cartesiano
    //int x,y,irange,az,iaz,az_min,az_max,cont;
    int x,y,iaz,az_min,az_max,cont;
    float az;
    CartData cd(max_bin);

    for(unsigned i=0; i<max_bin *2; i++)
        for(unsigned j=0; j<max_bin *2; j++)
            cart(i, j) = MISSING;

    LOG_INFO("Creo_cart - %u", max_bin);

    for(unsigned quad=0; quad<4; quad++)
        for(unsigned i=0; i<max_bin; i++)
            for(unsigned j=0; j<max_bin; j++)
            {
                unsigned irange = (unsigned)round(cd.range(i, j));
                if (irange >= max_bin)
                    continue;
                switch(quad)
                {
                    case 0:
                        x = max_bin + i;
                        y = max_bin + j;
                        az = cd.azimut(i, j);
                        break;
                    case 1:
                        x = max_bin + j;
                        y = max_bin - i;
                        az = cd.azimut(i, j) + 90.;
                        break;
                    case 2:
                        x = max_bin - i;
                        y = max_bin - j;
                        az = cd.azimut(i, j) + 180.;
                        break;
                    case 3:
                        x = max_bin - j;
                        y = max_bin + i;
                        az = cd.azimut(i, j)+270.;
                        break;
                }

                az_min = (int)((az - .45)/.9);
                az_max = ceil((az + .45)/.9);


                if(az_min < 0)
                {
                    az_min = az_min + NUM_AZ_X_PPI;
                    az_max = az_max + NUM_AZ_X_PPI;
                }
                cont=0;
                for(iaz = az_min; iaz<az_max; iaz++){
                    // Enrico: cerca di non leggere fuori dal volume effettivo
                    unsigned char sample = 0;
                    if (irange < cb.volume.scan(0).beam_size)
                        sample = DBtoBYTE(cb.volume.scan(0).get(iaz%NUM_AZ_X_PPI, irange));

                    if(cart(x, y) <= sample){
                        cart(x, y) = sample;
                        topxy(x, y)=cb.top(iaz%NUM_AZ_X_PPI, irange);
                        if (cb.do_quality)
                        {
                            if (irange < cb.volume.scan(cb.elev_fin[iaz%NUM_AZ_X_PPI][irange]).beam_size)
                                qual_Z_cart(x, y) = cb.qual->scan(cb.elev_fin[iaz%NUM_AZ_X_PPI][irange]).get(iaz%NUM_AZ_X_PPI, irange);
                            else
                                qual_Z_cart(x, y) = 0;
                            quota_cart(x, y)=cb.quota(iaz%NUM_AZ_X_PPI, irange);
                            dato_corr_xy(x, y)=cb.dato_corrotto(iaz%NUM_AZ_X_PPI, irange);
                            beam_blocking_xy(x, y)=cb.beam_blocking(iaz%NUM_AZ_X_PPI, irange);
                            elev_fin_xy(x, y)=cb.elev_fin[iaz%NUM_AZ_X_PPI][irange];
                            /*neve_cart(x, y)=qual_Z_cart(x, y);*/
                            if (cb.do_vpr)
                            {
                                neve_cart(x, y)=(cb.calcolo_vpr->neve(iaz%NUM_AZ_X_PPI, irange))?0:1;
                                corr_cart(x, y)=cb.calcolo_vpr->corr_polar(iaz%NUM_AZ_X_PPI, irange);
                            }
                        }
                        if (cb.do_class)
                        {
                            if (irange<cb.calcolo_vpr->x_size)
                                conv_cart(x, y)=cb.calcolo_vpr->conv[iaz%NUM_AZ_X_PPI][irange];
                        }
                    }
                    if (cb.do_zlr_media)
                    {
                        //if (volume.scan(0).get_raw(iaz%NUM_AZ_X_PPI, irange) > 0)
                        if (sample > 0)
                    {
                            cartm(x, y)=cartm(x, y)+BYTEtoZ(sample);
                            cont=cont+1;
                        }
                    }
                }

                if (cb.do_zlr_media)
                {
                    if (cont > 0) cartm(x, y)=cartm(x, y)/(float)(cont);
                }
                /*
                 *****  per scrivere in griglia cartesiana************
                 bloc_xy[MAX_BIN*2-y][x]=beam_blocking((int)((float)(az)/.9), irange);
                 elev_xy[MAX_BIN*2-y][x]=volume.elev_fin[(int)((float)(az)/.9)][irange];
                 dato_corrotto_xy(MAX_BIN*2-y, x)= dato_corrotto((int)((float)(az)/.9), irange);
                 elev_fin_xy[MAX_BIN*2-y][x]=first_level((int)((float)(az)/.9), irange);
                 dato_corrotto_xy(MAX_BIN*2-y, x)= dato_corrotto((int)((float)(az)/.9), irange); */
            }
}

CartLowris::CartLowris(unsigned cart_dim_zlr)
    : CART_DIM_ZLR(cart_dim_zlr),
      z_out(CART_DIM_ZLR),
      quota_1x1(CART_DIM_ZLR),
      beam_blocking_1x1(CART_DIM_ZLR),
      dato_corr_1x1(CART_DIM_ZLR),
      elev_fin_1x1(CART_DIM_ZLR),
      qual_Z_1x1(CART_DIM_ZLR),
      top_1x1(CART_DIM_ZLR),
      corr_1x1(CART_DIM_ZLR),
      neve_1x1(CART_DIM_ZLR),
      conv_1x1(CART_DIM_ZLR),
      cappi_1x1(CART_DIM_ZLR)
{
}

void CartLowris::creo_cart_z_lowris(const CUM_BAC& cb, const Cart& c)
{
    unsigned ZLR_OFFSET = cb.do_medium && cb.MyMAX_BIN != 1024 ? CART_DIM_ZLR/2 : 0;
    unsigned ZLR_N_ELEMENTARY_PIXEL = cb.do_medium && cb.MyMAX_BIN != 1024 ? 1 : 4;

    //tolta qui inizializzazione di z_out che era duplicata (già fatta all'inizio del main)
    // ciclo sui punti della nuova matrice. per il primo prenderò il massimo tra i primi sedici etc..
    for(unsigned i=0; i<CART_DIM_ZLR; i++)
        for(unsigned j=0; j<CART_DIM_ZLR; j++)
        {
            unsigned int cont=0;
            double zm = 0.;
            unsigned char z = 0;
            //reinizializzo tutte le variabili calcolate dentro la funzione .
            unsigned char q = 0;
            unsigned char nv = 0;
            unsigned char dc1x1=0;
            unsigned char el1x1=0;
            unsigned short q1x1=0;
            unsigned char c1x1=0;
            unsigned char bl1x1=0;
            unsigned char traw=0;
            for(unsigned x = 0; x < ZLR_N_ELEMENTARY_PIXEL; x++)
                for(unsigned y = 0; y < ZLR_N_ELEMENTARY_PIXEL; y++)
                    //ciclo a passi di 4 in x e y nella matrice a massima risoluzione, cercando il valore massimo di z tra i primi sedici e attribuendolo al primo punto della matrice a bassa risoluzione e poi i tra i secondi sedici e attribuendolo al secondo punto etc...
                {
                    unsigned src_x = i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET;
                    unsigned src_y = j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET;
                    if (src_x < c.max_bin*2 && src_y < c.max_bin*2 && c.cart(src_x, src_y) != MISSING)
                    {
                    //if (src_x >= cb.MyMAX_BIN*2) printf("X è fuori\n");
                    //if (src_y >= cb.MyMAX_BIN*2) printf("Y è fuori %d %d\n",src_y, CART_DIM_ZLR);
                    //if(c.cart(src_x, src_y) != MISSING)
                        if(c.cart(src_x, src_y) > z){
                            z= c.cart(src_x, src_y);
                            traw=c.topxy(src_x, src_y);
                            if (cb.do_quality)
                            {
                                q=c.qual_Z_cart(src_x, src_y);
                                q1x1=c.quota_cart(src_x, src_y);
                                dc1x1=c.dato_corr_xy(src_x, src_y);
                                el1x1=c.elev_fin_xy(src_x, src_y);
                                bl1x1=c.beam_blocking_xy(src_x, src_y);

                                if (cb.do_vpr)
                                {
                                    c1x1=c.corr_cart(src_x, src_y);
                                    nv= c.neve_cart(src_x, src_y);
                                }
                            }

                            if (cb.do_class)
                            {
                                conv_1x1(i, j)=c.conv_cart(src_x, src_y);
                            }

                            if (cb.do_zlr_media)
                            {
                                if (c.cartm(src_x, src_y) > 0) {
                                    zm = zm + c.cartm(src_x, src_y);
                                    cont=cont+1;
                                }
                            }
                        }
                    }
                    z_out(i, j)=z;
                    if (cb.do_quality)
                    {
                        qual_Z_1x1(i, j)=q;
                        quota_1x1(i, j)=128+(unsigned char)(q1x1/100);
                        dato_corr_1x1(i, j)=dc1x1;
                        elev_fin_1x1(i, j)=el1x1;
                        beam_blocking_1x1(i, j)=bl1x1;
                    }
                    top_1x1(i, j)=traw;

                    if (cb.do_vpr)
                    {
                        neve_1x1(i, j)=nv;
                        corr_1x1(i, j)=c1x1;
                    }

                    if (cb.do_zlr_media)
                    {
                        if (cont >0 ) {
                            z_out(i, j)=(unsigned char)round((10*log10(zm/(float)(cont))+20.)/80.*255);
                        }
                        if (cont == 0 ) z_out(i, j)=MISSING;
                    }
                }
        }
}


void CartLowris::write_out(const CUM_BAC& cb, Assets& assets)
{
    //-------------------scritture output -----------------------
    assets.write_image(z_out, "OUTPUT_Z_LOWRIS_DIR", ".ZLR", "file output 1X1");

    //-------------------scritture output -----------------------
    if (cb.do_quality)
    {
        //------------------ output qual  per operativo:qualità in archivio e elevazioni e anap in dir a stoccaggio a scadenza
        // temporanee
        // in archivio
        assets.write_image(qual_Z_1x1, "OUTPUT_Z_LOWRIS_DIR", ".qual_ZLR", "file qualita' Z");

        //------------------ stampe extra per studio
        if (cb.do_devel)
        {
            assets.write_image(beam_blocking_1x1, "DIR_QUALITY", ".bloc_ZLR", "file bloc");
            assets.write_image(dato_corr_1x1, "DIR_QUALITY", ".anap_ZLR", "file anap");
            assets.write_image(quota_1x1, "DIR_QUALITY", ".quota_ZLR", "file qel1uota");
            assets.write_image(elev_fin_1x1, "DIR_QUALITY", ".elev_ZLR", "file elev");
            assets.write_image(top_1x1, "DIR_QUALITY", ".top20_ZLR", "file top20");
        }

        //------------------ stampe correzioni da profili verticali in formato ZLR
        if (cb.do_vpr)
        {
            assets.write_image(corr_1x1, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
            //scrivo_out_file_bin(".neve","punti di neve",getenv("DIR_QUALITY"),sizeof(neve),neve);
            // scrivo_out_file_bin(".neve_ZLR","file presunta neve ",getenv("DIR_QUALITY"),sizeof(neve_1x1),neve_1x1);
        }

        //------------------se definita CLASS  stampo punti convettivi
        if (cb.do_class)
        {
            // in archivio?
            // scrivo_out_file_bin(".conv_ZLR","punti convettivi",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(conv_1x1),conv_1x1);
            assets.write_image(conv_1x1, "DIR_QUALITY", ".conv_ZLR", "punti convettivi");
        }
    }
}


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

    LOG_INFO("tempo parziale %ld ---- totale %ld", time2-time1, time2-time_tot);
    time1=time2;
    return;
}
