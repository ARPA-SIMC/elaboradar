#ifndef ARCHIVIATORE_CUM_BAC_CLASS_H
#define ARCHIVIATORE_CUM_BAC_CLASS_H

/**
 *  @file
 *  @defgroup progetto_cum_bac
 *  @brief progetto elaborazione dati radar per ottenere campi di Z da convertire in R
 *  @details elaborazione dati radar utilizzando dati da radiosondaggio, temperatura e dati in uscita da programma beam blocking che esegue un controllo di qualita' del volume, rimuove la propagazione anomala, corregge il beam blocking,   calcola il profilo verticale, calcola la  qualità, classifica le aree convettive  
*/


/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief codice principale di elaborazione dei volumi di riflettivita' radar usato per impulso corto
 *  @details questo codice contiene il main e alcune funzioni usate nell'elaborazione della riflettivita' radar 
 *  librerire necessarie: lib_SP20
*/

#include "logging.h"
#include "assets.h"
#include "volume.h"
#include <stdexcept>
#include <cmath>

//algoritmo
#include <MP_par.h>
#include <vpr_par.h>
#include <geo_par.h>


#define MAX_BIN 512

// Toglie NEL come costante
#ifdef NEL
#undef NEL
#endif

//Risoluzioni e limiti spaziali
// FIXME: rimuovere
#ifdef MEDIUM
const unsigned int NEL = 5;
#else
const unsigned int NEL = 15;  // n0 elevazioni massimo
#endif

//Dimensioni matrici statistica
#define DIM1_ST 16
#define DIM2_ST 13/*Cambiata dimensione a 13 per cambio dimensione raggio radar*/

namespace cumbac {

struct Site;

template<typename T>
struct PolarMap
{
    const unsigned beam_size;
    const unsigned beam_count;
    T* data;

    PolarMap(const PolarMap& pm)
        : beam_size(pm.beam_size), beam_count(pm.beam_count), data(new T[beam_size * beam_count])

    {
        for (unsigned i = 0; i < beam_size * beam_count; ++i)
            data[i] = pm.data[i];
    }

    PolarMap(unsigned beam_size=512, unsigned beam_count=400)
        : beam_size(beam_size), beam_count(beam_count), data(new T[beam_size * beam_count])
    {
        for (unsigned az = 0; az < beam_count; ++az)
            for (unsigned i = 0; i < beam_size; ++i)
                data[az * beam_size + i] = 0;
    }

    ~PolarMap()
    {
        delete[] data;
    }

    T* operator[](unsigned az) { return data + az * beam_size; }
    const T* operator[](unsigned az) const { return data + az * beam_size; }

    PolarMap& operator=(const PolarMap& pm)
    {
        if (beam_size != pm.beam_size or beam_count != pm.beam_count)
            throw std::runtime_error("PolarMap size mismatch");
        for (unsigned i = 0; i < beam_size * beam_count; ++i)
            data[i] = pm.data[i];
        return *this;
    }
};

template<typename T>
struct Image
{
    const unsigned SX;
    const unsigned SY;
    T* data;

    Image(unsigned sx, unsigned sy=0)
        : SX(sx), SY(sy == 0 ? sx : sy), data(new T[SX*SY])
    {
        for (unsigned y = 0; y < SY; ++y)
            for (unsigned x = 0; x < SX; ++x)
                data[y*SX+x] = 0;
    }
    ~Image()
    {
        delete[] data;
    }

    T* operator[](unsigned y) { return data + y * SX; }
    const T* operator[](unsigned y) const { return data + y * SX; }

    size_t size() const { return SX * SY; }

    T min() const
    {
        T res = data[0];
        for (unsigned i = 0; i < SX * SY; ++i)
            if (data[i] < res)
                res = data[i];
        return res;
    }

    T max() const
    {
        T res = data[0];
        for (unsigned i = 0; i < SX * SY; ++i)
            if (data[i] > res)
                res = data[i];
        return res;
    }

    T avg() const
    {
        double res = 0;
        for (unsigned i = 0; i < SX * SY; ++i)
            res += (double)data[i] / (double)(SX * SY);
        return (T)round(res);
    }
};

struct CalcoloVPR;

class CUM_BAC
{
public:
    log4c_category_t* logging_category;

    const Site& site;
    Assets assets;

    bool do_medium;

    /// Feature set required for this run
    bool do_clean;        // Clean and truncate input volume
    bool do_quality;
    bool do_beamblocking;
    bool do_declutter;
    bool do_bloccorr;
    bool do_vpr;
    bool do_class;
    bool do_zlr_media;
    bool do_devel;
    bool do_readStaticMap;  

    // dimensione matrice a 1x1 km
    const unsigned CART_DIM_ZLR;

    cumbac::Volume volume;

    CalcoloVPR* calcolo_vpr;

    /*-----------------------------------------------------------
      Variabili globali
      ------------------------------------------------------------*/

    // vol_pol riportato in cartesiano
    Image<unsigned char> cart;
    Image<double> cartm;  /* Z media dei bins adiacenti al punto */
    //se definita Z_LOWRIS, Z cartesiana al livello più basso
    Image<unsigned char> z_out;

    //variabili tempo per ottenere mese.. aggiunti nel main per leggere stagione dal nome file e ricavere MP coeff */
    struct tm *tempo;
    time_t Time;
    int month;
    char date[20];

    //coeff a e b relazione Z-R
    unsigned char MP_coeff[2]; /* a/10 e b*10 per scrivere come 2 byte */
    float aMP, bMP;   /*  coeff a e b relazione Z-R  */


    //metrici per statistiche
    int stat_anap[DIM1_ST][DIM2_ST]; /* statistica anaprop  */
    int stat_anap_tot[DIM1_ST][DIM2_ST]; /* contatore punti dentro ogni box per statistica */
    long int stat_bloc[DIM1_ST][DIM2_ST];   /* statistica beam blocking  */
    int stat_elev[DIM1_ST][DIM2_ST]; /* statistica cambio elevazione rispetto mappa statica  */

    //matrici first_level e first level da beam blocking e valore beam blocking
    PolarMap<unsigned char> first_level; //mappa dinamica complessiva
    PolarMap<unsigned char> first_level_static;//mappa statica

    PolarMap<unsigned char> bb_first_level;  /* mappa di elevazioni da beam blocking (input)*/
    PolarMap<unsigned char> beam_blocking;   /* mappa di beam blocking (input)*/

    //variabili legate a propagazione e beam blocking, da prog_bb
    // quota centro fascio in funzione della distanza e elevazione
    float hray[MAX_BIN][NEL];
    float hray_inf[MAX_BIN][NEL]; /*quota limite inferiore fascio in funzione della distanza e elevazione*/
    float dem[NUM_AZ_X_PPI][MAX_BIN]; /*dem in coordinate azimut range*/
    float dtrs;// distanza temporale radiosondaggio

    // attenuazione in formato cartesiano max risoluzione
    PolarMap<unsigned char> att_cart; /* matrice azimut-range di attenuazione */
    //quota centro fascio polare, cartesiana max risoluzione e cartesiana 1x1
    unsigned short quota_rel[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio relativa al suolo in prop da rsd e elevazioni nominali, in coordinate azimut range*/
    unsigned short quota[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio in prop standard e elev reali in coordinate azimut range*/
    Image<unsigned short> quota_cart;/*quota fascio in coordinate cart 1024*1024, risoluzione minima*/
    Image<unsigned char> quota_1x1;/* quota in formato 256*256 in centinaia di metri, risoluzione ZLR */
    //beam blocking cartesiano max resol e 1x1
    Image<unsigned char> beam_blocking_xy; //beamblocking cartesiano max resol
    Image<unsigned char> beam_blocking_1x1;//beam blocking cartesiano 1x1
    //uscite anaprop
    unsigned char dato_corrotto[NUM_AZ_X_PPI][MAX_BIN]; /*uscita controllo anaprop in coordinate azimut range */
    Image<unsigned char> dato_corr_xy; //uscite anap  cartesiano max resol
    Image<unsigned char> dato_corr_1x1; //uscite anap cartesiano  1x1
    Image<unsigned char> elev_fin_xy;
    Image<unsigned char> elev_fin_1x1;
    // metrici qualita' come sopra
    VolumeInfo<unsigned char>* qual; // qualita volume polare
    Image<unsigned char> qual_Z_cart; /* qualita della Z in formato 1024*1024, risoluzione minima */
    Image<unsigned char> qual_Z_1x1;/* qualita della Z in formato 256*256, risoluzione ZLR */
    // top, come sopra
    unsigned char top[NUM_AZ_X_PPI][MAX_BIN];
    Image<unsigned char> topxy;
    Image<unsigned char> top_1x1;

    // uscite  vpr: correzione VPR , come sopra
    Image<unsigned char> corr_cart;
    Image<unsigned char> corr_1x1;
    // uscite vpr: neve, come sopra
    Image<unsigned char> neve_cart;/* neve formato 1024*1024, risoluzione minima */
    Image<unsigned char> neve_1x1;/* neve in formato 256*256, risoluzione ZLR */

    //matrici per classificazione: cappi
    unsigned char cappi[NUM_AZ_X_PPI][MAX_BIN];
    // uscite: matrici class max resol e 1x1
    Image<unsigned char> conv_cart;
    Image<unsigned char> conv_1x1;
    //uscite:matrici cappi max resol e 1x1
    Image<unsigned char> cappi_cart;
    Image<unsigned char> cappi_1x1;

    /* variabili tolte perchè non presenti nel codice cum_bac... controllare che non richiamino qualcosa nelle funzioni
       struct tm *time_dbp;
       T_time, T_data, T_ora..*/


    CUM_BAC(const char* site_name, bool medium=false);
    ~CUM_BAC();

    bool read_sp20_volume(const char* nome_file, int file_type);
    bool read_odim_volume(const char* nome_file, int file_type);

    /**
     *  @brief funzione che verifica se il volume dati in input e' consistente con le esigenze del programma 
     *  @details  TESTA: A) risoluzione B) numero raggi per elevazione
     *  @param    tipofile  puo' essere 0(corto senza declutter) 1(corto con declutter) 2(short hail) o 3(medio)
     *  @return    0  in caso di errore   1  in caso di successo  .
     */
    bool test_file(int tipofile);

    void setup_elaborazione(const char* nome_file);
    /**
     *
     *  @brief funzione che elabora il dato radar rimuovendo anaprop e beam blocking
     *  @details  partendo dal livello della mappa dinamica esegue il test di
     *  continuità verticale e laddove non si verifica un cambio di elevazione
     *  corregge il beam blocking per ottenere infine un campo bidimensionale
     *  adeguato alla stima di R che ricopia su tutti i livelli del volume  a
     *  partire dallo 0 fino al livello della mappa dinamica . Memorizza
     *  elevazione finale usata per il campo bidimensionale e l'output della
     *  rimozione della propagazione anomala e quota al livello scelto per la
     *  stima di R.
     */
    void elabora_dato();

    /**
     *  @brief funzione che caratterizza i volumi polari tramite la qualita'
     *  @details utilizza i parametri di ouput dell'elaborazione e la PIA
     *  calcolata qui per calcolare un valore finale di qualita' del dato.
     *  Inoltre calcola il top dell'eco in base a soglia su ogni pixel.
     *  @param  sito identificativo del radar (spc o gat)
     *  @return non ritorna nulla
     */
    void caratterizzo_volume();

    /**
     *  ingresso:dbz in quel punto e attenuazione fin lì
     *  Doviak,Zrnic,1984 for rain as reported in cost 717 final document
     *
     *  @brief   funzione che calcola l'att enuazione totale 
     *  @details Ricevuto in ingresso il dato di Z in byte e l'attenuazione complessiva sul raggio fino al punto in considerazione, calcola l'attenuazione totale 
     *  @param[in]  DBZbyte valore in byte della Z nel pixel
     *  @param[in]  PIA attenuazione totale fino al punto
     *  @return att_tot attenuazione incrementata del contributo nel pixel corrente
     */
    double attenuation(unsigned char DBZbyte, double  PIA);

    /**
     *  @brief   funzione scrittura matrici statistica
     *  @details scrive le statistiche di beam blocking, anaprop, cambio di elevazione in un unsigined char DIM1_ST*DIM1_ST
     */
    void ScrivoStatistica();

    /**
     *
     *  @brief funzione che legge la mappa statica e la mappa di elevazioni da beam blocking e le condensa in un unica mappa
     *  @details crea una mappa di elevazioni scelte in ogni pixel del volume per ottenere successivamante la stima di R
     *  @return non ritorna nulla
     */
    void leggo_first_level();
    /**
     *  conversione da polare a cartesiano alta risoluzione
     *
     *  @brief funzione che crea l'output cartesiano dal polare
     *  @details cicla sui quadranti e su i e j, usando il range e l'azimut ottenuti tramite la funzione creo_matrice_conv()
     *  @return 
     */
    void creo_cart();

    /**
     *  funzioni di conversione cartesiana:associa a pixel matrice alta ris
     *  azimut e range, crea alta risoluzione e crea bassa risoluzione  
     *
     *  @brief funzione che calcola range e azimut su un quadrante centrato in 0,0
     *  @details
     */
    void creo_matrice_conv();

    /**
     *
     *  @brief funzione che  trasforma il dato da cartesiano alta risoluzione a cartesiano bassa risoluzione
     *  @details prende il massimo tra i punti considerati, il passo di ricerca è ZLR_N_ELEMENTARY_PIXEL, cioè il rapporto tra dimensioni ad alta risoluzione e le dimensioni bassa risoluzione.
     *  @return
     */
    void creo_cart_z_lowris();

    /**
     * funzione di scrittura matrici output binario
     *
     * @brief funzione che scrive in un file di output una matrice di byte di dimensione size
     * @details il formato del nome e' $dir/aaammgghhmm_$ext 
     * @param[in] ext estensione file output (aaammgghhmm_$ext)
     * @param[in] content contenuto del file 
     * @param[in] dir directory dove scrivere il file 
     * @param[in] size dimensione della matrice
     * @param[in] matrice matrice di dati da scrivere
     */
    void scrivo_out_file_bin(const char *ext, const char *content, const char *dir,size_t size, const void  *matrice);

    /**
     *  @brief funzione che legge la quota del centro fascio e del limite inferiore del fascio da file
     *  @details legge la quota del centro fascio e del limite inferiore del fascio da file e li memorizza nei vettori hray_inf e hray
     *  @return 
     */
    void leggo_hray();

    /**
     *  lettura dem e quote fascio 
     *
     *  @brief funzione che legge il dem
     *  @details legge il dem e lo memorizza nella matrice dem
     *  @return
     */
    void leggo_dem();

    /**
     * ingresso: elevazione e k di range bin
     *
     * @brief funzione  che calcola la quota in metri del centro del fascio
     * @details distanza=k*dimensionecella +semidimensionecella in metri .quota=f(distinkm, rstinkm, elevazinrad) in metri 
     * @param[in] elevaz elevazione
     * @param[in] k distanza in n0 bin
     * @return q_st quota standard
     */
    float quota_f(float elevaz, int k);
    void conversione_convettiva();
    bool esegui_tutto(const char* nome_file, int file_type);
// added function to calculate beamblocking correction
//
    float BeamBlockingCorrection(unsigned char bin_val, unsigned char beamblocking);

    // RtoDBZ calcolato su aMP e bMP
    float RtoDBZ(float rain) const;


    /**
     *  @brief funzione che a partire dal tempo in secondi arrotonda al NMIN-esimo minuto precedente o successivo
     *  @details 
     *  @param  time intero che rappresenta il numero di secondi
     *  @return restituisce l'arrotondamento in secondi al NMIN-esimo minuto precedente o successivo , se fallisce -1
     *
     * FIXME: this is a legacy from old procedures, code using it should be
     * reviewed to see if it is still needed.
     */
    time_t NormalizzoData(time_t time);

    void StampoFlag();
};

struct CilindricalColumn
{
    std::vector<float> data;

    CilindricalColumn(unsigned zsize, float def_val=-20)
    {
        data.resize(zsize, def_val);
    }

    float& operator[](unsigned z)
    {
        if (z >= data.size()) throw std::runtime_error("cil: fuori coordinata z");
        return data[z];
    }
};

struct CilindricalSlice
{
    std::vector<CilindricalColumn> columns;

    CilindricalSlice(unsigned xsize, unsigned zsize, float def_val=-20)
    {
        columns.reserve(xsize);
        for (unsigned i = 0; i < xsize; ++i)
            columns.push_back(CilindricalColumn(zsize, def_val));
    }

    CilindricalColumn& operator[](unsigned x)
    {
        if (x >= columns.size()) throw std::runtime_error("cil: fuori coordinata x");
        return columns[x];
    }
};

struct CilindricalVolume
{
    std::vector<CilindricalSlice> slices;

    void allocate(unsigned xsize, unsigned zsize)
    {
        slices.reserve(NUM_AZ_X_PPI);
        for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
            slices.push_back(CilindricalSlice(xsize, zsize));
    }

    CilindricalSlice& operator[](unsigned i)
    {
        if (i >= slices.size()) throw std::runtime_error("slices: fuori coordinata i");
        return slices[i];
    }
};

struct CalcoloVPR
{
    log4c_category_t* logging_category;

    CUM_BAC& cum_bac;
    long int area_vpr[NMAXLAYER]; /*area degli strati*/
    // ricampionamento del volume in coordinate cilindriche
    CilindricalVolume cil;
    long int gap; /* distanza temporale dall'ultimo file vpr */
    float zeroterm;//zerotermico
    float t_ground;
    //matrici che dicono se pixel convettivo secondo VIZ, STEINER, riassuntiva mette +50
    unsigned char *conv_VIZ[NUM_AZ_X_PPI],*conv_STEINER[NUM_AZ_X_PPI],*conv[NUM_AZ_X_PPI];
    unsigned char stratiform[NUM_AZ_X_PPI][MAX_BIN];
    float vpr[NMAXLAYER];/* vpr */
    int hvprmax; /* quota picco vpr */
    //elab classificazione: lista punti convettivi, iaz e ira, le dimensioni sono le massime possibili, in realtà i punti sono molti meno
    int lista_conv[NUM_AZ_X_PPI*MAX_BIN][2];
    //lista punti di background iaz e ira, le dimensioni sono le massime possibili, in realtà i punti sono molti meno
    int lista_bckg[NUM_AZ_X_PPI*MAX_BIN][2];
    // array di parametri, fisso , RES_HOR_CIL E RES_VERT_CIL
    float resol[2];
    int heating,livmin; /* variabile di riscaldamento e quota livello minimo calcolato*/
    int x_size,z_size;
    long int ncv,np;
    float *convective_radius;
    float htbb, hbbb;
    // array contenenti Z di background
    double *Z_bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in mm^6/m^3
    float *bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in dB
    unsigned char corr_polar[NUM_AZ_X_PPI][MAX_BIN];/*correzione vpr in byte 0-128 negativa 128-256 positiva, in coord az-ra*/
    unsigned char neve[NUM_AZ_X_PPI][MAX_BIN];/* matrice az-range che memorizza punti di neve*/
    int ier_vpr, ier_comb,ier_max,ier_stampa_vpr;/* flag d'errore su calcolo vpr istantaneo, combinazione vpr, funzione get_t_ground */
    // dati per vpr
    VolumeInfo<unsigned char>* flag_vpr; // punti del volume polare ok per calcolo VPR*/
    //obsol.
    float stdev;// obsol.
    // files vpr
    FILE *test_vpr;

    CalcoloVPR(CUM_BAC& cum_bac);
    ~CalcoloVPR();

    /**
     *  @brief   funzione che analizza il profilo 
     *  @details analizza il profilo usando : la temperatura al suolo, la quota del massimo, e una funzione di interpolazione 
     *  @param[out]  vpr_liq valore del profilo al livello liquido
     *  @param[out]  snow matrice che indica se è presente neve o no secondo l'analisi fatta
     *  @param[out]  hliq quota del livello liquido
     *  @return ier_ana valore che indica se tutto è andato a buon fine (0) o no (1)
     */
    int analyse_VPR(float *vpr_liq,int *snow,float *hliq);

    /**
     *  calcola riscaldamento in quarti d'ora
     *  @brief funzione che calcola quanto il profilo è 'caldo' restituendo la variabile heating
     *  @details calcola il 'riscaldamento' del profilo (numero di combinazioni successive) e stampa il file che contiene questo valore oltre al file contenente l'ultima data se il calcolo del profilo è andato ok. Se il numero è superiore a WARM passa direttamente al valore MEMORY che rappresenta la 'memoria del profilo'.
     *   heating=heating-gap; se il profilo non è stato aggiornato 
     *   heating=heating-gap+2; se il profilo è stato aggiornato 
     *   heating=MEMORY;   se heating raggiunge WARM resta costante finchè non inizia raffreddamento 
     *  @return heating , numero di combinazioni di riscaldamento
     */
    int profile_heating();

    /**
     *  trova il massimo del profilo
     *  @brief funzione che trova la quota del massimo valore del profilo
     *  @details trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca come massimo di almeno 5 dBZ più alto in quota
     *  @param[in] sito radar su cui calcolare profilo
     *  @return 0 se ok 1 se fallisce
     */
    int trovo_hvprmax(int *hmax);

    /**
     *  crea vpr istantaneo
     *  @brief funzione che calcola il profilo istantaneo  secondo il metodo di Germann e Joss (2003)
     *  @details   calcola il VPR istantaneo secondo il metodo di Germann e Joss (2003) 
     *  Per il calcolo si considerano i punti con Z>THR_VPR, qualità>QMIN_VPR, BeamBlocking<20 percento e clutter free all'interno del volume scelto.
     *  Il profilo è poi soggetto a quality check e eventualmente viene rigettato (return(1)) 
     *  @param[out] cv volume precipitante
     *  @param[out] ct volume totale 
     *  @param[out] vpr1 vettore vpr istantaneo
     *  @param[out] area_vpr vettore area di ogni stato 
     *  @param[in]  sito  sito radar
     *  @return  0 se ok 1 se fallisce
     */ 
    int func_vpr(long int *cv, long int *ct, float vpr1[], long int area_vpr[]);

    /**
     *
     *  @brief funzione che combina il profilo verticale corrente con quello precedente tramite il metodo di Germann
     *  @details  oltre a lanciare il calcolo del profilo istantaneo provvede alla combinazione del profilo calcolato con il precedente calcolato entro  un limite massimo di distanza temporale pari a 10 quarti d'ora.  restituisce un codice integer pari a 0 se ok 1 se fallisce 
     *  @param[in] sito radar corrente 
     *  @return 0 se combinazione ok 1 se fallisce
     */
    int combina_profili();

    /**
     *  calcola valore di background per individuare pixel convettivo
     *
     *  @brief funzione  che calcola il background 
     *  @details la classificazione di Steiner non ha bisogno di ricampionamento cilindrco perciò  uso direttamente la matrice polare
     */
    void calcolo_background();

    /**
     *
     *  @brief funzione  che classifica la precipitazione se stratiforme o convettiva
     *  @details esegue anche il ricampionamento cilindrico al suo interno
     * @return 
     */
    void classifica_rain();

    /**
     *  classifica tramite Vertical Integrated Z
     *  @brief funzione  che classifica secondo il metodo VIZ
     *  @details calcolo per ogni pixel polare l'integrale verticale esclusa la fascia della bright band
     */
    void classifico_VIZ();

    /**
     *  @brief funzione  che classifica secondo STEINER
     *  @details segna come convettivi i punti che hanno valore superiore a 40 dBZ e differenza col background elevata, quindi ingrandisce i nuclei di un raggio variabile
     */
    void classifico_STEINER();

    /**
     *  correzione vpr
     *  @brief funzione che corregge per il profilo verticale
     *  @details ciclando su tutti i bins della cartesiana polare scelta per la stima della pioggia,
     * - trovo la quota del centro del pixel
     * - correggo se: hbin>hliq , valore maggiore di soglia correzione (0 dBZ)
     * - se sefinita CLASS e pixel convettivo non correggo 
     * @param[in] sito radar su cui calcolare profilo
     * @return 0 se ok 1 se fallisce
     */
    int corr_vpr();

    /**
     *  @brief funzione  che ingrandisce i nuclei di Steiner
     *  @details ingrandisce i nuclei di Steiner di un valore pari al raggio convettivo
     *  @param[in] cr raggio convettivo
     *  @param[in] ja indice di azimut
     *  @param[in] kr indice di range
     */ 
    void ingrasso_nuclei(float cr,int ja,int kr);

    /**
     *  fa il merge dei metodi
     *  @brief funzione  che interseca i punti convettivi delle due classificazioni Viz e Steiner e sottrae quelli con  picco stratiforme
     *  @return non ritorna valori
     */
    void merge_metodi();

    // stampa profilo combinato
    int stampa_vpr();

    /**
     *  trova lo zero termico
     *  @brief funzione che legge il file contenente la quota dello 0 termico 
     *  @return ier codice di errore apertura file
     */
    int trovo0term();

    void esegui_tutto();
};

// Utility functions

/**
 *  combina livelli
 *
 *  @brief funzione che compone i singoli livelli del profilo v0 e v1 
 *  @details  result=((1.-peso)*v0+peso*v1)
 *  @param[in]  v0 valore del profilo vecchio nel punto 
 *  @param[in]  v1 valore del profilo nuovo nel punto 
 *  @param[in]  nodata valore dei nodata
 *  @param[in]  peso peso del profilo nuovo
 *  @return result :ritorna il valore combinato dei due profili e se uno dei due manca mette nodata
*/ 
float comp_levels(float v0, float v1, float nodata, float peso);

/**
 * funzione di controllo esistenza file e restituisce puntatore a file+messaggio ev errore o info su status apertura, utile? 
 *
 *  @brief funzione che restituisce un puntatore a file in lettura o scrittura dopo aver controllato esistenza e permessi 
 *  @details scrive un messaggio sul log in caso di errore durante l'apertura e l'accesso col permesso richiesto ed esce
 *  @param[in]  nome_file nome del file da aprire
 *  @param[in]  content contenuto del file (stringa esplicativa)
 *  @param[in]  mode modalita' di apertura (scrittura , lettura)
 *  @return file ritorna un puntatore a file
 */
FILE *controllo_apertura(const char *nome_file, const char *content, const char *mode);

}

#endif
