#ifndef ARCHIVIATORE_CUM_BAC_CLASS_H
#define ARCHIVIATORE_CUM_BAC_CLASS_H

#include "logging.h"
#include "assets.h"
#include <cum_bac_SP20_BB_VPR.h>
#include "volume.h"

#ifdef __cplusplus
extern "C" {
#endif
// libreria radar
#include <func_SP20read.h>
#ifdef __cplusplus
}
#endif


//algoritmo
#include <MP_par.h>
#include <vpr_par.h>
#include <geo_par.h>


// Definizioni e variabili che cambiano da corto a medio
#ifdef NEL // Se NEL è definito lo undefinisco perchè è ridefinito sotto
#undef NEL
#endif

#ifdef SHORT
//Risoluzioni e limiti spaziali
#define NEL 15                // n0 elevazioni massimo
#define RANGE_KM                 125  // distanza massima in km
#define CART_DIM_ZLR             256 // dimensione matrice a 1x1 km
#define ZLR_N_ELEMENTARY_PIXEL   4
#define ZLR_OFFSET               0
//Parametri passare da minuti del file a minuti standard arrotondando per difetto o eccesso ( prima si arrotondava al 5° ora si arrotonda al minuto )
#define NMIN 1 // cambiato da #define NMIN 5 a 1 dopo inserimento minuti maltempo, step di arrotondamento in minuti
#define MAX_TIME_DIFF 3 // massima differenza in minuti tra data acquisizione e standard per arrotondare per difetto
#endif
// v. parametri SHORT
#ifdef MEDIUM
#define NEL 5
#define RANGE_KM                 250
#define CART_DIM_ZLR             512
#define ZLR_N_ELEMENTARY_PIXEL   1
#define ZLR_OFFSET               CART_DIM_ZLR/2
#define NMIN 2
#define MAX_TIME_DIFF 1
#endif

//Dimensioni matrici statistica
#define DIM1_ST 16
#define DIM2_ST 13/*Cambiata dimensione a 13 per cambio dimensione raggio radar*/

namespace cumbac {

struct Site;

template<typename T, unsigned SX, unsigned SY=SX>
struct Image
{
    T data[SY][SX];

    Image()
    {
        for (unsigned y = 0; y < SY; ++y)
            for (unsigned x = 0; x < SX; ++x)
                data[y][x] = 0;
    }

    T* operator[](unsigned y) { return data[y]; }
    const T* operator[](unsigned y) const { return data[y]; }

    T min() const
    {
        T res = data[0][0];
        for (unsigned y = 0; y < SY; ++y)
            for (unsigned x = 0; x < SX; ++x)
                if (data[y][x] < res)
                    res = data[y][x];
        return res;
    }

    T max() const
    {
        T res = data[0][0];
        for (unsigned y = 0; y < SY; ++y)
            for (unsigned x = 0; x < SX; ++x)
                if (data[y][x] > res)
                    res = data[y][x];
        return res;
    }

    T avg() const
    {
        double res = 0;
        for (unsigned y = 0; y < SY; ++y)
            for (unsigned x = 0; x < SX; ++x)
                res += (double)data[y][x] / (double)(SX * SY);
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
    bool do_quality;
    bool do_beamblocking;
    bool do_declutter;
    bool do_bloccorr;
    bool do_vpr;
    bool do_class;
    bool do_zlr_media;

    cumbac::Volume volume;

    CalcoloVPR* calcolo_vpr;

    int MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT;/* differenza massima tra le due elevazioni successive perchè non sia clutter e valore minimo a quella superiore pe il primo e per i successivi (NEXT) bins*/


    /*-----------------------------------------------------------
      Variabili globali
      ------------------------------------------------------------*/

    //matrici per ricampionamento cartesiano
    float azimut[MAX_BIN][MAX_BIN];
    float range[MAX_BIN][MAX_BIN];
    Image<unsigned char, MAX_BIN*2> cart;
    Image<double, MAX_BIN*2> cartm;  /* Z media dei bins adiacenti al punto */
    //se definita Z_LOWRIS, Z cartesiana al livello più basso
    Image<unsigned char, CART_DIM_ZLR> z_out;

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
    unsigned char first_level[NUM_AZ_X_PPI][MAX_BIN]; //mappa dinamica complessiva
    unsigned char first_level_static[NUM_AZ_X_PPI][MAX_BIN];//mappa statica

    unsigned char bb_first_level[NUM_AZ_X_PPI][MAX_BIN];  /* mappa di elevazioni da beam blocking (input)*/
    unsigned char beam_blocking [NUM_AZ_X_PPI][MAX_BIN];   /* mappa di beam blocking (input)*/

    //variabili legate a propagazione e beam blocking, da prog_bb
    float hray[MAX_BIN][NEL];  /*quota centro fascio in funzione della distanza e elevazione*/
    float hray_inf[MAX_BIN][NEL]; /*quota limite inferiore fascio in funzione della distanza e elevazione*/
    float dem[NUM_AZ_X_PPI][MAX_BIN]; /*dem in coordinate azimut range*/
    float dtrs;// distanza temporale radiosondaggio

    // attenuazione in formato cartesiano max risoluzione
    unsigned char att_cart[NUM_AZ_X_PPI][MAX_BIN]; /* matrice azimut-range di attenuazione */
    //quota centro fascio polare, cartesiana max risoluzione e cartesiana 1x1
    unsigned short quota_rel[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio relativa al suolo in prop da rsd e elevazioni nominali, in coordinate azimut range*/
    unsigned short quota[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio in prop standard e elev reali in coordinate azimut range*/
    Image<unsigned short, MAX_BIN*2> quota_cart;/*quota fascio in coordinate cart 1024*1024, risoluzione minima*/
    Image<unsigned char, CART_DIM_ZLR> quota_1x1;/* quota in formato 256*256 in centinaia di metri, risoluzione ZLR */
    //beam blocking cartesiano max resol e 1x1
    Image<unsigned char, MAX_BIN*2> beam_blocking_xy; //beamblocking cartesiano max resol
    Image<unsigned char, CART_DIM_ZLR> beam_blocking_1x1;//beam blocking cartesiano 1x1
    //uscite anaprop
    unsigned char dato_corrotto[NUM_AZ_X_PPI][MAX_BIN]; /*uscita controllo anaprop in coordinate azimut range */
    Image<unsigned char, MAX_BIN*2> dato_corr_xy; //uscite anap  cartesiano max resol
    Image<unsigned char, CART_DIM_ZLR> dato_corr_1x1; //uscite anap cartesiano  1x1
    Image<unsigned char, MAX_BIN*2> elev_fin_xy;
    Image<unsigned char, CART_DIM_ZLR> elev_fin_1x1;
    // metrici qualita' come sopra
    unsigned char qual[NEL][NUM_AZ_X_PPI][MAX_BIN]; /* qualita volume polare */
    Image<unsigned char, MAX_BIN*2> qual_Z_cart; /* qualita della Z in formato 1024*1024, risoluzione minima */
    Image<unsigned char, CART_DIM_ZLR> qual_Z_1x1;/* qualita della Z in formato 256*256, risoluzione ZLR */
    // top, come sopra
    unsigned char top[NUM_AZ_X_PPI][MAX_BIN];
    Image<unsigned char, MAX_BIN*2> topxy;
    Image<unsigned char, CART_DIM_ZLR> top_1x1;

    // uscite  vpr: correzione VPR , come sopra
    Image<unsigned char, MAX_BIN*2> corr_cart;
    Image<unsigned char, CART_DIM_ZLR> corr_1x1;
    // uscite vpr: neve, come sopra
    Image<unsigned char, MAX_BIN*2> neve_cart;/* neve formato 1024*1024, risoluzione minima */
    Image<unsigned char, CART_DIM_ZLR> neve_1x1;/* neve in formato 256*256, risoluzione ZLR */

    //matrici per classificazione: cappi
    unsigned char cappi[NUM_AZ_X_PPI][MAX_BIN];
    // uscite: matrici class max resol e 1x1
    Image<unsigned char, MAX_BIN*2> conv_cart;
    Image<unsigned char, CART_DIM_ZLR> conv_1x1;
    //uscite:matrici cappi max resol e 1x1
    unsigned char cappi_1x1[CART_DIM_ZLR][CART_DIM_ZLR],cappi_cart[MAX_BIN*2][MAX_BIN*2];

    /* variabili tolte perchè non presenti nel codice cum_bac... controllare che non richiamino qualcosa nelle funzioni
       struct tm *time_dbp;
       T_time, T_data, T_ora..*/


    CUM_BAC(const char* site_name);
    ~CUM_BAC();

    bool read_sp20_volume(const char* nome_file, int file_type);
    bool read_odim_volume(const char* nome_file, int file_type);
    bool test_file(int tipofile);
    void setup_elaborazione(const char* nome_file);
    int elabora_dato();
    void caratterizzo_volume();
    /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/
    double attenuation(unsigned char DBZbyte, double  PIA);
    void ScrivoStatistica();
    void leggo_first_level();
    void creo_cart();
    void creo_matrice_conv();
    void creo_cart_z_lowris();
    void scrivo_out_file_bin(const char *ext, const char *content, const char *dir,size_t size, const void  *matrice);
    void leggo_hray();
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
    void class_conv_fixme_find_a_name();
    bool esegui_tutto(const char* nome_file, int file_type);
// added function to calculate beamblocking correction
//
    float BeamBlockingCorrection(unsigned char bin_val, unsigned char beamblocking);

    // RtoDBZ calcolato su aMP e bMP
    float RtoDBZ(float rain) const;
};

struct CalcoloVPR
{
    log4c_category_t* logging_category;

    CUM_BAC& cum_bac;
    long int area_vpr[NMAXLAYER]; /*area degli strati*/
    // ricampionamento del volume in coordinate cilindriche
    float **cil[NUM_AZ_X_PPI];
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
    long int ncv,ncs,np;
    float *convective_radius;
    float htbb, hbbb;
    // array contenenti Z di background
    double *Z_bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in mm^6/m^3
    float *bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in dB
    unsigned char corr_polar[NUM_AZ_X_PPI][MAX_BIN];/*correzione vpr in byte 0-128 negativa 128-256 positiva, in coord az-ra*/
    unsigned char neve[NUM_AZ_X_PPI][MAX_BIN];/* matrice az-range che memorizza punti di neve*/
    int ier_vpr, ier_comb,ier_max,ier_stampa_vpr;/* flag d'errore su calcolo vpr istantaneo, combinazione vpr, funzione get_t_ground */
    // dati elab vpr
    float chisqfin; //???puo' essere def in anal
    float rmsefin;
    // dati per vpr
    unsigned char flag_vpr[NEL][NUM_AZ_X_PPI][MAX_BIN];/* punti del volume polare ok per calcolo VPR*/
    //obsol.
    float stdev;// obsol.
    // files vpr
    FILE *test_vpr;

    CalcoloVPR(CUM_BAC& cum_bac);
    int analyse_VPR(float *vpr_liq,int *snow,float *hliq);
    int profile_heating();
    int trovo_hvprmax(int *hmax);
    int func_vpr(long int *cv, long int *ct, float vpr1[], long int area_vpr[]);
    int combina_profili();
    void calcolo_background();
    /**
     *
     *  @brief funzione  che classifica la precipitazione se stratiforme o convettiva
     *  @details esegue anche il ricampionamento cilindrico al suo interno
     * @return 
     */
    void classifica_rain();
    void classifico_VIZ();
    void classifico_STEINER();
    int interpola_VPR(float a[], int ma);
    int corr_vpr();
    void ingrasso_nuclei(float cr,int ja,int kr);
    void merge_metodi();
    int stampa_vpr();
    int trovo0term();

    void esegui_tutto();
};

// Utility functions

/// Linear gauss
void lineargauss(float x, float a[], float *y, float dyda[],int na);

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
 *
 *  @brief funzione che restituisce un puntatore a file in lettura o scrittura dopo aver controllato esistenza e permessi 
 *  @details scrive un messaggio sul log in caso di errore durante l'apertura e l'accesso col permesso richiesto ed esce
 *  @param[in]  nome_file nome del file da aprire
 *  @param[in]  content contenuto del file (stringa esplicativa)
 *  @param[in]  mode modalita' di apertura (scrittura , lettura)
 *  @return file ritorna un puntatore a file
 */
FILE *controllo_apertura(const char *nome_file, const char *content, const char *mode);


/**
 *  testa i parametri del fit in modo che abbiano significato fisico 
 *
 *  @brief   funzione che testa il fit dell'interpolazione del profilo
 *  @details verifica che i parametri del fit del profilo abbiano senso
 *  @param[in] a[] vettore dei parametri della funzione
 *  @param[in] chisq  chiquare
 *  @return codice di uscita 0
 *
 */
int testfit(float a[], float chisq, float chisqin);


}

/**
 * Check if time is inbetween a ??? interval, and return it rounded to 5-minute
 * intervals
 *
 * FIXME: this is a legacy from old procedures, code using it should be
 * reviewed to see if it is still needed.
 */
time_t NormalizzoData(time_t time);

#endif
