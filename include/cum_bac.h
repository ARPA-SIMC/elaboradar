#ifndef ARCHIVIATORE_CUM_BAC_CLASS_H
#define ARCHIVIATORE_CUM_BAC_CLASS_H

#include "logging.h"
#include "assets.h"
#include <cum_bac_SP20_BB_VPR.h>


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
//Array delle elevazioni
static const int elev_array_spc[NEL]={6,15,26,36,46,57,80,108,148,159,170,180,190,200,210};//GLI ULTIMI 5 fittizi: ANNA 30-03-2011
static const int elev_array_gat[NEL]={6,16,27,37,45,55,65,76,85,95,105,126,149,174,201};//105,126,149,174,201 è da completare NEL=15:ANNA 30-03-2011
#endif
// v. parametri SHORT
#ifdef MEDIUM
#define NEL 5
#define RANGE_KM                 250
#define CART_DIM_ZLR             512
#define ZLR_N_ELEMENTARY_PIXEL   1
#define ZLR_OFFSET               CART_DIM_ZLR/2
static const int elev_array_spc[NEL]={6,16,26,36,47};//ANNA 30-03-2011
static const int elev_array_gat[NEL]={6,16,27,36,47};//ANNA 30-03-2011
#define NMIN 2
#define MAX_TIME_DIFF 1
#endif

// dimensioni cella a seconda del tipo di acquisizione
static const float size_cell[]={62.5,125.,250.,500.,1000.,2000.};

//Dimensioni matrici statistica
#define DIM1_ST 16
#define DIM2_ST 13/*Cambiata dimensione a 13 per cambio dimensione raggio radar*/

extern int elev_array[NEL];

class Volume
{
public:
    time_t acq_date;

    Volume()
        : acq_date(0)
    {
    }
};

class CUM_BAC
{
public:
    log4c_category_t* logging_category;

    Assets assets;

    /// Feature set required for this run
    bool do_quality;
    bool do_beamblocking;
    bool do_declutter;
    bool do_bloccorr;
    bool do_vpr;
    bool do_class;

    Volume volume;

    int MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT;/* differenza massima tra le due elevazioni successive perchè non sia clutter e valore minimo a quella superiore pe il primo e per i successivi (NEXT) bins*/


    /*-----------------------------------------------------------
      Variabili globali
      ------------------------------------------------------------*/

    T_MDB_ap_beam_header  old_beam_header;  /* queste due dichiarazioni le metto qui perchè sono sfaticato*/
    T_MDB_data_header   old_data_header;

    char errori[256];

    //dato di base volume polare, struttura definita in libSP20
    struct VOL_POL vol_pol[NEL][NUM_AZ_X_PPI];

    //numero raggi per elevazione
    int nbeam_elev[NEL];

    //matrici per ricampionamento cartesiano
    float azimut[MAX_BIN][MAX_BIN];
    float range[MAX_BIN][MAX_BIN];
    unsigned char cart[MAX_BIN*2][MAX_BIN*2];
    double  cartm[MAX_BIN*2][MAX_BIN*2];  /* Z media dei bins adiacenti al punto */
    //se definita Z_LOWRIS, Z cartesiana al livello più basso
    unsigned char z_out[CART_DIM_ZLR][CART_DIM_ZLR];

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
    float zeroterm;//zerotermico

    // attenuazione in formato cartesiano max risoluzione
    unsigned char att_cart[NUM_AZ_X_PPI][MAX_BIN]; /* matrice azimut-range di attenuazione */
    //quota centro fascio polare, cartesiana max risoluzione e cartesiana 1x1
    unsigned short quota_rel[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio relativa al suolo in prop da rsd e elevazioni nominali, in coordinate azimut range*/
    unsigned short quota[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio in prop standard e elev reali in coordinate azimut range*/
    unsigned short quota_cart[MAX_BIN*2][MAX_BIN*2];/*quota fascio in coordinate cart 1024*1024, risoluzione minima*/
    unsigned char quota_1x1[CART_DIM_ZLR][CART_DIM_ZLR];/* quota in formato 256*256 in centinaia di metri, risoluzione ZLR */
    //beam blocking cartesiano max resol e 1x1
    unsigned char beam_blocking_xy[MAX_BIN*2][MAX_BIN*2]; //beamblocking cartesiano max resol
    unsigned char beam_blocking_1x1[CART_DIM_ZLR][CART_DIM_ZLR];//beam blocking cartesiano 1x1
    //uscite anaprop
    unsigned char dato_corrotto[NUM_AZ_X_PPI][MAX_BIN]; /*uscita controllo anaprop in coordinate azimut range */
    unsigned char dato_corr_xy[MAX_BIN*2][MAX_BIN*2]; //uscite anap  cartesiano max resol
    unsigned char dato_corr_1x1[CART_DIM_ZLR][CART_DIM_ZLR]; //uscite anap cartesiano  1x1
    //elevazioni finali come sopra
    unsigned char elev_fin[NUM_AZ_X_PPI][MAX_BIN]; /* elevazione finale in coordinate azimut range  */
    unsigned char elev_fin_xy[MAX_BIN*2][MAX_BIN*2];
    unsigned char elev_fin_1x1[CART_DIM_ZLR][CART_DIM_ZLR];
    // metrici qualita' come sopra
    unsigned char qual[NEL][NUM_AZ_X_PPI][MAX_BIN]; /* qualita volume polare */
    unsigned char qual_Z_cart[MAX_BIN*2][MAX_BIN*2]; /* qualita della Z in formato 1024*1024, risoluzione minima */
    unsigned char qual_Z_1x1[CART_DIM_ZLR][CART_DIM_ZLR];/* qualita della Z in formato 256*256, risoluzione ZLR */
    // top, come sopra
    unsigned char top[NUM_AZ_X_PPI][MAX_BIN];
    unsigned char topxy[MAX_BIN*2][MAX_BIN*2];
    unsigned char  top_1x1[CART_DIM_ZLR][CART_DIM_ZLR];

    // uscite  vpr: correzione VPR , come sopra
    unsigned char  corr_polar[NUM_AZ_X_PPI][MAX_BIN];/*correzione vpr in byte 0-128 negativa 128-256 positiva, in coord az-ra*/
    unsigned char  corr_cart[MAX_BIN*2][MAX_BIN*2];
    unsigned char corr_1x1[CART_DIM_ZLR][CART_DIM_ZLR];
    // uscite vpr: neve, come sopra
    unsigned char neve[NUM_AZ_X_PPI][MAX_BIN];/* matrice az-range che memorizza punti di neve*/
    unsigned char neve_cart[MAX_BIN*2][MAX_BIN*2];/* neve formato 1024*1024, risoluzione minima */
    unsigned char neve_1x1[CART_DIM_ZLR][CART_DIM_ZLR];/* neve in formato 256*256, risoluzione ZLR */
    // dati per vpr
    unsigned char flag_vpr[NEL][NUM_AZ_X_PPI][MAX_BIN];/* punti del volume polare ok per calcolo VPR*/
    float vpr[NMAXLAYER];/* vpr */
    long int gap; /* distanza temporale dall'ultimo file vpr */
    long int area_vpr[NMAXLAYER]; /*area degli strati*/
    int ier_vpr, ier_comb,ier_max,ier_stampa_vpr;/* flag d'errore su calcolo vpr istantaneo, combinazione vpr, funzione get_t_ground */
    int hvprmax; /* quota picco vpr */
    int heating,livmin; /* variabile di riscaldamento e quota livello minimo calcolato*/
    float t_ground;
    // dati elab vpr
    float chisqfin; //???puo' essere def in anal
    float rmsefin;
    // files vpr
    FILE *test_vpr;
    //obsol.
    float stdev;// obsol.

    //matrici per classificazione: cappi
    unsigned char cappi[NUM_AZ_X_PPI][MAX_BIN];
    //matrici che dicono se pixel convettivo secondo VIZ, STEINER, riassuntiva mette +50
    unsigned char *conv_VIZ[NUM_AZ_X_PPI],*conv_STEINER[NUM_AZ_X_PPI],*conv[NUM_AZ_X_PPI];
    // uscite: matrici class max resol e 1x1
    unsigned char conv_cart[MAX_BIN*2][MAX_BIN*2], conv_1x1[CART_DIM_ZLR][CART_DIM_ZLR];
    //uscite:matrici cappi max resol e 1x1
    unsigned char cappi_1x1[CART_DIM_ZLR][CART_DIM_ZLR],cappi_cart[MAX_BIN*2][MAX_BIN*2],stratiform[NUM_AZ_X_PPI][MAX_BIN];
    //elab classificazione: lista punti convettivi, iaz e ira, le dimensioni sono le massime possibili, in realtà i punti sono molti meno
    int lista_conv[NUM_AZ_X_PPI*MAX_BIN][2];
    //lista punti di background iaz e ira, le dimensioni sono le massime possibili, in realtà i punti sono molti meno
    int lista_bckg[NUM_AZ_X_PPI*MAX_BIN][2];
    // array contenenti Z di background
    double *Z_bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in mm^6/m^3
    float *bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in dB
    // ricampionamento del volume in coordinate cilindriche
    float **cil[NUM_AZ_X_PPI];
    // array di parametri, fisso , RES_HOR_CIL E RES_VERT_CIL
    float resol[2];
    int x_size,z_size;
    long int ncv,ncs,np;
    float *convective_radius;
    float htbb, hbbb;

    /* variabili tolte perchè non presenti nel codice cum_bac... controllare che non richiamino qualcosa nelle funzioni
       struct tm *time_dbp;
       T_time, T_data, T_ora..*/


    CUM_BAC();

    bool read_sp20_volume(const char* nome_file, const char* sito, int file_type);
    bool read_odim_volume(const char* nome_file, const char* sito, int file_type);
    bool test_file(int tipofile);
    void setup_elaborazione(const char* nome_file, const char* sito);
    int elabora_dato();
    void caratterizzo_volume();
    /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/
    double attenuation(unsigned char DBZbyte, double  PIA);
    void ScrivoStatistica();
    void leggo_first_level();
    void creo_cart();
    void creo_matrice_conv();
    void creo_cart_z_lowris();
    void scrivo_out_file_bin(const char *ext,char *content,char *dir,size_t size, void  *matrice);
    FILE *controllo_apertura(const char *nome_file, const char *content, const char *mode);
    void leggo_hray();
    void leggo_dem();
    int func_vpr(long int *cv, long int *ct, float vpr1[], long int area_vpr[], const char *sito);
    float comp_levels(float v0, float v1, float nodata, float peso);
    int combina_profili(const char *sito);
    int profile_heating();
    int stampa_vpr();
    int corr_vpr(const char *sito);
    int trovo_hvprmax(int *hmax);
    int analyse_VPR(float *vpr_liq,int *snow,float *hliq, const char *sito);
    int interpola_VPR(float a[], int ma);
    int testfit(float a[], float chisq, float chisqin);
    float quota_f(float elevaz, int k);
    void classifica_rain();
    void classifico_VIZ();
    void classifico_STEINER();
    void calcolo_background();
    void ingrasso_nuclei(float cr,int ja,int kr);
    void merge_metodi();
    int trovo0term();
    bool esegui_tutto(const char* nome_file, int file_type, const char* sito);
// added function to calculate beamblocking correction
//
    float BeamBlockingCorrection(unsigned char bin_val, unsigned char beamblocking);

};

// Utility functions

/**
 * Check if time is inbetween a ??? interval, and return it rounded to 5-minute
 * intervals
 *
 * FIXME: this is a legacy from old procedures, code using it should be
 * reviewed to see if it is still needed.
 */
time_t NormalizzoData(time_t time);

/// Linear gauss
void lineargauss(float x, float a[], float *y, float dyda[],int na);

#endif
