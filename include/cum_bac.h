#ifndef ARCHIVIATORE_CUM_BAC_CLASS_H
#define ARCHIVIATORE_CUM_BAC_CLASS_H

/**
 *  @file
 *  @ingroup elaboradar
 *  @brief codice principale di elaborazione dei volumi di riflettivita' radar usato per impulso corto
 *  @details questo codice contiene il main e alcune funzioni usate nell'elaborazione della riflettivita' radar 
 *  librerire necessarie: lib_SP20
*/

#include <radarelab/logging.h>
#include <radarelab/volume.h>
#include <radarelab/elev_fin.h>
#include <radarelab/algo/anaprop.h>
#include <radarelab/algo/dbz.h>
#include <radarelab/algo/vpr.h>
#include <radarelab/cylindrical.h>
#include "assets.h"
#include <stdexcept>
#include <cmath>

//algoritmo
#include <radarelab/vpr_par.h>

// TODO: prima o poi arriviamo a far senza di questi define
#define NUM_AZ_X_PPI 400

#define MAX_BIN 512

#define RES_VERT_CIL 0.25
#define RES_HOR_CIL 0.25

/** 
 * namespace per gli algoritmi di elaborazione
 */
namespace radarelab {

namespace algo {
struct CalcoloSteiner;
struct CalcoloVIZ;
}

}

/** 
 * name space generale del programma
 */
namespace elaboradar {

struct Site;

struct CalcoloVPR;
struct Cart;
struct CartLowris;
struct CartProducts;

/**
 *  Classe principale del programma. Deriva ancora dal vecchio nome del programma. Sarà da sistemare in futuro
 */
class CUM_BAC
{
public:
    /**
     * Read from SP20 data file 
     * @param [out] volume  - Container for radar data
     * @param [out] site - Radar site object descriptor
     * @param [in] nome_file - file to be read
     * @param [in] do_clean - flag to abilitate cleaning procedure 
     * @param [in] do_medium - flag to force processing as medium range data  
     */
    static void read_sp20_volume(radarelab::Volume<double>& volume, const Site& site, const char* nome_file, bool do_clean=false, bool do_medium=false);
    /**
     * Read from ODIM data file 
     * @param [out] volume  - Container for radar data
     * @param [out] site - Radar site object descriptor
     * @param [in] nome_file - file to be read
     * @param [in] do_clean - flag to abilitate cleaning procedure 
     * @param [in] do_medium - flag to force processing as medium range data  
     */
  static void read_odim_volume(radarelab::Volume<double>& volume, const Site& site, const char* nome_file, char* fuzzypath, bool do_clean=false, bool do_medium=false, bool set_undetect=false);

    log4c_category_t* logging_category; ///< logging category 

    unsigned MyMAX_BIN;         ///< maximum number of beam size
    const Config& cfg;          ///< Configuration object
    const Site& site;           ///< site information object
    Assets assets;          ///< others

    bool do_medium;         ///< medium processing flag

    /// Feature set required for this run
    bool do_quality = false;        ///< quality flag
    bool do_beamblocking = false;   ///< beamblocking corretion 
    bool do_bloccorr = false;       ///< bloccorrection
    bool do_declutter = false;      ///< use only static declutter map
    bool do_class = false;      ///< Convective-stratiform classification 
    bool do_zlr_media = false;      ///< Compute ZLR map using averaging
    bool do_devel = false;      ///< Produce additional output
    bool do_readStaticMap = false;  ///< Read Static clutter map
    bool do_anaprop=false;      ///< anaprop correction

    // add set_undetect
  bool set_undetect=false;///Set to Z undetect value the Zpixels classified as non-meteo echoes

    radarelab::Volume<double>& volume;      ///< Polar volume of Reflectivity 
    radarelab::Volume<double> SD_Z6;        ///< Polar volume of standard deviation of reflectivity over 6 km length
    radarelab::CylindricalVolume cil;       ///< Volume resampled as a cylindrical volume

    radarelab::algo::DBZ dbz;           ///< ????
    radarelab::Volume<unsigned char> flag_vpr;

    CalcoloVPR* calcolo_vpr = 0;    ///< Oggetto per calcolare e correggere con VPR

    /*-----------------------------------------------------------
      Variabili globali
      ------------------------------------------------------------*/

    // Data del volume che abbiamo letto
    char date[20];          ///<  Acquisition date 

    //matrici first_level e first level da beam blocking e valore beam blocking
    radarelab::PolarScan<unsigned char> first_level;        ///< mappa dinamica complessiva
    radarelab::PolarScan<unsigned char> first_level_static; ///< mappa statica

    radarelab::PolarScan<unsigned char> bb_first_level;         ///< mappa di elevazioni da beam blocking (input)
    radarelab::PolarScan<unsigned char> beam_blocking;          ///< mappa di beam blocking (input)

    radarelab::algo::Anaprop<double> anaprop;           ///< Oggetto per correzione ANAPRO

    //variabili legate a propagazione e beam blocking, da prog_bb
    radarelab::PolarScan<float> dem;                ///< dem in coordinate azimut range

    // metrici qualita' come sopra
    radarelab::Volume<unsigned char> qual;             ///< qualita volume polare
    // top, come sopra
    radarelab::PolarScan<unsigned char> top;            ///< Echo top a ???? dBZ [hm]

    /* variabili tolte perchè non presenti nel codice cum_bac... controllare che non richiamino qualcosa nelle funzioni
       struct tm *time_dbp;
       T_time, T_data, T_ora..*/

/**
 * Constructor
 * @param volume - reflectivity data
 * @param cfg - Configuration object
 * @param site - Site information
 * @param medium - force processing as medium range 
 * @param max_bin - maximum beam_size length 
 */
    CUM_BAC(radarelab::Volume<double>& volume, const Config& cfg, const Site& site, bool medium=false, unsigned max_bin=512);
    ~CUM_BAC();

    /**
     * Call this just after creating the CUM_BAC object, to signal that VPR
     * should also be computed
     */
    void want_vpr();

    /**
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
    void declutter_anaprop();

/// Esegue tutta la catena vpr (e classificazione) se richiesta
    void vpr_class();

    /**
     *  @brief funzione che caratterizza i volumi polari tramite la qualita'
     *  @details utilizza i parametri di ouput dell'elaborazione e la PIA
     *  calcolata qui per calcolare un valore finale di qualita' del dato.
     *  Inoltre calcola il top dell'eco in base a soglia su ogni pixel.
     */
    void caratterizzo_volume();

    /**
     *  @brief   funzione scrittura matrici statistica
     *  @details scrive le statistiche di beam blocking, anaprop, cambio di elevazione in un unsigined char DIM1_ST*DIM1_ST
     */
    void ScrivoStatistica(const radarelab::algo::anaprop::GridStats&);

    /**
     *  @brief funzione che legge la mappa statica e la mappa di elevazioni da beam blocking e le condensa in un unica mappa
     *  @details crea una mappa di elevazioni scelte in ogni pixel del volume per ottenere successivamante la stima di R
     *  @return non ritorna nulla
     */
    void leggo_first_level();

    /**
     *  funzioni di conversione cartesiana:associa a pixel matrice alta ris
     *  azimut e range, crea alta risoluzione e crea bassa risoluzione  
     *
     *  @brief funzione che calcola range e azimut su un quadrante centrato in 0,0
     *  @details
     */
    void creo_matrice_conv();

/// Nei punti convettivi ricalcola la Z come MP classica, partendo dalla stima di R (convettiva) 
    void conversione_convettiva();
    void generate_maps(CartProducts& products);
// added function to calculate beamblocking correction

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
};

/**
 * Struttara per il calcolo del VPR
 */
struct CalcoloVPR
{
    log4c_category_t* logging_category;

    CUM_BAC& cum_bac;   ///< oggeto CUM_BAC di riferimento
    radarelab::algo::InstantaneousVPR inst_vpr;
    long int gap;   ///< distanza temporale dall'ultimo file vpr [numero acquisizioni intercorse dall'ultimo vpr ?)
    float t_ground; ///< 2m temperature
    //matrici che dicono se pixel convettivo secondo VIZ, STEINER, riassuntiva mette +50
    radarelab::PolarScan<unsigned char> conv;   /// Informa se il pixel è convettivo
    radarelab::algo::VPR vpr;     ///< vpr 
    int hvprmax;            ///< quota picco vpr
    //elab classificazione: lista punti convettivi, iaz e ira, le dimensioni sono le massime possibili, in realtà i punti sono molti meno
    //int lista_conv[NUM_AZ_X_PPI*MAX_BIN][2];
   
    float resol[2];             ///< array di parametri, fisso , RES_HOR_CIL E RES_VERT_CIL
    int heating;            ///< variabile di riscaldamento
    int livmin = 0;             ///< quota livello minimo calcolato
    double htbb;                ///< altezza top brightband
    double hbbb;                ///< altezza bottom brightband
    radarelab::PolarScan<unsigned char> corr_polar; ///< correzione vpr in byte 0-128 negativa 128-256 positiva, in coord az-ra
    radarelab::PolarScan<unsigned char> neve;       ///< matrice az-range che memorizza punti di neve
    int ier_max;                ///< flag d'errore su calcolo quota max 
    int ier_stampa_vpr;             ///< flag d'errore su stampa profilo

    unsigned MyMAX_BIN;             ///< LUNGHEZZA MASSIMA 

/**
 * Constructor
 * @param [in] cum_bac - CUM_BAC object 
 */
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
    int profile_heating(bool has_inst_vpr);

    /**
     *  trova il massimo del profilo
     *  @brief funzione che trova la quota del massimo valore del profilo
     *  @details trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca come massimo di almeno 5 dBZ più alto in quota
     *  @param [out] hmax  - altezza massimo
     *  @return 0 se ok 1 se fallisce
     */
    int trovo_hvprmax(int *hmax);

    /**
     *  @brief funzione che combina il profilo verticale corrente con quello precedente tramite il metodo di Germann
     *  @details  oltre a lanciare il calcolo del profilo istantaneo provvede alla combinazione del profilo calcolato con il precedente calcolato entro  un limite massimo di distanza temporale pari a 10 quarti d'ora.  restituisce un codice integer pari a 0 se ok 1 se fallisce 
     *  @return 0 se combinazione ok 1 se fallisce
     */
    int combina_profili(const radarelab::algo::InstantaneousVPR& inst_vpr);

    /**
     *  @brief funzione  che classifica la precipitazione se stratiforme o convettiva
     *  @details esegue anche il ricampionamento cilindrico al suo interno
     */
    void classifica_rain();

    /**
     *  correzione vpr
     *  @brief funzione che corregge per il profilo verticale
     *  @details ciclando su tutti i bins della cartesiana polare scelta per la stima della pioggia,
     * - trovo la quota del centro del pixel
     * - correggo se: hbin>hliq , valore maggiore di soglia correzione (0 dBZ)
     * - se sefinita CLASS e pixel convettivo non correggo 
     * @return 0 se ok 1 se fallisce
     */
    int corr_vpr();

    /**
     *  fa il merge dei metodi
     *  @brief funzione  che interseca i punti convettivi delle due classificazioni Viz e Steiner e sottrae quelli con  picco stratiforme
     *  @param [in] steiner - oggetto che contiene classificazione secondo Steiner
     *  @param [in] viz - oggetto che contiene classificazione viz
     */
    void merge_metodi(const radarelab::algo::CalcoloSteiner& steiner, const radarelab::algo::CalcoloVIZ& viz);

    /// stampa profilo combinato
    /// @return codice errore
    int stampa_vpr();

/**
 * Metodo che lancia tutte le operazioni per il calcolo e la correzione del vpr
 */
    void esegui_tutto();
};

/**
 * Struttura per gestire output cartesiano di un singolo campo
 */ 
struct SingleCart
{
    const unsigned max_bin;     ///< dimensione matrice

    /// vol_pol interpolated in a cartesian map
    radarelab::Image<unsigned char> cart;       
/** 
 * Constructor
 * @param [in] max_bin - dimensione matrice
 */
    SingleCart(unsigned max_bin);

    /**
     *  conversione da polare a cartesiano alta risoluzione
     *
     *  @brief funzione che crea l'output cartesiano dal polare
     *  @details cicla sui quadranti e su i e j, usando il range e l'azimut ottenuti tramite la funzione creo_matrice_conv()
     *  NOTA Questa funzione mappa il campo in input in modo brutale su un unsigned char. DA SISTEMARE.
     *  @param [in] volume - volume passato per grafica
     *  @param [in] el_index - indice elevazione richiesta da graficare
     */
    void creo_cart(const radarelab::Volume <double> & volume, unsigned int el_index);

/**
 * Produce in output le immagini PNG dei campi in $DIR_DEBUG
 * @param [in] assets
 * @param [in] tagname  - Nome file di output 
 * @param [in] format - formato grafico richiesto
 */
    void write_out(Assets& assets, const std::string tagname, const std::string format="PNG");
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

}

#endif
