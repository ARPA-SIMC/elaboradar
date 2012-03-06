/*
  COMSTART cum_bac_sviluppo
  ----------------------------------------------------------------------------------------------
  | CUMULO SU BACINI
  |
  | Questo programma legge i file volumi polari DATAMAT 
  | scopo del programma e' eliminare il clutter tramite una maschera che può essere stagionale 
  | implementare un algoritmo di rimozione della propagazione anomala
  | un algoritmo di correzione e riduzione del beam_blocking (per il corto)
  | la correzione del VPR (se attivata)
  | calcolare la qualità del dato
  | come presentato al RADME 96 -ERAD 2004- ERAD 2006 -ERAD 2008 etc.
  | Il programma provvede inoltre a convertire il PPI, ripulito, a minor
  | elevazione in formato cartesiano, testare la classificazione stratiforme convettiva..
  | estrarre se presente la precipitazione sui bacini regionali,
  | e a  mediare la Z 
  | Sui file polari prima di essere utilizzati viene verificata la consistenza
  | con i seguenti criteri :
  |  - deve essere presente la rimozione del clutter tramite filtro doppler
  |  - devono essere presenti almeno NUM_MIN_BEAM per ognuna delle prime
  |    4 elevazioni
  |  - il file deve essere piu' recente dell'ultimo file processato
  |    ( questa condizione e' stata aggiunta per facilitare l'utilizzo 
  |       operativo a S.P.C.)
  |-----------------------------------------------------------------------------------------------
  | Nella nuova versione (2006) MODIFICHE PER L'IMPULSO CORTO: 
  |  -la maschera oltre al clutter climatologico puo' tenere 
  |   in cons. il beam blocking, e la traiettoria del lim inf. del fascio. 
  | - Inoltre se compilato con QUALITY il programma produce mappe di indici di qualita'
  | - il controllo anaprop è diverso a seconda che ci si trovi alla prima o successive elevazioni
  |----------------------------------------------------------------------------------------------- 
  |  AUTORI: Alberoni PierPaolo  S.M.R.  (modifiche A. Gioia, A. Fornasiero)
  |-----------------------------------------------------------------------------------------------
  |cc -o cum_bac_SP20_BSA cum_bac_SP20_BB.c  -DANAPROP -DTIME  -DZ_AVERAGE -DZ_LOWRIS  
  |       -DSHORT -DSTATIC  -DBEAMBLOCKING -DQUALITY  -I$INCLUDEDIR -lSP20_utility -L$HOME_BB/lib -lm
  |-------------------------------------------------------------------------------------------------
  | MODIFICHE
  |
  |
  | 31-07-1996 - Aggiunta la costruzione e la scrittura di una matrice di 
  |              riflettivita' alla risoluzione di 1Kmx1Km in output
  |              l'opzione e' attivabile tramite la dichiarazione di una variabile 
  |		in fase di 
  |
  | 06-08-1996 - Scrittura file di log
  |
  | 15-10-1996 - trovato errore nell'estrazione precipitazione sui bacini
  |		la matrice era stata trasposta
  |
  | 02-12-1996 - Aggiunta la routine write_xdr che scrive il vettore delle precipitazioni
  |   		medie sui bacini in formato xdr.
  |
  | 07-02-1997 - Ristretto il caso su cui viene fatto il test di anaprop:
  |                     solo se (bin_low > fondo_scala && bin_high >= fondo_scala)
  |              Trovato errore in statistica anaprop scriveva solo 0 o 1 ora calcola 
  |              la percentuale e la scrive come intero (0-100)
  |
  | 10-08-1998 - Corretto errore quando veniva fatto il test su anaprop.
  |              Versione precedente:
  |			if (bin_low > fondo_scala && bin_high >= fondo_scala)
  |			{
  |				----
  |			}
  |			else
  |			{
  |				----
  |			}
  |              Versione corretta:
  |			if (bin_low > fondo_scala && bin_high >= fondo_scala)
  |			{
  |				----
  |			}
  |			else if (bin_low < fondo_scala)
  |			{
  |				----
  |			}
  |			else if (bin_low = fondo_scala)
  |			{
  |				----
  |			}
  |		in questo modo se il raggio inferiore e` mancante riporta giu` il raggio
  |		superiore, mentre se e` uguale a fondo_scala riporto fondo_scala su tutti
  |		i livelli inferiori
  |
  | 01-08-2003 - Sono state definite due nuove funzioni,NewReadHeader e NewReadBeam,
  |              le quali permettono di leggere i dati relativi all'header e al raggio
  |              nel nuovo formato, reso compatibile con il vecchio formato (DATAMAT MDB_2.0).
  |
  | 17-09-2003 - Sono state introdotte le seguenti variabili di compilazione ANAPROP BACINI SHORT MEDIUM
  |              questo permette di attivare la rimozione del  clutter e della propagazione anomala,
  |              del calcolo della cumulata a livello di bacino e permette di gestire i file a 
  |              impulso corto e impulso medio.
  |              Vedi sotto per una descrizione completa delle variabili di compilazione
  !
  | 05-02-2004 - Introdotta l'opzione di compilazione BEAMBLOCKING per tenere conto delle mappe di BB
  |              A. Fornasiero
  | 09-03-2004 - Introdotta l'opzione di compilazione STATIC per scegliere se leggere la mappa statica
  |              A. Fornasiero
  | 08-2007      Introdotta stima e correzione VPR
  |              A. Fornasiero
  |-------------------------------------------------------------------------------------------------------------
  | 
  |  NB $HOME_BB_VPR è la directory dove sta il pacchetto---$WORKDIR quella di lavoro (per la versione studio e'
  |                  /scratch/fornasiero/DATI_OUT/CUM_BAC/cum_bac_tmp_BB_VPR/[SPC,GAT]_short, per quella ope è 
  |                  $HOME_BB_VPR/cum_bac_tmp_BB_VPR/[SPC,GAT]_short)
  |  Files input: 
  |            $HOME_BB_VPR/DB/DBP2_[aaaammgghhmm]_BOLOGNA(GATTATICO): file vol. polare
  |            $HOME_BB_VPR/MAPPE_STATICHE/FIRST_LEVEL_[....]: file mappa statica, il nome cambia a seconda della stagione
  |            $WORKDIR/anap_last.$$: indica l'ultima data analizzata (LAST_FILE)
  |            $HOME_BB_VPR/lista_date
  |            
  |          *** da elaborazione beam blocking ******
  |            $HOME_BB_VPR/PP+BLOC/output/[aaaammgghhmm]mat_el.bin:file mappa dinamica da beam blocking
  |            $HOME_BB_VPR/PP+BLOC/output/[aaaammgghhmm]mat_bloc.bin: file mappa beam blocking residuo
  |            $HOME_BB_VPR/PP+BLOC/dem_SanPi(Gatta).txt: file dem dell'area radar in coord. az-range (su 125 km raggio)
  |            $HOME_BB_VPR/PP+BLOC/[aaaammgghhmm]hray.txt: file traiettorie del centro fascio e diverse elev
  |            $HOME_BB_VPR/VPR/vpr_par.h: parametri calcolo VPR
  |  Files output: 
  |            $WORKDIR/log_file: file di log.. importante perchè memorizza anche le variabili di compilazione
  |            $WORKDIR/stat_anap(ANAP_STAT_FILE) $WORKDIR/stat_elev(ELEV_STAT_FILE)  $WORKDIR/stat_bloc(BLOC_STAT_FILE):
  |            v. variabili ambiente: files statistiche
  |            $HOME_BB_VPR/INPUT_ZLR_BB/[aaaammgghhmm].ZLR: file matrice cart. unsign.char valori riflettività [0,255] sta per [-20,80]
  |            $HOME_BB_VPR/INPUT_ZLR_BB/[aaaammgghhmm].qual_ZLR:file matrice cart. unsign.char valori qualita' ZLR
  |            $HOME_BB_VPR/INPUT_ZLR_BB/[aaaammgghhmm].qual_RLR:file matrice cart. unsign.char neve ( 0=preci e 1=neve)
   |            $WORKDIR/QUALITY/[aaammgghhmm].elev: file matrice unsigned char elevazioni usate 
  |            $WORKDIR/QUALITY/[aaammgghhmm].quota: file matrice unsigned short quota fascio sul suolo 
  |            $WORKDIR/QUALITY/[aaammgghhmm].bloc:  file matrice unsigned char % beam blocking
  |            $WORKDIR/QUALITY/[aaammgghhmm].corrpt: file matrice unsigned char classific. anaprop [0 ok, 1 anap, 2 no dato, 3 no test]
  |            $HOME_BB_VPR/VPR/last_vpr_[GAT,SPC]: file data ultimo vpr
  |            $HOME_BB_VPR/VPR/vpr_heat_[GAT,SPC]: file riscaldamento vpr
  |            $HOME_BB_VPR/VPR/vpr_[GAT,SPC]: file vpr 
  |            $HOME_BB_VPR/VPR/vpr_[GAT,SPC]_int: file vpr interpolato (per plot)
  |------------------------------------------------------------------------------
  | Vengono utilizzate le seguenti variabili di ambiente :
  |
  | LISTA_FILE          - nome del file contenete i file dati da utilizzare
  | LAST_FILE           - nome del file contenente l'ultima data processata
  | TIPO_FILE           - numero [0,3] indica il tipo di dato in ingresso (0=corto dec. 1= corto full volume 2=corto grandine  3= medio)
  | ANAP_STAT_FILE      - nome del file contenente la statistica sull'anaprop
  | BLOC_STAT_FILE      - nome del file contenente la statistica sulla perc. beam blocking
  | ELEV_STAT_FILE      - nome del file contenente la statistica sul cambio elevazione rispetto mappa statica
  | FIRST_LEVEL_FILE    - nome del file contenente la mappa statica del primo livello 
  |                       sopra il clutter
  | DIR_OUT_PP+BLOC     - directory file propagaz. e beam blocking
  | FILE_DEM(_SPC o _GAT)- file contenente il dem
  | MATRICE_BACINI      - nome del file contenente la matrice per la 
  |                       suddivisione dei bacini
  | BACINI_AUS_FILE     - nome file di lavoro per i bacini
  | BACINI_DIR          - nome directory di output per i file contenenti 
  |                       i dati sui bacini
  | BACINI_HISTORY_FILE - nome del file contenente la storia delle 
  |                       precipitazioni sui bacini
  | OUTPUT_Z_LOWRIS_DIR - directory di output per i file di riflettivita'
  |    			 su un grigliato di 1 Km x 1 Km
  |
  | LOG_FILE            - file di log
  |
  | BACINI_BOLOGNA      - File di uscita della precipitazione cumulata sui bacini
  |			 in formato xdr pronto per essere inviato a Bologna
  |
  | DIR_QUALITY         -directory di output degli indici di qualita'
  | 
  |--------------------------------------------------------------------------------
  | VARIABILI DI COMPILAZIONE
  | 
  | In fase di compilazione e' possibile richiedere che vengano eseguite
  | soltanto alcune operazioni tramite la definizioni di alcune variabili :
  |
  | ANAPROP    -  Attiva cancellazione clutter e propagazione anomala
  |
  | WRITE_DBP  -  permette di scrivere in output il volume polare 
  |               dopo la rimozione della propagazione anomala 
  |               attivando questa opzione l'elaborazione
  |               del volume polare si ferma dopo la rimozione
  |
  | TIME       -  vengono stampate sullo standard output informazioni 
  |               riguardanti il tempo impiegato da ogni singola operazione
  |               nella porcessazione ... appunto, porcessazione..
  |
  | WRITE_DBP_REORDER -  permette di scrivere in output il volume polare
  |                      dopo la fase di rilettura e di riordino ad azimut
  |                      fissi ed elevazioni fisse
  |
  | DECLUTTER -  permette di scrivere in output il volume polare
  |              dopo la fase di rimozione del clutter senza attivare la 
  |              rimozione della propagazione anomala
  |
  | Z_LOWRIS  -  abilita la costruzione e la scrittura di una matrice di
  |              riflettivita' a bassa risoluzione (1km)
  |
  | SIDME     -  abilitava la scrittura del file per procedura sidme, deve essere definita
  |              la variabile d'ambiente NOME_SIDME
  |
  | BACINI    -  abilita il calcolo e las crittura (ole') in output della cumulata a livello di bacini
  |
  | SHORT     -  configura il programma all'esecuzione per dati provenienti da impulso corto (125 km range; 250 m ris)
  |
  | MEDIUM    -  configura il programma all'esecuzione per dati provenienti da impulso medio (250 km range; 1000 m ris)
  | BEAMBLOCKING - abilita la lettura delle matrici di Beam Blocking e la correzione della misura
  |
  | QUALITY      - abilita la generazione e scrittura variabili con un un contenuto di qualità:
  |               bb residuo, hrelativa, dati 'corrotti'(anaprop,manca dato,aseenza controllo anaprop)
  |
  | STATIC       - abilita la lettura della mappa statica
  | 
  |
  | BLOCSTAT    - con quest'opzione la correzione si limita a mappa statica + rimozione ANAPROP, come nella prima versione.
  |               Viene attivata  quando si desidera  memorizzare il valore del beam blocking senza che venga corretto. 
  |               Si tiene memoria del Beam Blocking inferiore al 50%, per valori superiori si pone BB=51%
  |             
  | VPR         - attiva calcolo e correzione VPR: agosto 2007 deve essere abilitato anche QUALITY
  |-----------------------------------------------------------------------------
  |
  | Queste varialibi (sig.ra Fletcher!) di compilazione vengono momentaneamente a perdere di validità 17-09-2003
  | 
  | Deve essere inoltre definita una delle seguenti variabili che individua
  | l'ambiente di esecuzione del programma e utilizza i file di include 
  | Datamat corretti:
  |
  | BOLOGNA - per compilare il programma a Bologna
  |
  | SPC     - per compilare il programma a San Pietro Capofiume
  |
  |-----------------------------------------------------------------------------
  |
  | In fase di link deve essere specificata la libreria matematica
  |
  |-----------------------------------------------------------------------------
  |
  comend */
/*----------------------------------------------------------------------------*/
/*	INCLUDE file						       	      */
/*----------------------------------------------------------------------------*/
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
#include <func_SP20read.h>
#include <radar_parameter.h>
#include <errno.h> /* 15-2-2006 anna messa per evitare errore /lib/libc.so.6: could not read symbols: Bad value */
#include <MP_par.h>/**/
#include <vpr_par.h>
#include <geo_par.h> /* lat e lon radar */

#include <nrutil.h>
#include <nr.h>


#ifdef CLASS
#include <par_class.h>
#endif
#define AMPLITUDE 0.9 /* esternalizzo?*/
#define DTOR  M_PI/180. /* esternalizzo?*/
#define CONV_RAD 360./4096.*DTOR
#define DIM1_ST 16 /**/
#define DIM2_ST 13/*CAMBIATA DIMENSIONE  A 13 PER NUOVA ESTENSIONE RAGGIO RADAR*/
#define MAX_BIN 512
#define MAX_DIM 512
#define LIMITE_ANAP 240/* LIMITE per cambiare controllo anaprop*/
#define FATT_MOLT_EL ((double) 360./(double)4096.)
#define FATT_MOLT_AZ ((double) 360./(double)4096.)
#define NSCAN 6        /* numero elevazioni */
#define OVERBLOCKING 51 /* minimo BB non accettato*/
#define THRES_ATT 0 /* minimo valore di Z in dBZ per calcolare att rate */

// Ridefinisco NEL
#ifdef NEL
#undef NEL
#endif

#ifdef SHORT
#define NEL 10


#define RANGE_KM                 125

#define CART_DIM_ZLR             256
#define ZLR_N_ELEMENTARY_PIXEL   4
#define ZLR_OFFSET               0
#define NMIN 5
#define MAX_TIME_DIFF 3
int elev_array[NEL];

int  elev_array_spc[NEL]={6,15,26,36,46,57,80,108,148,159};//l'ultimo è fittizio:ANNA 30-03-2011
int  elev_array_gat[NEL]={6,16,27,37,45,55,65,76,85,95};//105,126,149,174,201 è da completare NEL=15:ANNA 30-03-2011
#endif

#ifdef MEDIUM
#define NEL 5
#define RANGE_KM                 250

#define CART_DIM_ZLR             512
#define ZLR_N_ELEMENTARY_PIXEL   1
#define ZLR_OFFSET               CART_DIM_ZLR/2
int elev_array[NEL];
int elev_array_spc[NEL]={6,16,26,36,47};//ANNA 30-03-2011
int elev_array_gat[NEL]={6,16,27,36,47};//ANNA 30-03-2011
#define NMIN 2
#define MAX_TIME_DIFF 1
#endif



/*   #include <datamat_file.h> */

/*----------------------------------------------------------------------------*/
/*	DICHIARATIVE GLOBALI						      */
/*----------------------------------------------------------------------------*/
T_MDB_ap_beam_header  old_beam_header;  /* queste due dichiarazioni le metto qui perchè sono sfaticato*/
T_MDB_data_header   old_data_header;    /* P.P. Sei un mito! */

/*extern char *sys_errlist[];  ANNA 17-2-2006 sostituito con strerror che sta in string.h */
extern int errno;
char errori[256];
                                                               

#define MAX_DIF_OR 30            /* differenzio limiti controllo anap      */
#define MIN_VALUE_OR -10         /* a seconda che sia alla prima o success.*/
#define MAX_DIF_NEXT_OR 15       /*/ elevazione                             */
#define MIN_VALUE_NEXT_OR 0 
#define MAX_DIF_LIMIT 55     /*  soglia diff tra valori bin per anaprop in caso il bin si trovi oltre LIMITE_ANAP (60 KM) */
#define MIN_VALUE_LIMIT -20  /*  soglia valore minimo bin per anaprop oltre LIMITE_ANAP (60 KM) */
int MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT;/* differenza massima tra le due elevazioni successive perchè non sia clutter e valore minimo a quella superiore pe il primo e per i successivi (NEXT) bins*/
#define MISSING 0 /*valore mancante*/
#define NUM_MIN_BEAM 200
#define PIOVE     1
#define NON_PIOVE 0
#define INIZIO_EVENTO PIOVE
#define FINE_EVENTO   NON_PIOVE
#define DIFF_TIME 1*60*60         /*--- tempo utile per chiudere evento ---*/
#define RAIN_MIN  0.5             /*--- soglia minima presenza rain ----*/
#define STEP_STAT_ANAP_RANGE  40 /*dim range griglia per stat anap*/
#define STEP_STAT_ANAP_AZ     25 /*dim azim griglia per stat anap*/
#define N_MIN_BIN           500 /*--- numero minimo di celle presenti in un 
				  settore per la statistica            ---*/
#define SHORT_DEC         0
#define SHORT_FULL_VOLUME 1
#define SHORT_HAIL        2
#define MEDIUM_PULSE      3 

/*-----------------------------------------------------------
  Variabili globali
  ------------------------------------------------------------*/
int argc;
char **argv;

struct VOL_POL vol_pol[NEL][NUM_AZ_X_PPI]; 

char *nome_fl,*nome_dem;
unsigned char first_level[NUM_AZ_X_PPI][MAX_BIN];
unsigned char first_level_static[NUM_AZ_X_PPI][MAX_BIN];
#ifdef BEAMBLOCKING
unsigned char bb_first_level[NUM_AZ_X_PPI][MAX_BIN];  /* mappa di elevazioni da beam blocking (input)*/
unsigned char beam_blocking [NUM_AZ_X_PPI][MAX_BIN];   /* mappa di beam blocking (input)*/
#endif

#ifdef QUALITY  
float hray[MAX_BIN][NEL];  /*quota centro fascio in funzione della distanza e elevazione*/
float hray_inf[MAX_BIN][NEL]; /*quota limite inferiore fascio in funzione della distanza e elevazione*/
float dem[NUM_AZ_X_PPI][MAX_BIN]; /*dem in coordinate azimut range*/
unsigned short quota_rel[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio in coordinate azimut range*/
unsigned short quota[NUM_AZ_X_PPI][MAX_BIN]; /*quota fascio in coordinate azimut range*/
unsigned short quota_cart[MAX_BIN*2][MAX_BIN*2];/*quota fascio in coordinate cart 1024*1024, risoluzione minima*/
unsigned char quota_1x1[CART_DIM_ZLR][CART_DIM_ZLR];/* quota in formato 256*256 in centinaia di metri, risoluzione ZLR */

unsigned char dato_corrotto[NUM_AZ_X_PPI][MAX_BIN]; /*uscita controllo anaprop in coordinate azimut range */
unsigned char elev_fin[NUM_AZ_X_PPI][MAX_BIN]; /* elevazione finale in coordinate azimut range  */ 
unsigned char qual[NEL][NUM_AZ_X_PPI][MAX_BIN]; /* qualita volume polare */
unsigned char qual_Z_cart[MAX_BIN*2][MAX_BIN*2]; /* qualita della Z in formato 1024*1024, risoluzione minima */
unsigned char qual_Z_1x1[CART_DIM_ZLR][CART_DIM_ZLR];/* qualita della Z in formato 256*256, risoluzione ZLR */

float zeroterm;
float dtrs;
unsigned char att_cart[NUM_AZ_X_PPI][MAX_BIN]; /* matrice azimut-range di attenuazione */
#ifdef VPR
unsigned char  corr_polar[NUM_AZ_X_PPI][MAX_BIN];/*correzione vpr in byte 0-128 negativa 128-256 positiva*/
unsigned char  corr_cart[MAX_BIN*2][MAX_BIN*2];
unsigned char corr_1x1[CART_DIM_ZLR][CART_DIM_ZLR];
unsigned char flag_vpr[NEL][NUM_AZ_X_PPI][MAX_BIN];/* punti del volume polare ok per calcolo VPR*/
unsigned char neve_cart[MAX_BIN*2][MAX_BIN*2];/* qualita della RAIN RATE in formato 1024*1024, risoluzione minima */
unsigned char neve_1x1[CART_DIM_ZLR][CART_DIM_ZLR];/* qualita del Rain Rate in formato 256*256, risoluzione ZLR */
float vpr[NMAXLAYER];/* vpr */
long int area_vpr[NMAXLAYER]; /*area degli strati*/
int ier_vpr, ier_comb,ier_stampa_vpr;/* flag d'errore su calcolo vpr istantaneo, combinazione vpr */
int hvprmax; /* quota picco vpr */
/*unsigned char ***vol_vpr;  non usato  */
unsigned char neve[NUM_AZ_X_PPI][MAX_BIN];/* matrice az-range che memorizza punti di neve*/
int heating,livmin=-1; /* variabile di riscaldamento e quota livello minimo calcolato*/
FILE *log_vpr,*test_vpr;
float t_ground;
int data[3],T_data[3],T_ora[2],ora[2],nstaz,nvar;
struct tm *T_tempo;
float stdev;
#endif
#endif
float azimut[MAX_BIN][MAX_BIN];
float range[MAX_BIN][MAX_BIN];
unsigned char cart[MAX_BIN*2][MAX_BIN*2];
double  cartm[MAX_BIN*2][MAX_BIN*2];  /* Z media dei bins adiacenti al punto */
int stat_anap[DIM1_ST][DIM2_ST]; /* statistica anaprop  */
int stat_anap_tot[DIM1_ST][DIM2_ST]; /* contatore punti dentro ogni box per statistica */
int long stat_bloc[DIM1_ST][DIM2_ST];   /* statistica beam blocking  */
int stat_elev[DIM1_ST][DIM2_ST]; /* statistica cambio elevazione rispetto mappa statica  */     
int n_elev; /* cambiato n elevazioni */ 

struct tm *time_dbp;

float size_cell[]={62.5,125.,250.,500.,1000.,2000.};
int nbeam_elev[NEL];

#ifdef Z_LOWRIS
unsigned char z_out[CART_DIM_ZLR][CART_DIM_ZLR];
#endif
unsigned char MP_coeff[2]; /* a/10 e b*10 per scrivere come 2 byte */
float aMP, bMP;   /*  coeff a e b relazione Z-R  */
long int gap; /* distanza temporale dall'ultimo file vpr */
struct tm *tempo; /* aggiunti nel main per leggere stagione dal nome file e ricavere MP coeff */ 
time_t Time;
int month;
float calibr=0.0; /* fattore di calibrazione radar in DB */

/*----------------------------------------------------------------------------*/
/*	FUNCTION PROTOTYPE						      */
/*----------------------------------------------------------------------------*/
float RtoDBZ();
float BYTEtoDB();
float  DBZtoR();
int elabora_dato();
int write_dbp();
void leggo_first_level();
void creo_matrice_conv();
void creo_cart();
float BYTE_to_mp_func();  /* trasforma dato byte in pioggia */

unsigned char DBtoBYTE(); /* trasforma dato di riflettivita' dBZ in valore byte 0-255 */
float BYTEtoZ();          /* trasforma dato byte in Z (non dBZ) */
float Z_to_mp_func();     /* trasforma Z (non dBZ) in pioggia */
time_t 	NormalizzoData();
int test_file();
void bacini();
void prendo_tempo();
void AggiornoHistoryFile();
void ScrivoBacini();
void ScrivoStatistica();
#ifdef Z_LOWRIS
void creo_cart_z_lowris();
void scrivo_out_file_bin();
#endif
char *PrendiOra();
void ScrivoLog();
void   write_xdr();
FILE *controllo_apertura();

#ifdef QUALITY
void leggo_dem();      
void leggo_hray();   
double attenuation(unsigned char DBZbyte, double  PIA); 
void caratterizzo_volume();
float func_q_Z();
float func_q_R();
float quota_f(float elevaz, int k);
int trovo0term();

#ifdef VPR
int func_vpr();
float comp_levels();
int corr_vpr();
long int profile_gap();
int profile_heating();
int stampa_vpr();
int combina_profili(char *sito, char *sito_ad);
int get_t_ground(float *t_gr);
int analyse_VPR();
int interpola_VPR();
int testfit(float a[], float chisq, float chisqin);
int trovo_hvprmax(int *hmax);
/*int n_getfev1_(); */
int n_close_(); 
void lineargauss(float x, float a[], float *y, float dyda[],int na);

time_t Time,T_Time;
int ier_t ;
float chisqfin=100; //???puo' essere def in anal
#endif

#ifdef CLASS
unsigned char cappi[NUM_AZ_X_PPI][MAX_BIN],class_VIZ[NUM_AZ_X_PPI][MAX_BIN],*conv_VIZ[NUM_AZ_X_PPI],*conv_STEINER[NUM_AZ_X_PPI],*conv[NUM_AZ_X_PPI];
unsigned char conv_cart[MAX_BIN*2][MAX_BIN*2],cappi_cart[MAX_BIN*2][MAX_BIN*2], conv_1x1[CART_DIM_ZLR][CART_DIM_ZLR],cappi_1x1[CART_DIM_ZLR][CART_DIM_ZLR];
int lista_conv[NUM_AZ_X_PPI*MAX_BIN][2];
int lista_bckg[NUM_AZ_X_PPI*MAX_BIN][2];
double *Z_bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in mm^6/m^3
float *bckgr; // array contenente i valori della Z di background per ogni pixel precipitante in dB
float azimut[MAX_BIN][MAX_BIN];
float rangexy[MAX_BIN][MAX_BIN]; 

float **cil[NUM_AZ_X_PPI];
float resol[2];
int x_size,z_size;
int ncv=0,ncs=0,np=0;
float *convective_radius;
float htbb, hbbb;
//unsigned char CONV_VAL=(unsigned char)(200);

float THR_ZABB=1000.; //Zabb soglia

//unsigned char allert[MAX_BIN][MAX_BIN];

void classifica_rain();
void classifico_VIZ();
void calcolo_background();
void ingrasso_nuclei(float cr,int j,int k);
void classifico_STEINER();
void classifico_ZLR();
void merge_metodi();
#endif

#endif
int setstat();
int setwork();

/* ================================ */
int main (int argc, char **argv)
/* ================================ */
{
  char *nome_file , *nome_file_volume;
  int l,ier,ier_test, ier_main=0;
  int tipo_dati_richiesti = INDEX_Z;
   FILE  *output;
  char *sito;
  /*---------------------------------------
    definisco le variabili per estrarre 
    pezzi della matrice cartesiana
    
    ---------------------------------------*/

  int imin,imax;
  int jmin,jmax;
  int dx,dy;
  printf("sono in programma elaborazione\n");
  if (argv[1]==NULL || argv[2]==NULL || argv[3]==NULL) {
    printf("insufficiente numero di argomenti\n");
    exit(1);
  }

  nome_file_volume=argv[1];
  /*setto ambiente lavoro*/

  setwork(argv[3]);
    nome_file=nome_file_volume;
  ScrivoLog(0,nome_file);
  ScrivoLog(1,nome_file);
  ScrivoLog(2,nome_file);

  ScrivoLog(18,argv[1]); /*nome_file*/
  ScrivoLog(18,argv[2]);   /* TIPO_FILE */
  ScrivoLog(18,argv[3]);  /* SITO */
  ScrivoLog(18,argv[4]); /* SITO_AD */
   
   sito=argv[3];
  if (!(strcmp(sito,"SPC")) ) {
    printf("sito  SPC elevazioni \n" );
    for(l=0; l<NEL; l++){
      elev_array[l]=elev_array_spc[l];
      printf(" %i " , elev_array[l]);
      
    }
  }

  if (!(strcmp(sito,"GAT")) ) {
    printf("sito  GAT elevazioni:\n" );
    for(l=0; l<NEL; l++){
      elev_array[l]=elev_array_gat[l];
      printf(" %i " , elev_array[l]);
    }
  }

 
#ifdef TIME
  prendo_tempo();
#endif      
  /* pulisco tutte le matrici    */
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
  memset(elev_fin,0,sizeof(elev_fin));        
  memset(dato_corrotto,2,sizeof(dato_corrotto));                 
  memset(att_cart,DBtoBYTE(0.),sizeof(att_cart));
  memset(qual,0,sizeof(qual));
  memset(quota_rel,0,sizeof(quota_rel));
  memset(quota,0,sizeof(quota));
  memset(quota_cart,0,sizeof(quota_cart));
  memset(quota_1x1,0,sizeof(quota_1x1));
  memset(dato_corrotto,0,sizeof(dato_corrotto));
  memset(elev_fin,0,sizeof(elev_fin));
  memset(qual,0,sizeof(qual));
  memset(qual_Z_cart,0,sizeof(qual_Z_cart));
  memset(qual_Z_1x1,0,sizeof(qual_Z_1x1));
  #ifdef VPR
  memset(corr_1x1,0,sizeof(corr_1x1));
  memset(neve_cart,0,sizeof(neve_cart));
  memset(neve_1x1,0,sizeof(neve_1x1));
  memset(neve,0,sizeof(neve));
  #endif
#endif
#ifdef BEAMBLOCKING
  memset(bb_first_level,0,sizeof(bb_first_level)); 
#endif 
     
  nome_fl=(char *)malloc(200*sizeof(char));
  nome_dem=(char *)malloc(200*sizeof(char));

  ScrivoLog(4,nome_file);
  printf("processo file dati: %s\n",nome_file);
  ier = read_dbp_SP20(nome_file,vol_pol,&old_data_header,
		      tipo_dati_richiesti,nbeam_elev);

  printf("fatta lettura\n");      
  ScrivoLog(5,nome_file);

       
#ifdef TIME
  prendo_tempo();
#endif
       
  ier_test=test_file(argv[2]);
  printf("ier -- test  %d  %d\n",ier,ier_test);     
  if(ier == OK  && ier_test)
    {
#ifdef WRITE_DBP_REORDER
	  
      /*------------------------------------------------------------------
	| eventuale scrittura del volume polare ad azimut e range fissi  |
	-----------------------------------------------------------------*/ 
      strcat(nome_file,"_reorder");	  
      ScrivoLog(6,nome_file);
      ier=write_dbp(nome_file);
#endif	  
      if( NormalizzoData(old_data_header.norm.maq.acq_date) != -1 )  
	{
	  /*------------------------------------------
	    | rimozione propagazione anomala e clutter |
	    ------------------------------------------*/
	  ScrivoLog(7,nome_file);	   
	  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
	  tempo = gmtime(&Time);
	  month=tempo->tm_mon+1;
	  setstat(argv[3],month, nome_dem, nome_fl);
	  printf ("nome dem %s nome fl %s\n",nome_dem, nome_fl);
	  ScrivoLog(1,nome_file);
	  /*definisco i coeff MP in base alla stagione*/
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
	      
	  /*fine definizione coeff MP in base alla stagione*/
#ifdef ANAPROP
	  printf ("inizio rimozione anaprop e beam blocking \n") ;
	  ier = elabora_dato();  
#endif
#ifdef QUALITY
	  printf ("calcolo Q3D \n") ;    
	  caratterizzo_volume();


#ifdef VPR
	  log_vpr=fopen(getenv("LOG_VPR"),"a+");
	  nstaz=4;
	  nvar=1;
	  data[2]= tempo->tm_year+1900;
	  data[1]= tempo->tm_mon+1;
	  data[0]= tempo->tm_mday;
	  ora[0]=  tempo->tm_hour;
	  ora[1]= 00;/* cerco i dati all' inizio dell'ora */ 
	  fprintf(log_vpr," data ora  %i %i %i %i %i  %i %i\n",data[0],data[1],data[2],ora[0],ora[1],nstaz,nvar);
	  t_ground=NODATAVPR;
	  ier_t=get_t_ground(&t_ground);
	  fclose(log_vpr);	
#endif

#ifdef CLASS
	  log_vpr=fopen(getenv("LOG_VPR"),"a+");
	  classifica_rain();
	  fclose(log_vpr);
#endif

#ifdef VPR
	  log_vpr=fopen(getenv("LOG_VPR"),"a+");
	  test_vpr=fopen(getenv("TEST_VPR"),"a+");
	  fprintf(log_vpr,"processo file dati: %s\n",nome_file); 
	  printf ("calcolo VPR \n") ;
	  hvprmax=INODATA;
	  ier_comb=combina_profili(argv[3],argv[4]);
	  printf ("exit status calcolo VPR istantaneo: (1--fallito 0--ok)  %i \n",ier_vpr) ; // debug
	  printf ("exit status combinaprofili: (1--fallito 0--ok) %i \n",ier_comb) ; // debug
	  printf ("heating %i \n", heating);
	  heating=profile_heating();	      	     
	  if (!ier_comb && heating >= WARM){

	    /*CORREZIONE */
	    ier=corr_vpr(argv[3]);
	     printf ("exit status correggo vpr: (1--fallito 0--ok) %i \n",ier) ; // debug
	    if ( ! ier ) 
	      ier_stampa_vpr=stampa_vpr();
	  }
	  
	  fclose(log_vpr);

#endif


#endif
	      
#ifdef TIME
	  prendo_tempo();
#endif	      
#ifdef WRITE_DBP
	  /*-------------------------------------------------
	    | eventuale scrittura del volume polare ripulito e |
	    -------------------------------------------------*/ 
#ifdef DECLUTTER
	  strcat(nome_file,"_decl");
#else
	  strcat(nome_file,"_anap");
#endif
	  /*exit (1);*/
	  ScrivoLog(8,nome_file);
	  ier=write_dbp(nome_file);
#endif
#ifndef WRITE_DBP
	      
	  /*--------------------------------------------------
	    | conversione di coordinate da polare a cartesiana |
	    --------------------------------------------------*/
	  ScrivoLog(9,nome_file);
	  creo_cart();
#ifdef TIME
	  prendo_tempo();
#endif

	      
#ifdef Z_LOWRIS

	  ScrivoLog(13,nome_file);
	  creo_cart_z_lowris();
	  printf(" dopo z_lowris\n");
	  ScrivoLog(14,nome_file);
	  /*scrivo_z_lowris();*/
	  scrivo_out_file_bin(".ZLR","dati output 1X1",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(z_out),z_out);
	  sprintf(nome_file,"%s/MP_coeff",getenv("OUTPUT_Z_LOWRIS_DIR"));

	  output=controllo_apertura(nome_file,"file coeff MP","w");
	  
	  fwrite(MP_coeff,sizeof(MP_coeff),1,output);
	  fclose(output);
	  printf(" dopo scrivo_z_lowris\n");	      
#endif
	      
	      
#ifdef TIME
	  prendo_tempo();
#endif	
	      
#endif
	      
#ifdef QUALITY      
	  scrivo_out_file_bin(".corrpt","file anap",getenv("DIR_QUALITY"),sizeof(dato_corrotto),dato_corrotto); 
	  scrivo_out_file_bin(".pia","file PIA",getenv("DIR_QUALITY"),sizeof(att_cart),att_cart);
	  scrivo_out_file_bin(".bloc","file bloc",getenv("DIR_QUALITY"),sizeof(beam_blocking),beam_blocking); 
	  //scrivo_out_file_bin(".quota","file quota",getenv("DIR_QUALITY"),sizeof(quota),quota);
	  scrivo_out_file_bin(".elev","file elevazioni",getenv("DIR_QUALITY"),sizeof(elev_fin),elev_fin);  
	  scrivo_out_file_bin(".quota_ZLR","file quota",getenv("DIR_QUALITY"),sizeof(quota_1x1),quota_1x1);// m/100 +128

#ifdef VPR  
	  scrivo_out_file_bin(".corr_ZLR","file quota",getenv("DIR_QUALITY"),sizeof(corr_1x1),corr_1x1);
	  scrivo_out_file_bin(".neve","punti di neve",getenv("DIR_QUALITY"),sizeof(neve),neve); 
	  scrivo_out_file_bin(".neve_ZLR","file presunta neve ",getenv("DIR_QUALITY"),sizeof(neve_1x1),neve_1x1); 
#endif	
#ifdef CLASS  
	  scrivo_out_file_bin(".conv_ZLR","punti convettivi",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(conv_1x1),conv_1x1); 
#endif	   
	  scrivo_out_file_bin(".qual_ZLR","file qualita' Z",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(qual_Z_1x1),qual_Z_1x1); 
#endif
	}
    }
  else ier_main=1;
  //}  COMMENTO FINE LOOP*/ 

  
  /*-----------------------------
    | fine loop sui volumi polari |
    -----------------------------*/
  //fclose(file2);    
  
  ScrivoLog(15,nome_file);
  
  return ier_main;

}

/*-------------------------------------------
  | FUNCTION : test_file		     |
  | 					     |
  -------------------------------------------
  | Verifica se il volume dati in input e'   |
  | consistente con le esigenze del programma|
  | la funzione ritorna :                    |
  |   0  in caso di errore                   |
  |   1  in caso di successo                 |
  -------------------------------------------*/
int test_file(char *tipofile )
{
  FILE *f_aus;
  time_t last_time;
  int k;
  // int n_elev,resolution;
  int resolution;
  int file_type;

  sscanf(tipofile,"%d",&file_type);
  printf("tipo file %1d\n",file_type);
  switch (file_type)
    {
    case SHORT_DEC:  
      if(!old_data_header.norm.maq.declutter_rsp ) 
	{
	  sprintf(errori,"File Senza Declutter Dinamico");
	  ScrivoLog(16,errori);
	  return 0;
	}
      resolution=2;
      n_elev=4;
      break;
    case SHORT_FULL_VOLUME:
      if(old_data_header.norm.maq.declutter_rsp ) 
	{
	  sprintf(errori,"File con Declutter Dinamico");
	  ScrivoLog(16,errori);
	  return 0;
	}
      resolution=2;
      n_elev=4;
      break;
    case SHORT_HAIL:
      resolution=2;
      n_elev=3;
      printf("CASO SHORT_HAIL\n");
      break;
    case MEDIUM_PULSE:
      resolution=4;
      n_elev=4;
      break;
    }                                                       /* end switch */

  if(old_data_header.norm.maq.resolution != resolution)
    {
      sprintf(errori,"File Risoluzione Sbagliata %1d", old_data_header.norm.maq.resolution);
      ScrivoLog(16,errori);
      return 0;
    }
 

  for(k=0; k<n_elev; k++) /* testo solo le prime 4 elevazioni */
    {

      printf(" numero beam presenti : %4d  -- elevazione%2d\n",
	     nbeam_elev[k],k);

      if(nbeam_elev[k] <  NUM_MIN_BEAM)
	{
	  sprintf(errori,"Trovati Pochi Beam Elevazione %2d - num.: %3d",k,nbeam_elev[k]);
	  ScrivoLog(16,errori);
	  return 0;
	}
    }                                                             /*end for*/
  sprintf(errori,"primi test passati");
  ScrivoLog(16,errori);
  /*-------------------------------------------------------------------
    | verifico la presenza del file contenente l'ultima data processata |
    -------------------------------------------------------------------*/
  if(access(getenv("LAST_FILE"),6) == 0)
    {
      /*--------------------------------------
	|  il file e' presente, leggo la data  |
	|  e la confronto con la data del      |
	|  volume dati in esame, se il dbp e'  |
	|  piu' giovane continuo altrimenti    |
	|  esco dalla function con status di   |
	|  errore (0).                         |
	-------------------------------------*/
      f_aus = fopen(getenv("LAST_FILE"),"r+");
      fread(&last_time,4,1,f_aus);
      if(old_data_header.norm.maq.acq_date <= last_time)
	{
	  fclose(f_aus);
	  sprintf(errori,"File Vecchio");
	  ScrivoLog(16,errori);
	  //return 0;
	}
      else{
      /*----------------------------
	|  aggiorno la data nel file |
	----------------------------*/
      rewind(f_aus);
      fwrite(&old_data_header.norm.maq.acq_date,4,1,f_aus);
      fclose(f_aus);
      }
    }
  else
    {
      /*--------------------------------------
	|  il file non e' presente, scrivo la  |
	|  data del volume in esame            |
	--------------------------------------*/
      f_aus = fopen(getenv("LAST_FILE"),"w");
      fwrite(&old_data_header.norm.maq.acq_date,4,1,f_aus);
      fclose(f_aus);
    }
  return 1;
}                                                       /*end funzione test_file()*/



/*===============================================*/
int elabora_dato()
/*===============================================*/
{   
  int i,l,k;
  int el_inf,el_up;
  int test_ridotto,test_vertZ, test_prec_anap,test_anap,test_nodatalow; // varibili di controllo sull'andamento del test anap
  float bin_low,bin_high, bin_low_low; 
  float fondo_scala,elevaz,quota_true;
  unsigned char flag_anap;
  float PIA;


  /*---------------------
    | calcolo fondo scala |
    ---------------------*/
  flag_anap = 1;

  fondo_scala = BYTEtoDB(flag_anap);/*-19.7 dBZ*/
  printf("leggo first level \n");
  leggo_first_level();

#ifdef QUALITY
  leggo_dem();
  leggo_hray();
#endif

/*per quando non rimuovo anap e voglio riscrivere volume polare*/
#ifdef DECLUTTER
  for(i=0; i<NUM_AZ_X_PPI; i++)
    {
      PIA=0.0
	for(k=0; k<vol_pol[0][i].b_header.max_bin; k++)
	  {
	    el_inf = first_level_static[i][k]; 	
	 
	    for(l=0; l<=el_inf; l++) 
	      {
		vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
#ifdef BEAMBLOCKING
#ifndef BLOCSTAT
		vol_pol[l][i].ray[k]=vol_pol[l][i].ray[k]+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5);
#endif
#endif
	      }
#ifdef QUALITY
	    elev_fin[i][k]=el_inf; 

#endif
	  }
    }
  return 0;
#endif  /* fine DECLUTTER */

  /*----------------------------------------
    | azzero matrici per statistica anaprop |
    ---------------------------------------*/
  memset(stat_anap_tot,0,sizeof(stat_anap_tot));
  memset(stat_anap,0,sizeof(stat_anap));
  memset(stat_bloc,0,sizeof(stat_bloc));
  memset(stat_elev,0,sizeof(stat_elev));
  printf("inizia controllo\n");

  for(i=0; i<NUM_AZ_X_PPI; i++)
    {
      flag_anap = 0;  
      PIA=0.0;
      for(k=0; k<vol_pol[0][i].b_header.max_bin; k++)
	{ 	  
	  stat_anap_tot[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
	  el_inf = first_level[i][k];
#ifdef QUALITY
	    elev_fin[i][k]=el_inf;
#endif
	  el_up = el_inf +1;
	  bin_low_low=fondo_scala+1;
	  if (el_inf>=1) bin_low_low=BYTEtoDB(vol_pol[el_inf-1][i].ray[k]);
	  bin_low  = BYTEtoDB(vol_pol[el_inf][i].ray[k]);
	  bin_high = BYTEtoDB(vol_pol[el_up][i].ray[k]);
     
	  /*------------------------------------------------
	    | correggo il test considerando che :            |
	    |   il bin inferiore deve essere superiore       |
	    |   al fondo scala;                              |
	    |   mentre quello superiore puo' essere maggiore |
	    |   o uguale al fondo scala.                     |
	    | Esempio (10., -19.7 ) ok faccio test           |
	    |         (-19.7, -19.7) non faccio test         |
	    |------------------------------------------------|
	    | vecchio codice                                 |
	    |	if(bin_low != -20. && bin_high != -20 )    |
	    ------------------------------------------------*/
	  /* 26-5-2004 : se sono alla 1 o successive elevazioni
	     e range > 80 km cambio il riconoscimento anaprop, in modo
	     da evitare di riconoscere come anaprop una pioggia a bassa quota
	     (la diff. tra Z elev. corrente e Z alla successiva puo' essere
	     grande perchè sopra non c'è niente).
	     Il criterio diventa: - se la differenza tra Z all'elevazione più bassa della
	     corrente e la Z corrente è <10 dbZ allora sono quasi certa di non avere 
	     anaprop, percio' cambio i limiti di riconoscimento anaprop, altrimenti li
	     lascio identici. */

	  MAX_DIF=MAX_DIF_OR;                        
	  MAX_DIF_NEXT=MAX_DIF_NEXT_OR;              
	  MIN_VALUE=MIN_VALUE_OR;                    
	  MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;

	  if((el_inf>=1)&&(k>LIMITE_ANAP)&&(bin_low_low-bin_low<10)) /*se sono a elev almeno 2, oltre 80km e test applicato sotto non da anaprop =>*/
	    {   
	      test_ridotto=1;
	      MAX_DIF_NEXT=MAX_DIF_LIMIT;
	      MAX_DIF=MAX_DIF_LIMIT;
	      MIN_VALUE=MIN_VALUE_LIMIT;                        
	      MIN_VALUE_NEXT=MIN_VALUE_LIMIT;     }

	  if(bin_low > fondo_scala && bin_high >= fondo_scala )// ho qualcosa sia sotto che sopra
	    {
	      test_vertZ=1;
	      if(flag_anap)
		{		  
		  test_prec_anap=1;// anaprop precedentemente riscontrata nel raggio (riduco le soglie di test)
		  if(bin_low-bin_high >= MAX_DIF_NEXT || bin_high <= MIN_VALUE_NEXT ||(bin_low_low==fondo_scala && bin_high==fondo_scala)) 
		    {
		      test_anap=1;// anaprop sì
		      for(l=0; l<el_up; l++) {
			vol_pol[l][i].ray[k]=1;	//metto a fondoscala i valori di vol_pol sul pixel fino a el_inf  (no preci nella ZLR finale)  
		      }
		      flag_anap = 1;
		      stat_anap[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++; //incremento la statitica anaprop
		      if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++; //incremento la statitica cambio elevazione
	      
#ifdef QUALITY
		      dato_corrotto[i][k]=1;/* matrice risultato test: propagazione anomala*/

#endif
#ifdef BEAMBLOCKING
		      beam_blocking[i][k]=0;/* beam blocking azzerato ( ho cambiato elevazione, non posso determinarlo)*/
#endif
	     
		    }
		  else 
		    {
		      test_anap=0;// anaprop no
		      for(l=0; l<=el_inf; l++){
			vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k]; // assegno a tutti i bin sotto el_inf il valore a el_inf (preci/Z a el_inf nella ZLR finale)
#ifdef BEAMBLOCKING
#ifndef BLOCSTAT
			vol_pol[l][i].ray[k]=vol_pol[l][i].ray[k]+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5); //correggo beam blocking	
			stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]= stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]+ beam_blocking[i][k]; // incremento statistica beam blocking
#endif
#endif
		      } 
#ifdef QUALITY		   
		      dato_corrotto[i][k]=0;/* matrice risultato test: no propagazione anomala*/	
#endif	    
		      if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;

		      flag_anap = 0;
		    }
		} 
	      else    /* se invece flag_anap == 0 cioè non ho ancora trovato anaprop nel raggio */
		{
		  test_prec_anap=0;// anaprop precedentemente non riscontrata nel raggio	
		  if(bin_low-bin_high >= MAX_DIF || bin_high <= MIN_VALUE ||(bin_low_low==fondo_scala && bin_high==fondo_scala)) 
		    {
		      test_anap=1;// anaprop sì
		      for(l=0; l<el_up; l++){
			vol_pol[l][i].ray[k]=1;
		      }
		      flag_anap = 1;
		      stat_anap[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
#ifdef QUALITY
		      dato_corrotto[i][k]=1;/*matrice risultato test: propagazione anomala*/
#endif
		      if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
#ifdef BEAMBLOCKING
		      beam_blocking[i][k]=0;
#endif
			    } 
		  else 
		    {

		      test_anap=0;// anaprop no		     
		      for(l=0; l<=el_inf; l++)
			{
			  vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
#ifdef BEAMBLOCKING
#ifndef BLOCSTAT
			  vol_pol[l][i].ray[k]=vol_pol[l][i].ray[k]+ceil(-3.1875*10.*log10(1.-(float)beam_blocking[i][k]/100.)-0.5);  
			  stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]= stat_bloc[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]+ beam_blocking[i][k];
#endif
#endif
			}

#ifdef QUALITY           
		      dato_corrotto[i][k]=0;	
#endif
		      if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
		      flag_anap = 0;
		      /*printf(" non corretta - 2\n");*/
		    } /*endif test anaprop*/
		}/*endif flaganap--anap precedente*/
	    }/*endif bin_low > fondo_scala && bin_high >= fondo_scala*/
	  else if (bin_low < fondo_scala)  
	    {
	      test_vertZ=0;
	      test_nodatalow=1;

	      for(l=0; l<el_up; l++) //riempio con i valori di el_up tutte le elevazioni sotto (ricostruisco il volume)
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
#ifdef QUALITY
	      if (bin_high<fondo_scala)   dato_corrotto[i][k]=2;/*manca dato sotto e sopra*/
	      if (bin_high>fondo_scala)  dato_corrotto[i][k]=3;/*manca controllo (sotto non ho nulla)*/
	      if (bin_high==fondo_scala) dato_corrotto[i][k]=0;/*non piove*/
              elev_fin[i][k]=el_up; 
#endif
	      if (el_up > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
#ifdef BEAMBLOCKING
	      beam_blocking[i][k]=0;
#endif   
	    }
	  else if (bin_low == fondo_scala || bin_high < fondo_scala)/* quel che resta da (bin_low > fondo_scala && bin_high >= fondo_scala) e (bin_low < fondo_scala) basterebbe un else*/
	    {
	      test_vertZ=0;
	      test_nodatalow=0;
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
#ifdef QUALITY
	      dato_corrotto[i][k]=0; //non so perchè..

#endif
	      if (el_inf > first_level_static[i][k]) stat_elev[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
	    }
 /*-----------------------------------------------------------fine-----------*/
	  /*CALIBRAZIONE*/

	  for (l=0;l<NEL;l++)
	    {
	      if (vol_pol[l][i].ray[k]>1) vol_pol[l][i].ray[k]=(DBtoBYTE(BYTEtoDB(vol_pol[l][i].ray[k])+calibr)>0)?DBtoBYTE(BYTEtoDB(vol_pol[l][i].ray[k])+calibr):1;
	    }                         

#ifdef QUALITY
	  elevaz=(float)(vol_pol[elev_fin[i][k]][i].teta_true)*CONV_RAD;	
	  quota_true=quota_f(elevaz,k);	     
	  quota_rel[i][k]=(unsigned short)(hray[k][elev_fin[i][k]]-dem[i][k]);/*quota sul suolo */  	 
	  quota[i][k]=(unsigned short)(quota_true);
#endif
	  }                     
    } 	

  printf ("fatta elabora_dato\n");
  ScrivoStatistica();

  return 0;
}                          /*end funzione elabora_dato()*/
#ifdef QUALITY
/*
  comstart caratterizzo_volume
  idx calcola qualita' volume polare
  calcola qualita' volume polare 
  NB il calcolo è fatto considerando q=0 al di sotto della mappa dinamica.
  per ora drrs=dist nche nel caso di Gattatico, mentre dtrs è letto da file
  si puo' scegliere tra qualita' rispetto a Z e rispetto a R.
 
  double PIA;  path integrated attenuation
  float dh: dimensione verticale bin calcolata tramite approcio geo-ottico
  float drrs: distanza radiosondaggio,
  float dtrs: tempo dal radiosondaggio 
  float dist: distanza dal radar
  float quo:  quota bin
  float el:   elevazione
  float rst:  raggio equivalente in condizioni standard
  float dZ:   correzione vpr
  float sdevZ:  stand. dev. cprrezione vpr
  unsigned char bb: beam blocking
  unsigned char cl: indice clutter da anaprop

  comend
*/
/*===============================================*/
void caratterizzo_volume() 
/*===============================================*/

{   
  int i,l,k;
  double PIA;
  float dh=1.,dhst=1.,drrs=1.,dist=1.,quo,el,elevaz;
  unsigned char bb=0,cl=0 ; 


  /* a questo punto servono: bb, cl,  PIA, dtrs e drrs radiosond, quota, hsup e hinf beam*/

  for (l=0; l<NSCAN; l++)/*ciclo elevazioni*/
    { 
      for (i=0; i<NUM_AZ_X_PPI; i++)/*ciclo azimuth*/
	{
	  el=l*0.9+0.5;
	  el=el*M_PI/180.;   
	  elevaz=(float)(vol_pol[l][i].teta_true)*CONV_RAD;
	  PIA=0.;
	  for (k=0; k<vol_pol[0][i].b_header.max_bin; k++)/*ciclo range*/
	    { 
	      dist= k*size_cell[old_data_header.norm.maq.resolution]+size_cell[old_data_header.norm.maq.resolution]/2.;/*distanza radar M*/
	      drrs=dist;
	      if (l == elev_fin[i][k]) att_cart[i][k]=DBtoBYTE(PIA);
	      PIA=attenuation(vol_pol[l][i].ray[k],PIA);	     
	      //dhst=dist*AMPLITUDE*M_PI/180.;
	      dhst =quota_f(elevaz+0.45*DTOR,k)-quota_f(elevaz-0.45*DTOR,k);
	      if (l<NSCAN-1   ) {
		dh=hray_inf[k][l+1]-hray_inf[k][l]; /* differenza tra limite sup e inf lobo centrale secondo appocci o geo-ott*/	
	      }
	      else {
		dh=dhst; /* non ho le altezze oltre nscan-1 pero' suppongo che a tali elevazioni la prop. si possa considerare standard*/
		hray[k][l]=(sqrt(pow(dist/1000.,2)+(rst*rst)+2.0*dist/1000.*rst*sin(el))-rst)*1000.; /*quota in prop. standard  */  	
	      } 

	      quo=quota_f(elevaz,k);  
	      if (l-elev_fin[i][k] <0) {
		qual[l][i][k]=0;  
		cl=2;
	      }
	      if (l == elev_fin[i][k] ) {
		cl=dato_corrotto[i][k];  /*cl al livello della mappa dinamica*/ 
		bb=beam_blocking[i][k];  /*bb al livello della mappa dinamica */
	      }
	      if (l-elev_fin[i][k] >0 ) {
		cl=0;	    /*per come viene scelta la mappa dinamica si suppone che al livello superiore cl=0 e bb=0*/
		bb=0;   
	      }
	      qual[l][i][k]=(unsigned char)(func_q_Z(cl,bb,drrs,dtrs,dh,dhst,PIA)*100);
	      if(qual[l][i][k]==0) qual[l][i][k]=1;
#ifdef VPR	
	      /* sezione PREPARAZIONE DATI VPR*/
	      if(cl==0 && bb<BBMAX_VPR )   /*pongo le condizioni per individuare l'area visibile per calcolo VPR*/		 
		flag_vpr[l][i][k]=1;	 
#endif                                                    
	    }

	}
    }
 
  return;
} 
#endif
/*===============================================*/
double attenuation(unsigned char DBZbyte, double  PIA)  /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/
/*===============================================*/

{ 

  double Zhh,att_rate,R;/* PIA diventa att_tot devo decidere infatti se PIA sarà 3d percio' temp. uso  nomi diversi*/
  double att_tot;
  
  
  att_tot=PIA;
  Zhh=(double)(BYTEtoZ(DBZbyte));
  if (10*log10(Zhh) > THRES_ATT )
    {   
      Zhh=pow(10., (log10(Zhh)+ 0.1*att_tot));
      R=pow((Zhh/aMP),(1.0/bMP));
      att_rate=0.0018*pow(R,1.05);
      att_tot=att_tot+2.*att_rate*0.001*size_cell[old_data_header.norm.maq.resolution];
      if (att_tot>BYTEtoDB(254)) att_tot=BYTEtoDB(254);
    }
  return att_tot;
} 



/*==============================================*/
void ScrivoStatistica()
/*==============================================*/
{
  FILE *f_stat;
  int az,ran;
  unsigned char statistica[DIM1_ST][DIM2_ST];
  unsigned char statistica_bl[DIM1_ST][DIM2_ST];
  unsigned char statistica_el[DIM1_ST][DIM2_ST];
  struct tm *tempo;
  time_t Time;
  char date[20];

  printf ("scrivo statistica ");
  memset(statistica,255,DIM1_ST*DIM2_ST);
  memset(statistica_bl,255,DIM1_ST*DIM2_ST);
  memset(statistica_el,255,DIM1_ST*DIM2_ST);
  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&Time);
  sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
	  tempo->tm_mday,tempo->tm_hour, tempo->tm_min);

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
  f_stat = fopen(getenv("ANAP_STAT_FILE"),"a");
  fwrite(date,12,1,f_stat);
  fwrite(statistica,DIM1_ST*DIM2_ST,1,f_stat);
  fclose(f_stat);

  f_stat = fopen(getenv("BLOC_STAT_FILE"),"a");
  fwrite(date,12,1,f_stat);
  fwrite(statistica_bl,DIM1_ST*DIM2_ST,1,f_stat);
  fclose(f_stat);

  f_stat = fopen(getenv("ELEV_STAT_FILE"),"a");
  fwrite(date,12,1,f_stat);
  fwrite(statistica_el,DIM1_ST*DIM2_ST,1,f_stat);
  fclose(f_stat);


  return ;
}


/*===============================================*/
int write_dbp(nome)
/*===============================================*/
     char *nome;
{
  int i,j,size_h;
  FILE *file1;
  
  size_h = sizeof(vol_pol[0][0].b_header);
  file1=fopen(nome,"w");
  if(file1 == NULL) return 1;
  
  if(fwrite(&old_data_header,sizeof(old_data_header),1,file1)== 0) return 2;
  
  for (i=0; i<NEL; i++)
    for (j=0; j<NUM_AZ_X_PPI; j++)
      {
	if(vol_pol[i][j].flag)
	  {
	    if(fwrite(&vol_pol[i][j].b_header,size_h,1,file1)== 0) return 3;
	    if(fwrite(vol_pol[i][j].ray,vol_pol[i][j].b_header.max_bin,1,file1) ==
	       0 ) return 4;
	  }
      }
  fclose(file1);
  return 0;
}


/*===============================================*/
void   leggo_first_level()
/*===============================================*/
{
  FILE *file;
  int i,j,dim;
  char first_level_bb_file[200],bb_value_file[200];
  struct tm *tempo;
  time_t time;

#ifdef STATIC
  /*-------------------
    Leggo mappa statica
    -------------------*/
    
  file=controllo_apertura(getenv("FIRST_LEVEL_DIM_FILE")," dimensioni mappa statica","r");
 
  fscanf(file,"%i ",&dim);
  fclose(file);
  file=controllo_apertura(nome_fl," mappa statica","r");
  for(i=0; i<NUM_AZ_X_PPI; i++)
    fread(&first_level_static[i][0],dim,1,file);
  memcpy(first_level,first_level_static,sizeof(first_level));	   
  fclose(file);
  printf ("letta mappa statica \n");
#endif

#ifdef BEAMBLOCKING
  /*----------------------------
    Leggo file elevazioni per BB  
    ----------------------------*/

#ifdef SHORT
  time = NormalizzoData(old_data_header.norm.maq.acq_date); /* arrotonda ai 5 minuti di precisione*/
  tempo = gmtime(&time);
#endif
#ifdef MEDIUM
  tempo = gmtime(&old_data_header.norm.maq.acq_date);
  time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&time);
#endif
        
  sprintf(first_level_bb_file,"%s/%04d%02d%02d%02d%02dmat_el.bin",
	  getenv("DIR_OUT_PP_BLOC"),
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);
  file=controllo_apertura(first_level_bb_file," elev BB ","r");	
  for(i=0; i<NUM_AZ_X_PPI; i++)
    fread(&bb_first_level[i][0],MAX_BIN,1,file);
  fclose(file);
  printf ("letta mappa elevazioni da prog beam blocking\n");
  /*------------------------
    Leggo file valore di BB  
    ------------------------*/
  sprintf(bb_value_file,"%s/%04d%02d%02d%02d%02dmat_bloc.bin",
	  getenv("DIR_OUT_PP_BLOC"),
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

  file=controllo_apertura(bb_value_file," elev BB ","r");

  /* Se elevazione clutter statico < elevazione BB, prendi elevazione BB,
     altrimeti prendi elevazione clutter statico e metti a 0 il valore di BB*/
  for(i=0; i<NUM_AZ_X_PPI; i++){   /*ciclo sugli azimut*/
    fread(&beam_blocking[i][0],MAX_BIN,1,file);

    for (j=0; j<MAX_BIN; j++) /*ciclo sul range  */      
      {  
#ifndef BLOCSTAT     
	if (first_level_static[i][j]<=bb_first_level[i][j])
	  first_level[i][j]=bb_first_level[i][j]; 
	else
	  {  beam_blocking[i][j]=0;
	    first_level[i][j]=first_level_static[i][j]; } 
# endif
#ifdef BLOCSTAT
	if (first_level_static[i][j]>bb_first_level[i][j])
	  beam_blocking[i][j]=0; 
        if (first_level_static[i][j]<bb_first_level[i][j])
	  beam_blocking[i][j]=OVERBLOCKING; 
#endif
      }
  }
  fclose(file);
  printf ("letta mappa beam blocking \n");
#endif

  /*-------------------------------
    patch per espandere il clutter
    -------------------------------*/
#ifdef MEDIUM
  {
    unsigned char first_level_tmp[NUM_AZ_X_PPI][MAX_BIN];
    int k;
    memcpy(first_level_tmp,first_level,sizeof(first_level));
    for (i=NUM_AZ_X_PPI; i<800; i++)
      {
	for (j=0; j<MAX_BIN; j++)
	  {
	    for (k=i-1; k<i+2; k++)
	      if(first_level[i%NUM_AZ_X_PPI][j] < first_level_tmp[k%NUM_AZ_X_PPI][j])
		first_level[i%NUM_AZ_X_PPI][j]=first_level_tmp[k%NUM_AZ_X_PPI][j];
	  }
      }
  }
#endif
  printf ("fine lettura mappe elevazioni \n");

	
}


/*===============================================*/
void creo_cart()
/*===============================================*/
{
  int i,j,quad,x,y,irange,az,iaz,az_min,az_max,cont;
  static int flag = 1;

  if(flag)
    {
      creo_matrice_conv();
      flag = 0;
    }

  for(i=0; i<MAX_BIN *2; i++)
    for(j=0; j<MAX_BIN *2; j++)
      cart[i][j] = MISSING;

  /*
    Stampa di controllo
    for(i=0; i<MAX_BIN*2; i++)
    for(j=MAX_BIN-10; j<MAX_BIN+10; j++)
    if (i<10 || i> MAX_BIN*2-10)
    printf("%3d  %3d  %3d\n",i,j,cart[i][j]);
  */

  for(quad=0; quad<4; quad++)
    for(i=0; i<MAX_BIN; i++)
      for(j=0; j<MAX_BIN; j++)
	{
	  irange = (int)range[i][j];
	  if(range[i][j] - irange >= .5) irange++;
	  if(irange < MAX_BIN)	      {
	    switch(quad)
	      {
	      case 0:
		x = MAX_BIN + i;
		y = MAX_BIN + j;
		az = azimut[i][j];
		break;
	      case 1:
		x = MAX_BIN + j;
		y = MAX_BIN - i;
		az = azimut[i][j] + 90.;
		break;
	      case 2:
		x = MAX_BIN - i;
		y = MAX_BIN - j;
		az = azimut[i][j] + 180.;
		break;
	      case 3:
		x = MAX_BIN - j;
		y = MAX_BIN + i;
		az = azimut[i][j]+270.;
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
	      if(cart[x][y]<=vol_pol[0][iaz%NUM_AZ_X_PPI].ray[irange]){
		cart[x][y] = vol_pol[0][iaz%NUM_AZ_X_PPI].ray[irange];
#ifdef QUALITY
		qual_Z_cart[x][y]=qual[elev_fin[iaz%NUM_AZ_X_PPI][irange]][iaz%NUM_AZ_X_PPI][irange];
		quota_cart[x][y]=quota[iaz%NUM_AZ_X_PPI][irange];
		/*neve_cart[x][y]=qual_Z_cart[x][y];*/
#ifdef VPR
		neve_cart[x][y]=(neve[iaz%NUM_AZ_X_PPI][irange])?0:1;  
		corr_cart[x][y]=corr_polar[iaz%NUM_AZ_X_PPI][irange];
#endif
#endif
#ifdef CLASS
		if (irange<x_size)
		  conv_cart[x][y]=conv[iaz%NUM_AZ_X_PPI][irange];
		cappi_cart[x][y]=cappi[iaz%NUM_AZ_X_PPI][irange];
#endif
	      }
#ifdef ZLR_MEDIA
	      if (vol_pol[0][iaz%NUM_AZ_X_PPI].ray[irange] > 0){
		cartm[x][y]=cartm[x][y]+BYTEtoZ(vol_pol[0][iaz%NUM_AZ_X_PPI].ray[irange]);
		cont=cont+1;
	      } 
#endif		
	    }
#ifdef ZLR_MEDIA
	    if (cont > 0) cartm[x][y]=cartm[x][y]/(float)(cont);
#endif	
	    /*  
	 *****  per scrivere in griglia cartesiana************
	 bloc_xy[MAX_BIN*2-y][x]=beam_blocking[(int)((float)(az)/.9)][irange];
	 elev_xy[MAX_BIN*2-y][x]=elev_fin[(int)((float)(az)/.9)][irange];  
	 dato_corrotto_xy[MAX_BIN*2-y][x]= dato_corrotto[(int)((float)(az)/.9)][irange];
	 elev_fin_xy[MAX_BIN*2-y][x]=first_level[(int)((float)(az)/.9)][irange]; 
	 dato_corrotto_xy[MAX_BIN*2-y][x]= dato_corrotto[(int)((float)(az)/.9)][irange]; */
	  }
	}
}


/*===============================================*/
void 	creo_matrice_conv()
/*===============================================*/
{
  int i,j;
        

  for(i=0; i<MAX_BIN; i++)
    for(j=0; j<MAX_BIN; j++)
      {
	range[i][j] = hypot(i+.5,j+.5);
	azimut[i][j] = 90. - atan((j+.5)/(i+.5)) * M_1_PI*180.;
      }
  return;
}




/*===============================================*/
time_t 	NormalizzoData(time)    /* time è in secondi, itime è un intero che rappresenta il numero intero di intervalli da 5 minuti*/
/*===============================================*/
     time_t time;
{
  int itime;

  itime = time/(NMIN*60);

  /*
    printf(" esco da Normalizzo %d %d %d \n",time,itime,(time - itime*NMIN*60));
    printf("%s\n",ctime(&time));
    printf("%d\n",(NMIN-MAX_TIME_DIFF)*60);
  */
  if(time - itime*NMIN*60 <MAX_TIME_DIFF*60) return (itime*NMIN*60); /* se la differenza è meno di tre minuti vado al 5° min. prec*/
  if(time - itime*NMIN*60 >(NMIN-MAX_TIME_DIFF)*60) return ((itime+1)*NMIN*60); /* se la differenza è meno di tre minuti vado al 5° min. successivo*/
  return -1;
}


/*===============================================*/
void bacini()
/*===============================================*/
{

  short tabnew[1000][1000];
  float bacini[500],noss_bacini[500];
  FILE *fp, *f_aus;
  struct tm *tempo;
  time_t Time;
  char date[150],dateold[150];
  int i,j,ind,flag;
  static int Piove = NON_PIOVE;
  struct {
    time_t last_ist,last_rain;
    char name[256];
  } aus_file;

  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&Time);
  sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
	  tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
  fp=fopen(getenv("MATRICE_BACINI"),"r");
  fread(tabnew,sizeof(tabnew),1,fp);

  fclose(fp);

  /*---------------------------
    |  azzero vettori di lavoro |
    ---------------------------*/
  memset(bacini,0,sizeof(bacini));
  memset(noss_bacini,0,sizeof(noss_bacini));

  /*---------------------------
    |  ciclo su tutto il cappi  |
    ---------------------------*/
  for (i=0;i<1000; i++)
    for(j=0; j<1000; j++)
      {
	/*---        ind = tabnew[i][j];  ---*/
        ind = tabnew[j][i];
        if(ind != -1 && cart[i+12][j+12] != MISSING)
	  {
	    bacini[ind] +=   BYTE_to_mp_func(cart[12+i][12+j],aMP,bMP);
	    if(bacini[ind] <0. )
	      {
		sprintf(errori,"ERRORE GRAVE Bacino con Precipitzaione Negativa");
		ScrivoLog(16,errori);	
		sprintf(errori,"Coordinate Punto %4d %4d  Precipitazione %f",i,j,
			BYTE_to_mp_func(cart[12+i][12+j],aMP,bMP));
		ScrivoLog(16,errori);	
	      }
	    noss_bacini[ind]++;
	  }
      }
  for(i=0; i<500; i++)
    if(noss_bacini[i] != 0)
      bacini[i] = bacini[i]/noss_bacini[i];
    else 
      bacini[i] = -1.;

  /*-----------------------------------
    | scrivo file output in formato xdr |
    -----------------------------------*/
  write_xdr(bacini,Time);

  /*----------------
    | testo se piove |
    ----------------*/
  flag = NON_PIOVE;
  for(i=0; i<500; i++)
    if(bacini[i] >= RAIN_MIN)
      flag = PIOVE;

  if(flag == PIOVE) printf(" data : %s PIOVE\n",date );
  else 
    printf(" data : %s NON PIOVE \n",date);

  switch(flag)
    {
    case PIOVE:
      if(access(getenv("BACINI_AUS_FILE"),F_OK) == 0)
	{
	  Piove = PIOVE;
	  f_aus = fopen(getenv("BACINI_AUS_FILE"),"r+");
	  fread(&aus_file,sizeof(aus_file),1,f_aus);
	  /*---------------------------------
	    | testo se presente buco nei dati |
	    ---------------------------------*/
	  /*-----------------------------------------
	    | Se esite aggiorno il file History       |
	    | con la chiusura dell' evento precedente |
	    | identificando come data finale          |
	    | l'ultimo record presente                |
	    | ridefinisco il nome del file            |
	    | contenente l'evento e assegno alla      |
	    | flag Piove il valore NON_PIOVE          |
	    -----------------------------------------*/
	  if(old_data_header.norm.maq.acq_date-aus_file.last_ist > DIFF_TIME)
	    { /*-- esiste un buco temporale nei dati --*/
	      sprintf(errori,"Chiudo il File di Hystory per Buco");
	      ScrivoLog(16,errori);	
	      Time = NormalizzoData(aus_file.last_ist);
	      tempo = gmtime(&Time);
	      sprintf(dateold,"%04d%02d%02d%02d%02d",
		      tempo->tm_year+1900, tempo->tm_mon+1,
		      tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
	      AggiornoHistoryFile(dateold,NON_PIOVE);
	      Time = NormalizzoData(old_data_header.norm.maq.acq_date);
	      tempo = gmtime(&Time);
	      sprintf(aus_file.name,"%s/%04d%02d%02d%02d%02d", 
		      getenv("BACINI_DIR"),
		      tempo->tm_year+1900, tempo->tm_mon+1,tempo->tm_mday,
		      tempo->tm_hour, tempo->tm_min);
	      Piove = NON_PIOVE;
	    }
	  aus_file.last_rain = old_data_header.norm.maq.acq_date;
	  aus_file.last_ist  = old_data_header.norm.maq.acq_date;
	  rewind (f_aus);
	  fwrite(&aus_file,sizeof(aus_file),1,f_aus);
	  fclose(f_aus);
	  /* devo scrivere */
	}
      else
	{ /*-- non esiste il file BACINI_AUS_FILE  --*/

	  f_aus = fopen(getenv("BACINI_AUS_FILE"),"w");
	  aus_file.last_rain = old_data_header.norm.maq.acq_date;
	  aus_file.last_ist  = old_data_header.norm.maq.acq_date;
	  sprintf(aus_file.name,"%s/%04d%02d%02d%02d%02d",
		  getenv("BACINI_DIR"),
		  tempo->tm_year+1900, tempo->tm_mon+1,tempo->tm_mday,
		  tempo->tm_hour, tempo->tm_min);
	  fwrite(&aus_file,sizeof(aus_file),1,f_aus);
	  fclose(f_aus);
	}
      ScrivoBacini(aus_file.name,date,bacini,466);
      if(Piove == NON_PIOVE)
	{
	  AggiornoHistoryFile(date,PIOVE);
	  Piove = PIOVE;
	}
      break;
    case NON_PIOVE:
      if(access(getenv("BACINI_AUS_FILE"),F_OK) != 0) break;
      f_aus = fopen(getenv("BACINI_AUS_FILE"),"r+");
      fread(&aus_file,sizeof(aus_file),1,f_aus);
      /*----------------------------------------------
	| e' presente un buco nei dati chiudo l'evento |
	----------------------------------------------*/
      if(old_data_header.norm.maq.acq_date-aus_file.last_ist > DIFF_TIME) 
	{
	  Time = NormalizzoData(aus_file.last_ist);
	  tempo = gmtime(&Time);
	  sprintf(dateold,"%04d%02d%02d%02d%02d",
		  tempo->tm_year+1900, tempo->tm_mon+1,
		  tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
	  Piove = NON_PIOVE;
	  AggiornoHistoryFile(dateold,NON_PIOVE);
	  fclose(f_aus);
	  unlink(getenv("BACINI_AUS_FILE"));
	  break;
	}
      /*---------------------------------
	| non piove da  DIFF_TIME secondi |
	---------------------------------*/
      else if(old_data_header.norm.maq.acq_date-aus_file.last_rain > DIFF_TIME)
	{
	  fclose(f_aus);
	  unlink(getenv("BACINI_AUS_FILE"));
	  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
	  tempo = gmtime(&Time);
	  sprintf(dateold,"%04d%02d%02d%02d%02d",
		  tempo->tm_year+1900, tempo->tm_mon+1,
		  tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
	  Piove = NON_PIOVE;
	  AggiornoHistoryFile(dateold,NON_PIOVE);
	  break;
	}
      aus_file.last_ist  = old_data_header.norm.maq.acq_date;
      rewind (f_aus);
      fwrite(&aus_file,sizeof(aus_file),1,f_aus);
      fclose(f_aus);
      ScrivoBacini(aus_file.name,date,bacini,466);
      break;
    }
  return;
}                   /*end funzione bacini()*/



/*===============================================*/
void ScrivoBacini(nome,data,array_bacini,quanti)
/*===============================================*/
     char *nome,*data;
     float *array_bacini;
     int quanti;
{
  FILE  *f_out;
  /*--------------------------
    |  scrivo i dati in output |
    --------------------------*/
  if(access(nome,F_OK) == 0 )
    f_out = fopen(nome,"a");
  else
    f_out = fopen(nome,"w");

  fwrite(data,strlen(data),1,f_out);
  fwrite(array_bacini,quanti*4,1,f_out);
  fclose(f_out);
  return;
}


/*===============================================*/
void prendo_tempo()
/*===============================================*/
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



/*===============================================*/
void AggiornoHistoryFile(data,flag)
/*===============================================*/
     char *data;
     int flag;
{
  FILE *fhistory;
  int ntot = 0;
  char data2[] = {"NNNNNNNNNNNN"};
  char scratch[80];

  if(access(getenv("BACINI_HISTORY_FILE"),6) != 0)
    /*--- il file di history non esiste ---*/
    {
      sprintf(errori,"NON ESISTE IL FILE : %s",getenv("BACINI_HISTORY_FILE"));
      ScrivoLog(16,errori);	
      printf("%s\n",strerror(errno));
      switch (flag)
	{
	case INIZIO_EVENTO:
	  fhistory = fopen(getenv("BACINI_HISTORY_FILE"),"w");
	  ntot = ntot+1;
	  fprintf(fhistory," Sono presenti %05d eventi\n",ntot);
	  fseek(fhistory,0,2);
	  fprintf(fhistory," Inizio: %s  Fine: %s\n",data,data2);
	  fclose(fhistory);
	  break;
	case FINE_EVENTO :
	  printf(" errore nel file di history\n");
	  break; 
	}
    }
  else
    /*--- il file di history esiste ---*/
    {
      sprintf(errori,"ESISTE IL FILE : %s",getenv("BACINI_HISTORY_FILE"));
      ScrivoLog(16,errori);	
      switch (flag)
	{
	case INIZIO_EVENTO:
	  fhistory = fopen(getenv("BACINI_HISTORY_FILE"),"r+");
	  fscanf(fhistory,"%14c%6d%8c",scratch,&ntot,scratch);
	  rewind(fhistory);
	  ntot = ntot+1;
	  fprintf(fhistory," Sono presenti %05d eventi\n",ntot);
	  fseek(fhistory,0,2);
	  fprintf(fhistory," Inizio: %s  Fine: %s\n",data,data2);
	  fclose(fhistory);
	  break;
	case FINE_EVENTO :
	  fhistory = fopen(getenv("BACINI_HISTORY_FILE"),"r+");
	  fscanf(fhistory,"%14c%6d%8c",scratch,&ntot,scratch);
	  fseek(fhistory,-42,2);
	  fscanf(fhistory,"%9c%12c",scratch,data2);
	  data2[12]=0;
	  sprintf(errori,"CASO DI FINE EVENTO : %s",data2);
	  ScrivoLog(16,errori);	
	  fseek(fhistory,-42,2);
	  fprintf(fhistory," Inizio: %s  Fine: %s\n",data2,data);
	  fclose(fhistory);
	  break;
	}
    }
}                   /*end funzione AggiornoHistoryFile(data,flag)*/



/*===============================================*/
void	creo_cart_z_lowris()
/*===============================================*/
{
  int i,j,x,y,cont;
  unsigned char z,q,nv,c1x1;
  unsigned short q1x1;
  float zm;
        

  memset(z_out,0,sizeof(z_out));
 
  for(i=0; i<CART_DIM_ZLR; i++)
    for(j=0; j<CART_DIM_ZLR; j++)
      {
	z = 0;
	q = 0;
	zm = 0.;
	q1x1=0;
	c1x1=0;
	cont=0;
	for(x = 0; x < ZLR_N_ELEMENTARY_PIXEL; x++)
	  for(y = 0; y < ZLR_N_ELEMENTARY_PIXEL; y++)
	    {
	      if(cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] != MISSING)
		if(cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] > z){ 
		  z= cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
#ifdef QUALITY
		  q= qual_Z_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
 		  q1x1=quota_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET]; 

#ifdef VPR
		  c1x1=corr_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET]; 
		  nv= neve_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET]; 
#endif	  
#endif 	 
#ifdef CLASS
		  conv_1x1[i][j]=conv_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
		  //  if (conv_1x1[i][j] >= CONV_VAL && BYTEtoDB(z)>THRES_ALLERT1) conv_1x1[i][j]=CONV_VAL+80;
		  // if (conv_1x1[i][j] >= CONV_VAL && BYTEtoDB(z)>THRES_ALLERT2) conv_1x1[i][j]=CONV_VAL+150;
		  cappi_1x1[i][j]=cappi_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
#endif    

#ifdef ZLR_MEDIA
		  if (cartm[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] > 0) {
		    zm = zm + cartm[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
		    cont=cont+1;
		  }
#endif
		}    
	      z_out[i][j]=z;
	      qual_Z_1x1[i][j]=q;

	      quota_1x1[i][j]=128+(unsigned char)(q1x1/100);	
#ifdef VPR
	      neve_1x1[i][j]=nv;     
	      corr_1x1[i][j]=c1x1;
#endif 	 
#ifdef ZLR_MEDIA
	      if (cont >0 ) {
		z_out[i][j]=(unsigned char)((10*log10(zm/(float)(cont))+20.)/80.*255); 
	      }
	      if (cont =0 ) z_out[i][j]=MISSING; 	  
#endif

	    }
      }
}

/*=======================================================================================*/
void scrivo_out_file_bin (char *ext,char *content,char *dir,size_t size, void  *matrice)
/*=======================================================================================*/
/* scrive in output matrice di byte di dimensione size*/
{
  char nome_file [150];
  FILE *output;
  struct tm *tempo;
  time_t time;
  /*----------------------------------------------------------------------------*/
  /*	apertura file dati di output	                                      */
  /*----------------------------------------------------------------------------*/

#ifdef SHORT
  time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&time);
#endif
#ifdef MEDIUM
  tempo = gmtime(&old_data_header.norm.maq.acq_date);
  time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&time);
#endif
  sprintf(nome_file,"%s/%04d%02d%02d%02d%02d%s",dir,
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min, ext);

  output=controllo_apertura(&nome_file,content,"w");
  printf("aperto file %s dimensione matrice %d\n",nome_file,size);
  fwrite(matrice,size,1,output);
  fclose(output);
  return;
}


/*=======================================================================================*/
FILE *controllo_apertura (char *nome_file, char *content,char *mode)
/*=======================================================================================*/
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
      sprintf(errori,"Errore Apertura %s %s",content, nome_file);
    
      ScrivoLog(16,errori);	
      exit(1);
    }
  if (file == NULL) {
      sprintf(errori,"Errore Apertura %s %s",content, nome_file);
     
      ScrivoLog(16,errori);	
      exit(1);    
  }
  return(file);
}

/*=======================================================================================*/
void ScrivoLog(i,stringa)
/*=======================================================================================*/
     int i;
     char *stringa;
{
  static FILE *log;
  switch (i)
    {
    case  0:
      log=fopen(getenv("LOG_FILE"),"a");
      if(log==NULL) 
	{
	  printf(" impossibile aprire il file di log \n");
	  exit(1);
	}
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
#ifdef BACINI
      fprintf(log," BACINI");
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
      fprintf(log,"MATRICE_BACINI = %s\n",getenv("MATRICE_BACINI"));
      fprintf(log,"BACINI_AUS_FILE = %s\n",getenv("BACINI_AUS_FILE"));
      fprintf(log,"BACINI_DIR = %s\n",getenv("BACINI_DIR"));
      fprintf(log,"BACINI_HISTORY_FILE = %s\n",getenv("BACINI_HISTORY_FILE"));
      fprintf(log,"OUTPUT_Z_LOWRIS_DIR = %s\n",getenv("OUTPUT_Z_LOWRIS_DIR"));
      fprintf(log,"BACINI_HISTORY_FILE = %s\n",getenv("BACINI_HISTORY_FILE"));
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
    case 10:
      fprintf(log,"%s -- Estrazione Precipitazione sui Bacini\n",PrendiOra());
      break;
    case 11:
      fprintf(log,"%s -- Estrazione Precipitazione 5X5\n",PrendiOra());
      break;
    case 12:
      fprintf(log,"%s -- Scrittura Precipitazione 5X5  %s\n",PrendiOra(),stringa);
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


/*===============================================*/
char *PrendiOra()
/*===============================================*/
{
  time_t clock;
  struct tm *tempo;         

  clock=time(0);

  tempo=gmtime(&clock);

  return   asctime(tempo);

}



/*===============================================*/
void   write_xdr(bacini,Time)
/*===============================================*/
     time_t Time;
     float *bacini;
{
  FILE *f_output;
  XDR xdr;

  f_output=fopen(getenv("BACINI_BOLOGNA"),"a");
  xdrstdio_create(&xdr ,f_output,XDR_ENCODE);
  xdr_long(&xdr,&Time);
  /*  if(xdr_vector(&xdr,bacini,500,4,xdr_float) == 1)*/
  xdr_destroy(&xdr);
  fclose(f_output);
  return;
}

#ifdef QUALITY

/*===============================================*/
void   leggo_hray( )
/*===============================================*/
{
  struct tm *tempo;
  time_t time;
  FILE *file;
  int i,j;
  char nome_file_hray [150];

  /*--------------------------
    Leggo quota centro fascio
    --------------------------*/
#ifdef SHORT
  time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&time);
#endif
#ifdef MEDIUM
  tempo = gmtime(&old_data_header.norm.maq.acq_date);
  time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&time);
#endif
  sprintf(nome_file_hray,"%s/%04d%02d%02d%02d%02dh_ray.txt",
	  getenv("DIR_OUT_PP_BLOC"),
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

  file=controllo_apertura(nome_file_hray,"File hray","r");
  fscanf(file,"%f ",&dtrs);
      
  for(i=0; i<MAX_BIN; i++){
    for(j=0; j<NSCAN;j++) 
      fscanf(file,"%f ",&hray[i][j]);
  }         
  fclose(file);
  sprintf(nome_file_hray,"%s/%04d%02d%02d%02d%02dh_rayinf.txt",
	  getenv("DIR_OUT_PP_BLOC"),
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

  file=controllo_apertura(nome_file_hray,"File hray inf","r");
  fscanf(file,"%f ",&dtrs);
  for(i=0; i<MAX_BIN; i++){
    for(j=0; j<NSCAN;j++) 
      fscanf(file,"%f ",&hray_inf[i][j]);
  }
  fclose(file);
  sprintf(errori,"lette hray\n");

  return  ;
}
/*===============================================*/
void   leggo_dem()              
/*===============================================*/
{
  FILE *file;
  int i,j;
  /*---------------------
    Leggo dem
    ---------------------*/
  //        file_dem=controllo_apertura(getenv("FILE_DEM"),"File dem","r");
       
  file=controllo_apertura(nome_dem,"File dem","r");
  for(i=0; i<MAX_BIN; i++){
    for(j=0; j<NUM_AZ_X_PPI;j++)   
      fscanf(file,"%f ",&dem[j][i]); 
  }
  fclose(file);
  printf("letto dem \n");
  return ;
}
#endif


#ifdef VPR

/*===============================================*/
/*
  comstart func_vpr
  idx funzione calcolo VPR istantaneo
  calcola il VPR istantaneo secondo il metodo di Germann e Joss (2003)
  Per il calcolo si considerano i punti con Z>THR_VPR, qualità>QMIN_VPR, BeamBlocking<20% e clutter free all'interno del volume scelto. 
  Il punto del vpr corrispondente ad un livello è calcolato come la somma dei volumi di pioggia di ciascuna cella polare presente nell'intervallo di quota, diviso l'area totale. Per ogni livello devo avere un'area minima coperta di pioggia,per i primi due livelli, se non si verifica, ricopio il valore del primo livello che raggiunge questa copertura.
  La parte del profilo alta, che sta cioè al di sopra dell'ultimo livello calcolato, viene riempite tirando una retta con coeff. amgolre negativo e costante, pari a -0.03 (valore passibile di discussione)
  Il profilo è poi soggetto a quality check e viene rigettato (return(1)) se:
  - estensione verticale troppo bassa
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
  cappi1100[1200*MAX_BIN]: vettore volumi buoni tra 1000 e 1200 m
  long int area_vpr[NMAXLAYER],areaqual[NMAXLAYER]; area totale usata per calcolo vpr; area totale pesata con le qualità(per ora non usata)

  comend
*/

/*===============================================*/
int func_vpr(long int *cv, long int *ct, float vpr1[], long int area_vpr[], char *sito)
//int func_vpr(long int *cv, long int *ct, float *vpr1, long int *area_vpr, char *sito)
/*===============================================*/


{
 
 
  int l,i,iA,k,ilay,il,ilast,iaz_min,iaz_max,icounter,naz,nra,strat;
  long int dist,vert_ext,vol_rain,cappi1100[1200*MAX_BIN];
  long int areaqual[NMAXLAYER];
  //float noval,area,stdev,somma,media;
  float noval,area,somma,media;
 /*  printf("1 - %d %d - ",sizeof(*cv),sizeof(*ct));
  printf("2 - %d ",sizeof(vpr1));
  printf("3 - %d ",sizeof(area_vpr));
  printf("4 - %d ",sizeof(*sito));
  printf("%s\n",sito);
  return 1;
*/
  somma=0;
  stdev = -1.;
  for (i=0;i<NMAXLAYER; i++ ) {
areaqual[i]=(long int )(0);

  }
  for (i=0;i<1200*MAX_BIN; i++ ) cappi1100[i]=(long int )(0);

  noval=NODATAVPR;
  icounter=0;

  /*SECTION A: cv and ct retrieval and calculation of current area- weighted vpr*/
 
  *cv=0;
  *ct=0;
  vert_ext=0;
  if (!(strcmp(sito,"SPC"))){  /* limiti azimuth per calcolo VPR sarebbe bello fare un CASE ma con stringhe non si puo'!*/
    iaz_min=IAZ_MIN_SPC;
    iaz_max=IAZ_MAX_SPC; }
  if (!(strcmp(sito,"GAT"))){
    iaz_min=IAZ_MIN_GAT;
    iaz_max=IAZ_MAX_GAT;}
  if (strcmp(sito,"SPC") && strcmp(sito,"GAT")) {
    fprintf(log_vpr,"errore, sito sconosciuto o non definito");
    exit (1);
  }
  naz=iaz_max-iaz_min;
  nra=(RMAX_VPR-RMIN_VPR)/size_cell[old_data_header.norm.maq.resolution];
  for (l=0; l<NEL; l++)//ciclo elevazioni
    { 
      for (k=0; k<MAX_BIN; k++)/*ciclo range*/
	{
	  dist=k*(int)(size_cell[old_data_header.norm.maq.resolution])+(int)(size_cell[old_data_header.norm.maq.resolution])/2.;
	   ilay=floor(hray[k][l]/TCK_VPR); // poichè in teoria quota indipendente da azimuth  
	   if (ilay <0 || ilay >= NMAXLAYER) {
	     fprintf(log_vpr,"ilay %d errore\n",ilay);
	     break;}
	  for (iA=iaz_min; iA<iaz_max; iA++)//ciclo azimuth         
	    { 
	      i=(iA+NUM_AZ_X_PPI)%NUM_AZ_X_PPI;
	      vol_rain=0;
	      if (dist<RMIN_VPR || dist> RMAX_VPR || abs(vol_pol[l][iA].teta_true - vol_pol[l][iA].b_header.teta) > 1)  
		flag_vpr[l][i][k]=0;
	      area=dist*AMPLITUDE*M_PI/180.*size_cell[old_data_header.norm.maq.resolution]*flag_vpr[l][i][k]/1000.; // divido per  mille per evitare nr troppo esagerato
	      *cv=*cv+(long int)(area);
	      /*definisco area precipitante*/
	      strat=1;
	      // if (BYTEtoDB(vol_pol[l][i].ray[k])> THR_VPR && qual[l][i][k]>QMIN_VPR && flag_vpr[l][i][k]>0 && strat)
	      if (BYTEtoDB(vol_pol[l][i].ray[k])> THR_VPR && qual[l][i][k]>QMIN_VPR && flag_vpr[l][i][k]>0 )
		{   
           	 
		  vol_rain=(long int)(BYTE_to_mp_func(vol_pol[l][i].ray[k],aMP,bMP)*area);//peso ogni cella con la sua area

		  *ct=*ct+(long int)(area);               

		  if (area_vpr[ilay]> 0) vpr1[ilay]=vpr1[ilay]+(float)(vol_rain); 
		  else vpr1[ilay]=(float)(vol_rain); 
		  area_vpr[ilay]=area_vpr[ilay]+area;  
		  /*	  ---------------------------------*/
		  areaqual[ilay]=areaqual[ilay]+(long int)((area)*qual[l][i][k]); 
	       
		   if (abs(hray[k][l]-1100)<TCK_VPR/2) {	 
		     //if (abs(hray[k][ind_ray]-1100)<TCK_VPR/2) {	    
		    cappi1100[icounter]=vol_rain;    
		    icounter=icounter+1;
		  }  
		}
	    }
	}
    }
  
  fprintf(log_vpr,"calcolati ct e cv ct= %li cv= %li\n",*ct,*cv);

  /*SECTION B: vpr quality checks and re-normalisation of vpr*/
  if ((*ct) > CT_MIN*(*cv)) { 
   
 
    ilast=0;  
    vert_ext=0;
    for  (ilay=0; ilay<NMAXLAYER; ilay++){ 
      fprintf(log_vpr,"  ilay %d area_vpr= %ld  ct= %ld  cv= %ld \n", ilay, area_vpr[ilay],*ct,*cv );
      if (area_vpr[ilay]>=MIN_AREA) {
	vert_ext=vert_ext+TCK_VPR;
	vpr1[ilay]=vpr1[ilay]/(float)(area_vpr[ilay]);
	/*vpr1[ilay]=vpr1[ilay]/(float)(areaqual[ilay]);*/
	// start=1;  
	ilast=ilay;
      }
      else{

	if (ilast > 0 && vert_ext>VEXTMIN_VPR){  
	  fprintf(log_vpr,"raggiunta cima profilo \n");
	  ilast=ilay;
	  ilay=NMAXLAYER; 
	}
	else
	  vpr1[ilay]=NODATAVPR;
	/*        } */
      }
     
    }
    if (vert_ext<VEXTMIN_VPR ){
      fprintf(log_vpr,"estensione profilo verticale troppo bassa\n");
      *ct=0;
      ilast=0;
      for  (il=0; il<NMAXLAYER; il++) vpr1[il]=NODATAVPR;
      return(1);
    }
   
  }
  else {
    fprintf(log_vpr,"volume precipitante troppo piccolo\n");
    *ct=0;
    ilast=0;
    for  (il=0; il<NMAXLAYER; il++) vpr1[il]=NODATAVPR;
    return(1);
  }  
 
 
  for (ilay=ilast; ilay<NMAXLAYER; ilay++) vpr1[ilay]= vpr1[ilay-1]-0.03*(float)TCK_VPR/1000.; /*USO GRADIENTE STANDARD PER PARTE ALTA DEL PROFILO..*/

  /*  calcolo la varianza del volume di R a 1100 m normalizzata rispetto al volume totale precipitante */


  if (area_vpr[5] > 0){
    media=vpr1[5]*(float)area_vpr[5];
    for (k=0; k<icounter; k++){
      somma=somma+pow(((float)cappi1100[k]-media),2);
    }
    printf("somma %f\n",somma);
    stdev=(sqrt(somma)/(icounter+1))/((float)area_vpr[5]*vpr1[5]);
  }
 
  fprintf(log_vpr,"stdevN= %f  \n",stdev);

  if (stdev>THR_STDEVR){ 
    for  (il=0; il<NMAXLAYER; il++) vpr1[il]=NODATAVPR; 
    fprintf(log_vpr,"stdevN maggiore di soglia");
    return (1);
  }
  else
    return(0);
}

/*
  comstart comp_levels
  idx esegue composizione dei profili istantaneo e 'ultimo'
  esegue composizione dei profili istantaneo e 'ultimo'
  se il profilo istantaneo è più 'corto', inserisco i valori dell'ultimo profilo negli strati'mancanti'  

  float v0  ultimo profilo
  float v1  profilo istantaneo
  float nodata  
  float peso     

  result=((1.-peso)*v0+peso*v1)  
  comend
*/
/*===============================================*/
 
float comp_levels(v0, v1, nodata, peso) 
 
     float v0,v1,nodata,peso;
{
  float result;
  if ((v0<nodata+1)&&(v1<nodata+1)) result=nodata;
  if (v0<nodata+1) result=v1;
  if (v1<nodata+1) result=v0;                            
  if ((v0>nodata+1) && (v1>nodata+1)  ) result=((1.-peso)*v0+peso*v1); /* in questa configurazione il vpr è di altezza costante  nel tempo ma un po' 'sconnesso' in alto*/
  return(result);
}
/*===============================================*/
/*
  comstart profile_gap
  idx calcola il no di quarti d'ora che intercorrono dall'ultimo profilo calcolato (combinato)
  calcola il no di quarti d'ora che intercorrono dall'ultimo profilo calcolato(combinato) memorizzato in 'LAST_VPR'         
  comend
*/
long int profile_gap(char nomefile[])
{ 
  FILE *file;
  time_t last_time;
  long int gap1;
  	

  file = fopen(nomefile,"r");
  if(file != NULL ){ /*contemplo la prima iterazione dopo installazione*/
    //file = controllo_apertura(nomefile," ultimo vpr ","r");
    fread(&last_time,4,1,file);
    fclose(file); 
    gap1=abs(old_data_header.norm.maq.acq_date-last_time)/900;	
    fprintf (log_vpr,"old_data_header.norm.maq.acq_date last_time gap %d %ld %ld \n",old_data_header.norm.maq.acq_date,last_time,gap1);
  }
  else{
    gap1=100 ;
  }
  return(gap1);
}

/*===============================================*/
/*
  comstart combina_profili
  idx combina il profilo verticale corrente con quello precedente
  combina il profilo verticale corrente con quello precedente tramite il metodo di Germann (2003)
  a) calcolo gap tra ultimo profilo e istante corrente
  b) se gap < MEMORY leggo profilo 
  c) se gap> MEMORY peso 0 il profilo storico e cerco oltre la data (per casi vecchi)
  d) faccio func_vpr
  f) cerco il profilo con cui combinare (->proprio, se gap<MEMORY ->dell'altro radar se gap_res<MEMORY e profile_heating_res=WARM)
  g) Combino livelli con peso sottostante
  Dati cv e ct, volume totale e volume precipitante il peso del vpr istantaneo è calcolato come segue:
  c0=2*(*cv);
  peso=(float)(*ct)/(c0+(*ct))
  long int c0,*cv,*ct; costanti di combinazione (v. ref.)
  h) trovo livello minimo, se livello minimo profilo combinato più alto del precedente calcolo la diff media e sommo al vecchio
  e) ricalcolo livello minimo 
  float vpr0[NMAXLAYER],vpr1[NMAXLAYER],vpr[NMAXLAYER]; profilo precedente, ultimo e combinato
  float alfat,noval; peso, nodata
  FILE *file;   
  int mode,ilay;  modalità calcolo profilo (0=combinazione, 1=istantaneo),indice di strato    
  comend
*/
int combina_profili(char *sito, char *sito_ad)
 
{ 
  long int c0,*cv,*ct;
  float vpr0[NMAXLAYER],vpr1[NMAXLAYER],vprmax,vpr_dbz;
  float alfat,noval;
  FILE *file;
  int mode,ilay,gap_res,heating_res,i,foundlivmin=0,il,ier_ap;
  long int  area[NMAXLAYER],ar;
  int n=0,diff=0;
  char nomefile[150],stringa[100];
  struct tm *T_tempo;
  time_t Time,T_Time;
  int ier_cp;

  mode=MOD_VPR;

  for (i=0;i<NMAXLAYER;i++) {
    area[i]=-9999;//inizializzo area
    area_vpr[i]=0;
    vpr0[i]=NODATAVPR;
    vpr1[i]=NODATAVPR;  
  }
 
  
  noval=NODATAVPR;
  vprmax=0.;
  /* questo per fare ciclo sul vpr vecchio*/
  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
    
  /* READ OLD VPR */
  gap=profile_gap(getenv("LAST_VPR"));
  //CONTROLO CHE IL FILE CI SIA
  file=fopen(getenv("VPR0_FILE"),"r");
  if(file == NULL ) {
    fprintf(log_vpr,"non esiste file vpr vecchio \n",getenv("VPR0_FILE"));
    gap=100;    
  }
  if (gap<=MEMORY){     
    printf("%s\n",getenv("VPR0_FILE")); 
    controllo_apertura(getenv("VPR0_FILE")," old VPR ","r"); 
    for(ilay=0; ilay<NMAXLAYER; ilay++){
     ier_cp=fscanf(file,"%f %li\n",&vpr0[ilay],&area[ilay]);
    }
    fprintf(log_vpr,"fatta lettura vpr\n");
    fclose(file);
 
  }
  /*  SE GAP > MEMORY FACCIO UN CICLO DI RICERCA DEI PROFILI DOPO LA DATA ENTRO UN NUMERO DI QUARTI D'ORA PARI A MEMORY */
  else {
    for (i=0;i<MEMORY;i++){
      T_Time=Time+i*900;
      T_tempo=gmtime(&T_Time);

      sprintf(nomefile,"%s/%04d%02d%02d%02d%02d_vpr_%s",getenv("DIR_STORE_VPR"),
	      T_tempo->tm_year+1900, T_tempo->tm_mon+1, T_tempo->tm_mday,
	      T_tempo->tm_hour, T_tempo->tm_min,getenv("SITO"));

      //printf("provo ad aprire %s\n",nomefile);	
      ier_ap=access(nomefile,R_OK);  	
      if  (!ier_ap){
	file=fopen(nomefile,"r");
	gap=0;
	heating=WARM;  //profilo caldo 
	fscanf(file," %s" ,stringa); 
	for (ilay=0;  ilay<NMAXLAYER; ilay++){
	  ier_cp=fscanf(file," %i %f %li", &il, &vpr_dbz, &ar);  
	  if (vpr_dbz>0){ 
	    vpr0[ilay] = DBZtoR(vpr_dbz,aMP,bMP);
	    area[ilay]=ar;
	  }
	  else
	    vpr0[ilay] = NODATAVPR;	    
	}
	break;
      }	       
    }
  }
  /* da mettere qui */
      
  
  /* SECTION C profiles combination using recoursive formula and new vpr printing*/
  cv=( long int *)malloc(sizeof(long int));
  ct=( long int *)malloc(sizeof(long int));
  if (ct == NULL) printf("malloc fallita per ct\n");
  if (cv == NULL) printf("malloc fallita per cv\n");

  ier_vpr=func_vpr(cv,ct,vpr1,area_vpr,sito); // ho fatto func_vpr, il profilo istantaneo
  fprintf(log_vpr,"fatta func vpr\n");
  c0=2*(*cv);
  if(mode == 0) {    /*modalità VPR combinato*/
    fprintf(log_vpr,"gap %li \n",gap);
    if (gap>MEMORY) 
      { 	

        gap_res=profile_gap(getenv("LAST_VPR_RES"));
	fprintf(log_vpr,"gap_res %i  \n",gap_res );
        if (gap_res<=MEMORY){
	  file=controllo_apertura(getenv("VPR_HEATING_RES")," riscaldamento vpr ","r");
	  fscanf(file,"%i ",&heating_res);
	  fprintf(log_vpr," heating_res %i \n",heating_res);
	  fclose(file); 
	  if (heating_res >= WARM){
            file=controllo_apertura(getenv("VPR0_FILE_RES")," ultimo vpr secondo radar","r");
            for(ilay=0; ilay<NMAXLAYER; ilay++){
	      fscanf(file,"%f %li \n",&vpr0[ilay],&area[ilay]);
	      gap=gap_res;
 	    }
	  }
	}
 	else  /////////qui va ristrutturato, cioè se non ho il vecchio ma ho il nuovo ..organizzare 
	  if (ier_vpr){
	    free(cv);
	    free(ct);
	    return (1);   
	  }
      }

    if (!ier_vpr) {  // se il calcolo dell'istantaneo è andato bene  (////va aggiunto il caso in cui il vecchio non c'è)
      alfat=(float)(*ct)/(c0+(*ct)); // peso vpr corrente per combinazione
      for (ilay=0;  ilay<NMAXLAYER; ilay++){
	vpr[ilay]=comp_levels(vpr0[ilay],vpr1[ilay],noval,alfat);// combino livelli 
	if (area_vpr[ilay] > MIN_AREA && !foundlivmin && ilay>=0) { //trovo livello minimo 
	  livmin=ilay*TCK_VPR+TCK_VPR/2; /* primo livello VPR (sopra il livello 0) */
	  foundlivmin=1;
	}
      }
    }
    else { // se il calcolo dell'istantaneo non è andato bene , ricopio l'altro vpr, la sua area(? discutibile) e trovo il livello minimo 

      for (ilay=0;  ilay<NMAXLAYER; ilay++){
	area_vpr[ilay]=area[ilay];
	vpr[ilay]=vpr0[ilay];
	if (area_vpr[ilay] > MIN_AREA && !foundlivmin && ilay>=0) {
	  livmin=ilay*TCK_VPR+TCK_VPR/2; /* primo livello VPR (sopra il livello 0) */
	  foundlivmin=1;
	}
      }
    }
    // inserisco un pezzo di prova che calcola la differenza media tra 2 profili successivi e sìstema i pti in basso usandola 
    // in pratica li calcolo dal livello minimo del vpr1 
    //

    for (ilay=livmin/TCK_VPR;  ilay<NMAXLAYER; ilay++){
      if ( vpr0[ilay]> NODATAVPR && vpr[ilay]>NODATAVPR ){ 
	diff=diff + vpr[ilay]-vpr0[ilay];
	n=n+1;
      }
    }
    if (n>0){
      diff=diff/n;
      for (ilay=0; ilay<livmin/TCK_VPR; ilay++){ 
	if (vpr[ilay]> NODATAVPR) {
	  vpr[ilay]=vpr[ilay]+diff;
	  
	}
      }
    }
    //}  
    /// questo pezzo è indipendente dal precedente ma sembra sensato fare così , una volta combinati i profili prendere il livello minimo del combinato
    // resetto livello minimo e lo calcolo sul combinato..togliere.. 
    livmin=0;
    foundlivmin=0;
    for (ilay=0; ilay<NMAXLAYER; ilay++){
      if (vpr[ilay]> NODATAVPR && !foundlivmin) {
	livmin=ilay*TCK_VPR+TCK_VPR/2;
	foundlivmin=1;
      }
    }
    /// fine pezzo di prova
    fprintf(log_vpr," livmin %i  \n", livmin);
    if (livmin>=(NMAXLAYER-1)*TCK_VPR+TCK_VPR/2  || !foundlivmin) return (1);

  }   
  /*fine VPR combinato*/  
  else for (ilay=0; ilay<NMAXLAYER; ilay++) vpr[ilay]=vpr1[ilay];

    livmin=TCK_VPR;
    foundlivmin=0;
    for (ilay=livmin/TCK_VPR; ilay<NMAXLAYER; ilay++){
      if (vpr[ilay]> NODATAVPR && !foundlivmin) {
	livmin=ilay*TCK_VPR+TCK_VPR/2;
	foundlivmin=1;
      }
    }
 
  free(cv);
  free(ct);

  file=controllo_apertura(getenv("VPR0_FILE")," ultimo vpr ","w");
  for (ilay=0;  ilay<NMAXLAYER; ilay++)
    fprintf(file," %10.3f %li\n",vpr[ilay], area_vpr[ilay]); 
  fclose(file);

   
  return(0);
}

/*===============================================*/
/*
  comstart profile_heating
  idx calcola quanto il profilo è 'caldo'
  calcola quanto il profilo è 'caldo' 
  heating=heating-gap; se il profilo non è stato aggiornato
  heating=heating-gap+2; se il profilo è stato aggiornato
  heating=MEMORY;   se heating raggiunge WARM resta costante finchè non inizia raffreddamento
  comend
*/

int profile_heating()
#include <vpr_par.h>

{
  int heating;
  FILE *file;
 
  // read  VPR HEATING

  //file=controllo_apertura(getenv("VPR_HEATING")," riscaldamento vpr ","r");

  if (!access(getenv("VPR_HEATING"),F_OK)){
    file=fopen(getenv("VPR_HEATING"),"r");
    if (file == NULL )  fprintf (log_vpr,"non ho il file di riscaldamento %s \n",getenv("VPR_HEATING"));
    fscanf(file,"%i ",&heating);
    fclose(file); 
  } /*contemplo la prima iterazione dopo installazione testando l'esistenza del file*/
  else heating=0;
  if (ier_vpr){
    heating=heating-gap; /*se il profilo non è stato aggiornato, ho raffreddamento, in caso arrivi sotto WARM riparto da 0, cioè serve riscaldamento  */     
  }
  else  { 
    heating=heating-gap+2; /*se il profilo è stato aggiornato, ho riscaldamento , in caso arrivi sopra WARM riparto da MEMORY  */  
    if (heating>=WARM) heating=MEMORY;  /* se heating raggiunge WARM allora lo pongo uguale a MEMORY     */
  }  
  if (heating<0) heating=0;   
  file=controllo_apertura(getenv("VPR_HEATING")," riscaldamento vpr ","w");  
  fprintf(file," %i \n",heating); 
  fclose(file);
  fprintf (log_vpr,"ier_vpr %i ier_comb %i \n",ier_vpr,ier_comb);
  if (!ier_vpr){
    file = controllo_apertura(getenv("LAST_VPR")," ultimo vpr ","w");
    fwrite(&old_data_header.norm.maq.acq_date,4,1,file);
    fclose(file); 
  }
  fprintf (log_vpr,"gap %li heating %i \n",gap,heating);
  return(heating);
}

/*===============================================*/
/*
  comstart stampa_vpr
  idx stampovpr in dbz con la data
  comend
*/
int stampa_vpr()
{
  float vpr_dbz;
  int ilay;
  FILE *file;

  file=controllo_apertura(getenv("VPR_ARCH")," ultimo vpr in dBZ per il plot","w");
    fprintf(file," QUOTA   DBZ    AREA PRECI(KM^2/1000)\n" ); 
  for (ilay=0;  ilay<NMAXLAYER; ilay++){
    if (vpr[ilay]>NODATAVPR) {
      vpr_dbz=RtoDBZ(vpr[ilay],aMP,bMP);
      fprintf(file," %i %10.3f %li\n", ilay*TCK_VPR+TCK_VPR/2, vpr_dbz, area_vpr[ilay]);
    }
    else
      fprintf(file," %i %10.3f %li\n", ilay*TCK_VPR+TCK_VPR/2, NODATAVPR, area_vpr[ilay]);
  }
  fclose(file);
  return 0;
}
/*
  =======================================================================================*/
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

int corr_vpr(char *sito)
#include <vpr_par.h>

{
 
  int ilray,ilref,ilay2,i,k,ier,snow,strat;
  float corr,vpr_liq,vpr_hray,hbin,hliq,elevaz;
  snow=0;
  vpr_liq=NODATAVPR;
  memset(neve,0,sizeof(neve));
  hliq=NODATAVPR;
  ier=analyse_VPR(&vpr_liq,&snow,&hliq,sito);

  fprintf (log_vpr,"ier_analisi %i \n",ier) ;
  if (ier) return 1; /* se analisi dice che non è il caso di correggere non correggo (NB in questo caso non riempio la matrice di neve)*/
  fprintf (log_vpr,"CORREGGO VPR \n") ;
  fprintf (log_vpr,"altezza bright band %i \n",hvprmax) ;
  for (i=0; i<NUM_AZ_X_PPI; i++){
    for (k=0; k<vol_pol[0][i].b_header.max_bin; k++){      
      corr=0.;
      //hbin=(float)hray[k][elev_fin[i][k]];
      elevaz=(float)(vol_pol[elev_fin[i][k]][i].teta_true)*CONV_RAD;      
      hbin=quota_f(elevaz,k);
      if (snow) neve[i][k]=1;
      strat=1;
#ifdef CLASS	
      if (conv[i][k] >= CONV_VAL){ 	
	strat=0;
      }
#endif     
      if (BYTEtoDB(vol_pol[0][i].ray[k])>THR_CORR && hbin > hliq  && strat){   
	ilray=floor((hbin)/TCK_VPR);
	if (ilray>= NMAXLAYER ) ilray=NMAXLAYER-2;
        ilref=(dem[i][k]>livmin)?(floor(dem[i][k]/TCK_VPR)):(floor(livmin/TCK_VPR));
	
	if ((int)hbin%TCK_VPR > TCK_VPR/2) ilay2=ilray+1;
	else ilay2=(int)(fabs(ilray-1));
	if (vpr[ilref]>0 && vpr[ilray]>0 ){ /*devo avere un dato DI PIOGGIA nel VPR!*/
	  vpr_hray=vpr[ilray]+((vpr[ilray]-vpr[ilay2])/(ilray*TCK_VPR-TCK_VPR/2-ilay2*TCK_VPR))*(hbin-ilray*TCK_VPR-TCK_VPR/2);	/*per rendere la correzione continua non a gradini */ 		 
	  if (dem[i][k]> hvprmax+HALF_BB-TCK_VPR || snow){ /*classifico neve*/            
	    neve[i][k]=1;
	  }
	  if(snow)
	    corr=RtoDBZ(vpr[ilref],aMP,bMP)-RtoDBZ(vpr_hray,aMP,bMP);
          else
	    corr=RtoDBZ(vpr_liq,aMP,bMP)-RtoDBZ(vpr_hray,aMP,bMP);/*riporto comunque al valore liquido anche se sono sopra la bright band*/
	  if (corr>MAX_CORR) corr=MAX_CORR; /*soglia sulla massima correzione*/
	  if (hbin<hvprmax && corr>0.) corr=0; /*evito effetti incrementi non giustificati*/
	  vol_pol[0][i].ray[k]=DBtoBYTE(BYTEtoDB(vol_pol[0][i].ray[k])+corr);
          if (vol_pol[0][i].ray[k] < 1 ) vol_pol[0][i].ray[k]=1; /*evito di creare dei no data*/
	  corr_polar[i][k]=(unsigned char)(corr)+128;
	  //printf(" i %i k %i hbin %f hliq %f corr %f \n",i,k, hbin,hliq,corr);

	}
      }
    }
  } 
  return(0);
}


int trovo_hvprmax(int *hvprmax)

{
  int i,imax,istart,foundlivmax;
  float vprmax,h0start,peak;

  istart=livmin/TCK_VPR+1;
  if (! ier_t ){
    if ( t_ground > T_MAX_ML){
      printf("trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca\n");
      h0start=t_ground/9.8*1000 ;
      istart=h0start/TCK_VPR -2;
      printf("istart %i  \n",istart);
    }
  }
  else {
    printf("trovo hvprmax  a partire da livmin\n");
  }

  /* trovo hvprmax e il suo livello a partire dal livello istart*/
  foundlivmax=0;
  hvprmax=INODATA; 
  vprmax=NODATAVPR;
  imax=istart;	    
  for (i=istart;i<NMAXLAYER;i++) //la ricerca è un po' diversa dall'originale.. trovo il picco + alto con valore  rispetto a 4 sopra > soglia 
    { 
      peak=10*log10(vpr[i]/vpr[i+4]);
      if (vpr[i]>vpr[i-1]  && peak> MIN_PEAK_VPR ) // se vpr(i) maggiore del massimo e picco sufficientemente alto
	{
	  imax=i;
	  vprmax=vpr[imax];
	}
    }
  if ( imax  > istart ){
    foundlivmax=1;
    peak=10*log10(vpr[imax]/vpr[imax+4]);
    hvprmax=imax*TCK_VPR+TCK_VPR/2;
    printf("trovato ilaymax %i %i \n",hvprmax,imax);
    printf (" picco in dbR %f  \n",peak);
  }

  printf("exit status trovo_hvprmax %i \n",foundlivmax);
  return (foundlivmax);
} 
   

/*
  comstart analyse_VPR
  idx analisi del profilo calcolato
  analisi del profilo calcolato: 

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


int analyse_VPR(float *vpr_liq,int *snow,float *hliq, char *sito)
    
{
  int i,ier=1,ier_ana=0,ier_ap,liv0,ier_max;
  int nstaz,nvar,tipo_profilo;
  int npar=5;
  float *a,*dyda,v1000,v1500,v600sottobb,vliq,vhliquid,vprmax; //*togliere gli ultimi tre*/;
  char date[20];
  struct tm *tempo;
  time_t Time;
  int ndata=15;
  FILE *file;
  // inizializzazioni

  nstaz=4;
  nvar=1;

  strcpy(date,"000000000000");
  hvprmax=INODATA;// inizializzo
  tipo_profilo=-1;
  stdev=NODATAVPR;
  chisqfin=100.;
  vliq=NODATAVPR;
  vhliquid=NODATAVPR;
  v600sottobb=NODATAVPR;
  v1000=NODATAVPR;
  v1500=NODATAVPR;
  vprmax=NODATAVPR;
  a=vector(1,npar);
  dyda=vector(1,npar); 
  for (i=1;i<=npar;i++){
    a[i]=NODATAVPR;
    dyda[i]=NODATAVPR;
  }
 
  ier_max=trovo_hvprmax(&hvprmax);
  
  if (ier_t) //1  se non ho nè T nè il massimo esco altrimenti tipo_profilo=0 
    {
      fprintf(log_vpr,"non ho T,...  \n");

      if ( ! ier_max ) {
	fprintf(log_vpr," non ho trovato hvprmax, nè T , esco \n");
	free_vector(a,1,npar);
	free_vector(dyda,1,npar);
	return 1;
      } 
      tipo_profilo=0;
    }
  else 
    {
     
      if (t_ground >= T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100.){ //1  se T > T_MAX_ML e non ho il massimo esco
	if ( ! ier_max ) {
	  fprintf(log_vpr," non ho trovato hvprmax, esco \n");
	  free_vector(a,1,npar);
	  free_vector(dyda,1,npar);
	  return 1;
	}
	tipo_profilo=0; 
      }

     
      if (t_ground >= T_MIN_ML  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100.)
      	{

	  if ( ! ier_max ) {
	    if (t_ground <= T_MAX_SN ){
		tipo_profilo=3;
	      }
	      else{
		fprintf(log_vpr," non ho trovato hvprmax, esco \n");
		free_vector(a,1,npar);
		free_vector(dyda,1,npar);
		return 1;
	      }
	    
	  }

	  liv0=livmin+HALF_BB;  
      	  if (hvprmax > liv0) fprintf(log_vpr," il livello %i è sotto la Bright band, ma T bassa  interpolo\n",livmin);
      	  else fprintf(log_vpr," il livello %i potrebbe essere dentro la Bright Band, interpolo\n",livmin);
      	  tipo_profilo=2;
      	}

      if (t_ground < T_MAX_SN )  {
	if ( ier_max ){ 		
	  fprintf(log_vpr," temperatura da neve e massimo in quota\n");
	  tipo_profilo=2;	 
	}

	//    }
      	else {
      	  fprintf(log_vpr," temperatura da neve e massimo al suolo, non interpolo\n");
      	  tipo_profilo=3;
        }
      }
    }
  

  switch
    (tipo_profilo)
    {
    case 0:
    case 1:
    case 2:
  
      ier=interpola_VPR(a,npar);
      if (ier){ 
	fprintf(log_vpr," interpolazione fallita \n");           
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
	    //*vpr_liq=vpr[(hvprmax+HALF_BB-100)/TCK_VPR];/*21 aprile 2008*/
	    *vpr_liq=vpr[(hvprmax+1000)/TCK_VPR]*2.15;/*21 aprile 2008*/
	    *hliq=0;		         
	    break;
	  }  
      }
      else{
  	fprintf(log_vpr," interpolazione eseguita con successo \n");    	
	/*calcolo valore di riferimento di vpr_liq per l'acqua liquida nell'ipotesi che a[2]=quota_bright_band e a[2]-1.5*a[3]=quota acqua liquida*/
        if (tipo_profilo == 2 ) { 	 
	  *hliq=(a[2]-2.1*a[3])*1000.;
	  lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);	 
         
	  *hliq=0;  /*con casi di bright band bassa.. cerco di correggere il più possibile*/
	  //*vpr_liq=vpr[(hvprmax+(int)(7.*a[3]*1000))/TCK_VPR]*2.15; /*23 aprile 2008 cioè vpr neve + 6 dB*/
	  *vpr_liq=vpr[(hvprmax+1000)/TCK_VPR]*2.15;
	}
        else { 	  
	  *hliq=(a[2]-2.1*a[3])*1000.;

	  lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);
	  if ( *hliq > livmin) {
	    *vpr_liq=vpr[(int)(*hliq/TCK_VPR)]; // ... SE HO IL VALORE VPR USO QUELLO.
	  }
	  else // altrimenti tengo il valore vpr neve + 6 dB* e metto tipo_profilo=2
	    {
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
  fprintf(log_vpr,"TIPO_PROFILO= %i vpr_liq %f hliq %f \n", tipo_profilo, *vpr_liq,*hliq );
 
  /* nome data */
  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&Time);
  sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
	  tempo->tm_mday,tempo->tm_hour, tempo->tm_min);

  /* parte di stampa test vpr*/

 
  if (! ier ) {
    if(*hliq > livmin +200 )
      vhliquid=RtoDBZ(vpr[(int)(*hliq)/TCK_VPR],aMP,bMP);
    vliq=RtoDBZ(*vpr_liq,aMP,bMP);
  }
  if (ier_max) {
    if ( hvprmax-600 >= livmin ) 
      v600sottobb=RtoDBZ(vpr[(hvprmax-600)/TCK_VPR],aMP,bMP);
    if ((hvprmax+1000)/TCK_VPR < NMAXLAYER )
      v1000=RtoDBZ(vpr[(hvprmax+1000)/TCK_VPR],aMP,bMP);
    if ((hvprmax+1500)/TCK_VPR < NMAXLAYER )
      v1500=RtoDBZ(vpr[(hvprmax+1500)/TCK_VPR],aMP,bMP);
    vprmax=RtoDBZ(vpr[(hvprmax/TCK_VPR)],aMP,bMP);
  }
  //fprintf(test_vpr,"%s %i %i %f %f %f  %f %f %f %f %f %f %f %f \n",date,hvprmax,tipo_profilo,stdev,chisqfin,*hliq,vliq,vhliquid,v600sottobb,v1000+6,v1500+6,v1000,v1500,vprmax);
 

  fprintf(test_vpr,"%s %i %i %f %f %f  %f %f %f %f %f %f %f %f \n",date,hvprmax,tipo_profilo,stdev,chisqfin,*hliq,vliq,vhliquid,v600sottobb,v1000+6,v1500+6,v1000,v1500,vprmax);	     
  free_vector(a,1,npar);
  free_vector(dyda,1,npar);
  /* fine parte di stampa test vpr*/
  fclose(test_vpr);
  /* fine parte di stampa test vpr*/

  ier_ap=access(getenv("VPR_HMAX"),R_OK);   
  file = fopen(getenv("VPR_HMAX"),"w");
  fprintf(file,"%d", hvprmax);
  fprintf(log_vpr,"fatta scrittura hmax vpr = %d \n",hvprmax);
  fclose(file);
  return (ier_ana);
}

int get_t_ground(float *t_gr)


# include <geo_par.h>
{
  int icount,ierr,ier_file_t;
  float t,media_t,lon,lat,radar_lat,radar_lon;
  FILE *file_t;
  char *sito;
  /*apertura file T*/
  ierr=1;     
  media_t=NODATAVPR-1;   
  ier_file_t=access(getenv("FILE_T"),R_OK);
  if (!ier_file_t) { 
    file_t=fopen(getenv("FILE_T"),"r");       
    media_t=0;
    icount=0;
    while (1) {
      if(fscanf(file_t,"%f %f %f \n",&lon,&lat,&t) == EOF) break;
      sito=getenv("SITO");
      if (!(strcmp(sito,"SPC"))) {
	radar_lat=SPC_LAT;
	radar_lon=SPC_LON; }
      if (!(strcmp(sito,"GAT"))) {
	radar_lat=GAT_LAT;
	radar_lon=GAT_LON; }
      if (fabs(radar_lat-lat)<=maxdlat && fabs(radar_lon-lon)<=maxdlon) {
	icount+=1;
	media_t+=t-273.15;
      }
    }
    fclose(file_t);     
  }
  else {
    fprintf(log_vpr,"non trovo file delle temperature al suolo vicino al radar %s \n",getenv("FILE_T"));
    return(ierr);
  }
  if (icount>0)
    {
      media_t/=(float)(icount);
      fprintf(log_vpr,"ho %i stazioni dati affidabili e la t media è  %f\n",icount,media_t);
      *t_gr=media_t;
      ierr=0;
    }
  else {
    ierr=1;
    fprintf(log_vpr,"non trovo dati di temperatura\n");
  }
   
  return(ierr);
}

/*


  comstart interpola_VPR
  idx interpola il profilo verticale tramite una funzione lingauss 
  interpola il profilo verticale tramite una funzione gaussiana + lineare del tipo 
    
  y= B*exp(-((x-E)/G)^2)+C+Fx

  usa la funzione mrqmin delle numerical recipes in C: tutti i vettori passati a mrqmin devono essere allocati e deallcocati usando le funzioni di NR (vector, matrix, free_vector, free_matrix.. etc) che definiscono vettori con indice a partire da 1.
  NB gli ndata dati considerati partono da 1000 m sotto il massimo (in caso il massimo sia più basso di 1000 m partono da 0 m)
  A ogni iterazione si esegue un test sui parametri. Se ritorna 1 si torna ai valori dell'iterazione precedente.
  A fine interpolazione si verifica che il chisquare non superi una soglia prefissata, in tal caso ritorna 1 e interpol. fallisce.

  INIZIALIZZAZIONE PARAMETRI:  
  a[1]=B=vpr(liv del massimo)-vpr(liv. del massimo+500m);
  a[2]=E= quota liv. massimo vpr (in KM);
  a[3]=G=semiampiezza BB (quota liv .massimo- quota massimo decremento vpr nei 600 m sopra il massimo ) in KM;
  a[4]=C=vpr(liv massimo + 700 m);
  a[5]=F=coeff. angolare segmento con estremi nel vpr ai livelli max+900m e max+1700m, se negativo =0.;

 
  float a[ma], int ma: vettore parametri e n0 parametri
  float *x, *y:  quote in KM e valori del vpr usati per l'interpolazione
  float *sig,alamda : vettore dev st. e variabile che decrementa al convergere delle iterazioni
  float *dyda: vettore derivate rispetto ai parametri 
  float B,E,C,G,F:  parametri da ottimizzare, first guess
  float chisq; scarto quadratico
  int i,in1,in2,in3,in4,*ia,ifit,ii,ndati_ok,k;
  int ndata=15;  numero di dati considerati
  float **covar,**alpha; matrice ovarianze, matrice alpha

  comend
*/
int interpola_VPR(float a[], int ma)

# include <vpr_par.h>

{
  float *x, *y,*sig,alamda,y1,*dyda,B,E,C,G,F,xint,qdist,*abest;
  float chisq=100.;
  float chisqold=0.0;
  float chisqin=0.0;
  int i,in1,in2,in3,in4,*ia,ifit,ii,ndati_ok,k,ier_int;
  //int ma=5;
  int ndata=12;  
  float **covar;
  float **alpha;
  char file_vprint[200];
  FILE *file;

  fprintf(log_vpr,"sono in interpola_vpr\n");
  ier_int=0;

  in2=(hvprmax+HALF_BB)/TCK_VPR; //indice del massimo + 500 m
  in1=(hvprmax-TCK_VPR/2)/TCK_VPR; //indice del massimo
  in3=in2+1;
  in4=in2+5; //indice del massimo + 1000 m
  fprintf(log_vpr,"in1 in2 %i %i %f %f \n",in1,in2,vpr[in1],vpr[in2]);

  if (in4 > NMAXLAYER-1) {
    ier_int=1;
    return ier_int;
  }

  /* inizializzazione vettore parametri */ 
  abest=vector(1,ma);
  x=vector(1,ndata);
  y=vector(1,ndata);
  sig=vector(1,ndata);
  ia=ivector(1,ma);
  covar=matrix(1,ma,1,ma);
  alpha=matrix(1,ma,1,ma);

  for (k=in1+2; k<=in3; k++)  
    { 
      ier_int=0;
      dyda=vector(1,ma);
    
      a[1]=B=vpr[in1]-vpr[in2];
      a[2]=E=hvprmax/1000.;
      a[3]=G=(k-in1-0.5)*TCK_VPR/1000.;
      a[4]=C=vpr[in2];
      a[5]=F=vpr[in4]<vpr[in3]?(vpr[in4]-vpr[in3])/((in4-in3)*TCK_VPR/1000.):0.;

      /*       dermax=0; */
      /*       der=fabs(vpr[k-1]-vpr[k]); */
      /*       if ( der > dermax) {  */
      /*       a[3]=G=(k-in1-0.5)*TCK_VPR/1000.;// almeno 0.3 km fino al max 0.5 */
      /*        dermax=der; */
      /*        } */
      /*         }  */
      
      /*        a[4]=C=vpr[in2]; */
      /*       a[5]=F=vpr[in4]<vpr[in3]?(vpr[in4]-vpr[in3])/((in4-in3)*TCK_VPR/1000.):0.; */
      fprintf(log_vpr,"B F %f %f\n",B,F);
      /* fine inizializzazione vettore parametri */ 
      
      alamda=-0.01;

      for (i=1;i<=ma;i++) ia[i]=1;
      y1=0;     
      ii=1;
      ndati_ok=ndata;
      qdist=0;
      for (i=1; i<=ndati_ok; i++)
	{ 
	  sig[ii]=0.5; 
	  //x[ii]= ((hvprmax-1000.)>livmin)? (i*TCK_VPR+(hvprmax-800)-TCK_VPR)/1000. : (livmin+(i-1)*TCK_VPR)/1000.;
	  //y[ii]= ((hvprmax-1000.)>livmin)? vpr[i+((hvprmax-800)-TCK_VPR)/TCK_VPR] : vpr[i-1+livmin/TCK_VPR]; 
	  x[ii]= ((hvprmax-800.)>livmin)? (i*TCK_VPR+(hvprmax-600)-TCK_VPR)/1000. : (livmin+(i-1)*TCK_VPR)/1000.;
	  y[ii]= ((hvprmax-800.)>livmin)? vpr[i+((hvprmax-600)-TCK_VPR)/TCK_VPR] : vpr[i-1+livmin/TCK_VPR]; 	  
	  lineargauss(x[ii], a, &y1, dyda, ndata);
	  qdist=(y1-y[ii])*(y1-y[ii]);
	  if (sqrt(qdist) < DIST_MAX) 
	    { 
	      ii+=1;
	      chisqin=qdist+chisqin;
	    }
	  else 
	    ndati_ok=ndati_ok+1;  
	  if    ( ndati_ok > 40 || ndati_ok - ndata > 2 ) 
	    {
	      fprintf (log_vpr,"  esco dal ciclo di fit , interpolazione fallisce \n");  
	      ier_int=1;
	      break;
	    }
	} 
      
      if (!ier_int)
	{ 
	  fprintf (log_vpr,"\n alamda %f   chisqin % f a[1]  % f a[2] % f a[3]  % f a[4]  % f a[5]  % f\n", alamda,chisqin,a[1], a[2], a[3],a[4],a[5]);   
	  ifit=0;
	  while (fabs(chisq-chisqold) > DCHISQTHR && ifit < 1)
	    {
	      chisqold=chisq;
	      B=a[1];
	      E=a[2];
	      G=a[3];
	      C=a[4];
	      F=a[5];
	      mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, &chisq, &lineargauss, &alamda);
	      fprintf (log_vpr,"alamda %f   chisq % f a[1]  % f a[2] % f a[3]  % f a[4]  % f a[5]  % f\n", alamda,chisq,a[1], a[2], a[3],a[4],a[5]);      
	      ifit=testfit(a,chisq,chisqin);  /*test sul risultato del fit */
	      if (ifit)
		{/*test sul risultato del fit */
		  a[1]=B; /*test sul risultato del fit */
		  a[2]=E; /*test sul risultato del fit */
		  a[3]=G; /*test sul risultato del fit */
		  a[4]=C; /*test sul risultato del fit */
		  a[5]=F; /*test sul risultato del fit */
		  chisq=chisqin;
		}
	    } 
      
	  /*ultimo test sul risultato del fit */
	  /* 	if (chisq>CHISQ_MAX){ */
	  /* 	  ier_int=1; */
	  /* 	} */
	  /* 	else {  */
	  /* 	  sprintf(file_vprint,"%s_int",getenv("VPR0_FILE")); */
	  /* 	  file=controllo_apertura(file_vprint," vpr interpolato ","w");  */
	  /* 	  for (i=1; i<=ndata+20; i++) */
	  /* 	    { */
	  /* 	      xint=(i*TCK_VPR-TCK_VPR/2)/1000.;/ */
	  /* 	      lineargauss(xint, a, &y1, dyda, ndata);/  */
	  /* 	      fprintf(file," %10.3f \n",y1);  */
	  /* 	    } */
	  /* 	  fclose(file); */
	  if (chisq <chisqfin) 
	    { 
	      chisqfin=chisq;
	      for (i=1;i<=ma;i++) abest[i]=a[i];
	    }
	}
    }
  if (chisqfin>CHISQ_MAX)
    {
      ier_int=1;
    }
  else { 
    for (i=1;i<=ma;i++) a[i]=abest[i];
    //sprintf(file_vprint,"%s_int",getenv("VPR0_FILE"));
    sprintf(file_vprint,"%s_int",getenv("VPR_ARCH"));
    file=controllo_apertura(file_vprint," vpr interpolato ","w"); 
    for (i=1; i<=NMAXLAYER; i++)
      {
	xint=(i*TCK_VPR-TCK_VPR/2)/1000.;
	lineargauss(xint, a, &y1, dyda, ndata); 
	fprintf(file," %f \n",RtoDBZ(y1,aMP,bMP)); 
      }
    fclose(file);
  } 
  free_vector(dyda,1,ma);  
  free_vector(abest,1,ma);
  free_ivector(ia,1,ma); 
  free_vector(x,1,ndata);
  free_vector(y,1,ndata);
  free_vector(sig,1,ndata);
  free_matrix(alpha,1,ndata,1,ndata);
  free_matrix(covar,1,ndata,1,ndata);
  
  return ier_int; 
} 

int testfit(float a[], float chisq, float chisqin)
{
  if (a[1]<0. || a[1] >15.) return 1;
  if (a[2] >10.) return 1;
  if (a[3]<0. ) return 1;
  if (a[4]<0. ) return 1;
  if (a[5]>0 ) return 1;
  if (chisq>chisqin ) return 1;
  return 0;
}

/* comstart lineargauss

   calcola derivate della gaussiana linare rispetto ai parametri e valore (*y) in un punto x

   y(x,a) is the sum of a gaussian and a linear function with amplitude B=a[1], center E=a[2] and width G=a[3] and a linear function with coefficients C(shift)=a[4] and F(slope)=a[5] 
   a is the parameters vector and dyda is the vector of derivatives respect to the different parameters
   comend
*/
void lineargauss(float x, float a[], float *y, float dyda[],int na)
 
{
  float fac, ex, arg;

  *y=0.0;
  arg=(x-a[2])/a[3];
  ex=exp(-arg*arg);
  fac=a[1]*ex*2.0*arg;
  *y+=a[1]*ex+a[4]+a[5]*x;
  dyda[1]=ex;
  dyda[2]=fac/a[3];
  dyda[3]=fac*arg/a[3];
  dyda[4]=1.;
  dyda[5]=x;     
}

#endif

float quota_f(float elevaz, int k) // quota funzione di elev(radianti) e range
{
float dist;
float quota_true;
  dist=k*(int)(size_cell[old_data_header.norm.maq.resolution])+(int)(size_cell[old_data_header.norm.maq.resolution])/2.;	
  quota_true=(sqrt(pow(dist/1000.,2)+(rst*rst)+2.0*dist/1000.*rst*sin(elevaz))-rst)*1000.; /*quota in prop. standard  */;

  return quota_true;
}

#ifdef CLASS
void classifica_rain()
{
  float a;
  float range[MAX_BIN];
  float zz[MAX_BIN][NEL];
  float xx[MAX_BIN][NEL];
  int  i_xx[MAX_BIN][NEL],i_zz[MAX_BIN][NEL],i_xx_min[MAX_BIN][NEL],i_xx_max[MAX_BIN][NEL],i_zz_min[MAX_BIN][NEL],i_zz_max[MAX_BIN][NEL];
  int  im[MAX_BIN][NEL], ix[MAX_BIN][NEL], jm[MAX_BIN][NEL], jx[MAX_BIN][NEL];
  int i,j,kx,kz,k,iel,iaz,ibin,imin,imax,jmin,jmax,RHI_ind[NEL][MAX_BIN],jbb;
  int wimin,wjmin,wimax,wjmax;
  int hmax=-9999, ier_ap,ier_0term;
 
  //float w_size[2]={3.,1.5}; //dimensione della matrice pesi
  float w_size[2]={3.,0.3}; //dimensione della matrice pesi
  float **rhi_cart,**rhi_weight,RHI_beam[NEL][MAX_BIN],*w_x,*w_z,**w_tot,**beamXweight[MAX_BIN]; // da inizializzare in fase di programma
  float range_min,range_max,xmin,zmin,xmax,zmax;
  int w_x_size,w_z_size,w_x_size_2,w_z_size_2;
  FILE *file;

  resol[0]=RES_HOR;
  resol[1]=RES_VERT;
  a=REARTH;


  /* ;---------------------------------- */
  /* ;          FASE 0 :  leggo informazioni di temperatura da modello*/
  /* ;---------------------------------- */

  //file0=fopen("./src/0termico.PREV","r");
  ier_0term=trovo0term();

  /* READ OLD VPR_hmax */ 
  gap=profile_gap(getenv("LAST_VPR"));
  if (gap<=MEMORY){    
    ier_ap=access(getenv("VPR_HMAX"),R_OK);
    if (!ier_ap) {  
      file = fopen(getenv("VPR_HMAX"),"r");
      //file=controllo_apertura(getenv("VPR_HMAX")," altezza hmax ultimo vpr ","r");
      fscanf(file,"%i", &hmax);
      fprintf(log_vpr,"fatta lettura hmax vpr = %i \n",hmax);
      fclose(file);
    }
  }
   
  if (hmax <0 && !ier_0term) {
   
    htbb=zeroterm/1000. + 0.4; // se non ho trovato il vpr allora uso un range più ristetto, prob caso convettivo
    hbbb=zeroterm/1000. - 0.6;
   
  }
  else
    {   
      //THR_ZABB=1000.; // se ho trovato il vpr faccio crescere la soglia per limitare contaminazioni BBand
      //hbbb=(hmax-700.)/1000.;
      //htbb=(hmax+900.)/1000.;
      hbbb=(hmax-500.)/1000.;
      htbb=(hmax+700.)/1000.;

    }
  if (t_ground<T_MAX_ML){ 
    hbbb=0;
  }
  if (hmax < 2000 ) THR_ZABB=3000.;// piccola modifica per i casi invernali dove si rischia interferenza brigth band
 
  /* ---------------------------------- */
  /*           FASE 1 */
  /* ---------------------------------- */
  /*    Costruzione matrice generale per indicizzare il caricamento dei dati */
  /*    da coordinate radar a coordinate (X,Z) */
  /*  Calcolo distanza per ogni punto sul raggio */
  /*  xx contiene la distanza sulla superfice terrestre in funzione del range radar e dell'elevazione */
  /*  Metodo 1 - -Calcolo le coordinate di ogni punto del RHI mediante utilizzo di cicli */
 
  /* #ifdef TIME */
  /*  prendo_tempo(); */
  /* #endif */
 
  // estremi x e z (si procede per rhi)
  range_min=0.5*size_cell[old_data_header.norm.maq.resolution]/1000.;
  range_max=(MAX_BIN-0.5)*size_cell[old_data_header.norm.maq.resolution]/1000.;
  fprintf(log_vpr,"range_min = %f range_max = %f\n",range_min,range_max );
  xmin=floor(range_min*cos(elev_array[n_elev-1]*CONV_RAD)); // distanza del primo punto
  zmin=pow(pow(range_min,2.)+pow(4./3*a,2.)+2*range_min*4./3.*a*sin(elev_array[0]*CONV_RAD),.5) -4./3.*a+h_radar; // quota del primo punto
  xmax=floor(range_max*cos(elev_array[0]*CONV_RAD)); // distanza massima
  zmax=pow(pow(range_max,2.)+pow(4./3*a,2.)+2*range_max*4./3.*a*sin(elev_array[n_elev-1]*CONV_RAD),.5) -4./3.*a+h_radar;//quota massima
 
  x_size=(xmax-xmin)/resol[0]; //dimensione orizzontale
  z_size=(zmax-zmin)/resol[1]; //dimensione verticale
 
  w_x_size=ceil((w_size[0]/resol[0])/2)*2+1; //dimensione x matrice pesi
  w_z_size=ceil((w_size[1]/resol[1])/2)*2+1; //dimensione z matrice pesi
 
  if (w_x_size < 3)  w_x_size=3; 
  if (w_z_size < 3 ) w_z_size=3;
 
  w_x_size_2=w_x_size/2;
  w_z_size_2=w_z_size/2;
        
  w_x=(float *)malloc(w_x_size*sizeof(float));
  w_z=(float *)malloc(w_z_size*sizeof(float));
  w_tot=(float **) malloc(w_x_size*sizeof(float *));
  for(k=0;k<w_x_size;k++)
    w_tot[k]=(float *) malloc(w_z_size*sizeof(float));
    
  for (i=0; i<MAX_BIN; i++){
    range[i]=(i+0.5)*size_cell[old_data_header.norm.maq.resolution]/1000.;
      
    for (k=0; k<n_elev; k++){
      zz[i][k]=pow(pow(range[i],2.)+pow(4./3*a,2.)+2.*range[i]*4./3.*a*sin(elev_array[k]*CONV_RAD),.5) -4./3.*a+h_radar;// quota
      xx[i][k]=range[i]*cos(elev_array[k]*CONV_RAD); // distanza
      i_zz[i][k]=floor((zz[i][k]-zmin)/resol[1]);// indice in z, nella proiezione cilindrica, del punto i,k 
      i_xx[i][k]=floor((xx[i][k]-xmin)/resol[0]);// indice in x, nella proiezione cilindrica, del punto i,k 
      RHI_ind[k][i]=i_xx[i][k]+i_zz[i][k]*x_size;
      //shift orizzontale negativo del punto di indice i_xx[i][k] per costruire la finestra in x
      // se l'estremo minimo in x della finestra è negativo assegno come shift il massimo possibile e cioè la distanza del punto dall'origine
      i_xx_min[i][k]=i_xx[i][k];
      if (i_xx[i][k]-w_x_size_2 >= 0) 
	i_xx_min[i][k]= w_x_size_2;
	
      //shift orizzontale positivo attorno al punto di indice i_xx[i][k] per costruire la finestra in x
      i_xx_max[i][k]=x_size-i_xx[i][k]-1;
      if (i_xx[i][k]+w_x_size_2 < x_size) 
	i_xx_max[i][k]= w_x_size_2;
	
      //shift verticale negativo attorno al punto di indice i_zz[i][k] per costruire la finestra in z
      i_zz_min[i][k]=i_zz[i][k];
      if (i_zz_min[i][k] - w_z_size_2 > 0) 
	i_zz_min[i][k] = w_z_size_2;
	
      //shift verticale positivo attorno al punto di indice i_zz[i][k] per costruire la finestra in z
      i_zz_max[i][k]=z_size-i_zz[i][k]-1;
      if (i_zz[i][k]+w_z_size_2 < z_size) 
	i_zz_max[i][k]= w_z_size_2;
	
      //indici minimo e massimo in x e z per definire la finestra sul punto	
      im[i][k]=i_xx[i][k]-i_xx_min[i][k];
      ix[i][k]=i_xx[i][k]+i_xx_max[i][k];
      jm[i][k]=i_zz[i][k]-i_zz_min[i][k];
      jx[i][k]=i_zz[i][k]+i_zz_max[i][k];
	
    }
  }
    
  /*
    ;------------------------------------------------------------------------------
    ;          FASE 2
    ;------------------------------------------------------------------------------
    ;   Costruzione matrice pesi
    ;   Questa matrice contiene i pesi (in funzione della distanza) per ogni punto.
    ;-----------------------------------------------------------------------------*/
    
  for (k=0;k<w_x_size;k++)
    w_x[k]=exp(-pow(k-w_x_size_2,2.)/pow(w_x_size_2/2.,2.)); 
  for (k=0;k<w_z_size;k++)
    w_z[k]=exp(-pow(k-w_z_size_2,2.)/pow(w_z_size_2/2.,2.));
    
  for (i=0;i<w_x_size;i++){
    for (j=0;j<w_z_size;j++){
      w_tot[i][j]=w_x[i]*w_z[j];
    }
  }
        
  /* ;----------------------------------------------------------- */
  /* ;   Matrici per puntare sul piano cartesiano velocemente */
    
  /* #ifdef TIME */
  /*     prendo_tempo(); */
  /* #endif */
    
  /* ;---------------------------------- *//* ;---------------------------------- *//* ;---------------------------------- */
  /* ;---------------------------------- */
  /* ;          FASE 3 */
  /* ;---------------------------------- */
  /* ; Selezione dati per formare RHI */
  /* ;---------------------------------- */
    
    
  for (i=0; i<NUM_AZ_X_PPI; i++){
    cil[i]=(float **) malloc (x_size*sizeof(float *));
    for (k=0; k<x_size; k++){
      cil[i][k]= (float *) malloc (z_size*sizeof(float));
      for (j=0;j<z_size;j++) 
	cil[i][k][j]= -20. ;    }
  }	
    
  rhi_cart= (float **) malloc (x_size*sizeof(float *));
  rhi_weight= (float **) malloc (x_size*sizeof(float *));      
  for (k=0; k<x_size; k++){
    rhi_cart[k]= (float *) malloc (z_size*sizeof(float ));
    rhi_weight[k]= (float *) malloc (z_size*sizeof(float ));
  }            
  for(k=0;k<MAX_BIN;k++){
    beamXweight[k]=(float **) malloc(w_x_size*sizeof(float *));
    for(i=0;i<w_x_size;i++)
      beamXweight[k][i]=(float *) malloc(w_z_size*sizeof(float));
  }
      
  for (iaz=0; iaz<NUM_AZ_X_PPI; iaz++){
    for (k=0; k<x_size; k++){
      for (j=0; j<z_size; j++){
	rhi_cart[k][j]=0. ;
	rhi_weight[k][j]=0. ;
      }		    
    }
   
    for (i=0;i<n_elev;i++){
      // vol_pol[i][iaz].b_header.max_bin sostituito da vol_pol[0][iaz].b_header.max_bin	  
      for (k=0;k<vol_pol[i][iaz].b_header.max_bin;k++){	 
	RHI_beam[i][k]=BYTEtoDB(vol_pol[i][iaz].ray[k]);	  
      }
      for (k=vol_pol[0][iaz].b_header.max_bin;k<MAX_BIN;k++)
	RHI_beam[i][k]=BYTEtoDB(0);
    }
	
    /* ;---------------------------------- */
    /* ;          FASE 4 */
    /* ;---------------------------------- */
    /* ;   Costruzione RHI */
    /* ;---------------------------------- */
   
    for (iel=0;iel<n_elev;iel++){
      for (ibin=0;ibin<vol_pol[0][iaz].b_header.max_bin;ibin++) {
     
	for(kx=0;kx<w_x_size;kx++){
	  for(kz=0;kz<w_z_size;kz++){
	    beamXweight[ibin][kx][kz]=RHI_beam[iel][ibin]*w_tot[kx][kz];	 
	  }	      		
	}
      }

      for (ibin=0;ibin<vol_pol[0][iaz].b_header.max_bin;ibin++) {
	imin=im[ibin][iel];
	imax=ix[ibin][iel];
	jmin=jm[ibin][iel];
	jmax=jx[ibin][iel];
     
	wimin=w_x_size_2-i_xx_min[ibin][iel];
	wimax=w_x_size_2+i_xx_max[ibin][iel];
	wjmin=w_z_size_2-i_zz_min[ibin][iel];
	wjmax=w_z_size_2+i_zz_max[ibin][iel];	       
	for (i=imin;i<=imax;i++) {
	  for (j=jmin;j<=jmax;j++) {		
	    rhi_cart[i][j]=rhi_cart[i][j] +  beamXweight[ibin][wimin+(i-imin)][wjmin+(j-jmin)];
	    rhi_weight[i][j]=rhi_weight[i][j]+w_tot[wimin+(i-imin)][wjmin+(j-jmin)];
	  }
	}  
      }  
    }
    for (i=0;i<x_size;i++) {    
      for (j=0;j<z_size;j++) {
	if (rhi_weight[i][j] > 0.0) 
	  rhi_cart[i][j]=rhi_cart[i][j]/rhi_weight[i][j];
	else {
	  rhi_cart[i][j]=missing_value;	 

	} 
	cil[iaz][i][j]=rhi_cart[i][j];
	jbb=ceil((htbb)/RES_HOR);       
	if (j == jbb ) cappi[iaz][i] = DBtoBYTE(cil[iaz][i][j]);
       
      }  
    }   
  }
  
  //output = fopen("CAPPI","w");
  //fwrite(cappi,sizeof(cappi),1,output);
  //fclose(output);


  classifico_VIZ();
  if (hmax > 2000.) {// per evitare contaminazioni della bright band
    calcolo_background();
    classifico_STEINER();
  }
  merge_metodi();
  //creo_cart();
  //creo_cart_z_lowris();
  // classifico_ZLR();

  //output = fopen("CAPPI_cart","w");
  //fwrite(z_out,sizeof(z_out),1,output);
  //fclose(output);

  //output = fopen("CONV_cart","w");
  //fwrite(conv_1x1,sizeof(conv_1x1),1,output);
  //fclose(output);

  return ;
}

/*===============================================*/
void 	classifico_VIZ()
/*===============================================*/
{
  int i,j,k,kbb,ktbb;
  float cil_Z,base;
  double *Zabb[NUM_AZ_X_PPI], ext;
  float LIM_VERT= 8.;//questo l'ho messo io
  FILE *file_Zabb;
  ///* DETERMINAZIONE DEL VIL*/
  ///* determino top bright band  htbb*/
  ///* determino top echo..come? 0 dBZ?*/
  ///* ciclo su volume cilindrico a partire da quota htbb fino al top dell' eco , calcolo Zabb e se > thr segno Convettivo*/
   
  kbb=ceil(htbb/resol[1]);
  ktbb=ceil(LIM_VERT/resol[1]); 
  if (t_ground < T_MAX_ML) ktbb=0;/////se t suolo dentro t melting layer pongo ktbb=00 e in tal modo non classifico
  printf(" %i\n", kbb);
  if (kbb>z_size) kbb=z_size;
  printf ("kbb= %i \n  z_size= %i \n ",kbb,z_size);
  for (k=0; k<NUM_AZ_X_PPI; k++){
    Zabb[k]= (double *) malloc(x_size*sizeof(double )); 
    conv_VIZ[k]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;
    conv_STEINER[k]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;
    conv[k]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;	
    for (j=0; j<x_size; j++){
      Zabb[k][j]=0.;
      conv_VIZ[k][j]=MISSING;
      conv_STEINER[k][j]=MISSING;
      conv[k][j]=MISSING;
    }
  }	
  
  for(i=0; i<NUM_AZ_X_PPI; i++){
    for(j=0; j<x_size; j++)
      {
	ext=0.;
	//	    for(k=kbb; k<z_size; k++)
	for(k=kbb; k<ktbb; k++)
	  {
	    if (cil[i][j][k] > -20.){	
	      base=(cil[i][j][k])/10.;
	      cil_Z=pow(10.,base);
	      Zabb[i][j] = Zabb[i][j] + resol[1]*cil_Z;
	      ext=resol[1]+ext;
	      //if (j == 400) printf(" %i  %i %f  %f  %f \n", i, j, cil[i][j][k],Zabb[i][j],ext );
	    }
	  }
	if (ext>0.5) {
	  Zabb[i][j]=Zabb[i][j]/ext;
	  class_VIZ[i][j]=1;

	  if (Zabb[i][j] > THR_ZABB){
	    conv_VIZ[i][j]=CONV_VAL;
	    lista_conv[ncv][0]= i;
	    lista_conv[ncv][1]= j;	      
	    ncv=ncv+1;
	  }	     
	}	    	   
      }
  } 
  
  printf("numero nuclei VIZ = %i \n",ncv);  
  file_Zabb= fopen("file_Zabb","w");
  for(i=0; i<NUM_AZ_X_PPI; i++)
    for(j=0; j<x_size; j++)
      {
	fprintf(  file_Zabb, " %f   ",Zabb[i][j] );
	fprintf(  file_Zabb, " \n   ");	  
      }	  
  fclose(file_Zabb);
  return;
}

void 	classifico_STEINER() 
{   
  int i,j,k;
  
  unsigned char BYTE;
  float diff_bckgr;
 
  for(i=0; i<np-1; i++){
    j=lista_bckg[i][0];
    k=lista_bckg[i][1];
    BYTE=vol_pol[elev_fin[j][k]][j].ray[k];
    diff_bckgr=BYTEtoDB(BYTE)-bckgr[i];
    if (((BYTEtoDB(BYTE)>40.)||
	 (bckgr[i]<0 && diff_bckgr > 10) ||
	 (bckgr[i]<42.43 &&  bckgr[i]>0 &&  diff_bckgr > 10. - bckgr[i]*bckgr[i]/180. )||
	 (bckgr[i]> 42.43 &&  diff_bckgr >0)) && (class_VIZ[j][k]==0) ) {
   
      conv_STEINER[j][k]=CONV_VAL;
      ingrasso_nuclei(convective_radius[i],j,k);
      ncs=ncs+1;


    }
  }
  return; 
} 

void calcolo_background() // sui punti precipitanti calcolo bckgr
{
  int i,j,k,jmin,jmax,J_MAX,J_BBB,J_ABB,kmin,kmax,d_az,npoints=0; 
  float z,zmax,range,delta_r;
  int N_PUNTI_MIN=100;
  
  // calcolo la J max sotto bright band
  for (j=0;j<x_size;j++){
    range=RES_HOR/2+j*RES_HOR;
    z=pow(pow(range,2.)+pow(4./3*REARTH,2.)+2*range*4./3.*REARTH*sin(elev_array[0]*CONV_RAD),.5) -4./3.*REARTH+h_radar;
    zmax=pow(pow(range,2.)+pow(4./3*REARTH,2.)+2*range*4./3.*REARTH*sin(elev_array[NEL-1]*CONV_RAD),.5) -4./3.*REARTH+h_radar;
    if (z<hbbb) J_BBB=j;
    if (zmax<htbb) J_ABB=j;
  } 
  if (J_ABB< J_BBB ) J_MAX=J_ABB;
  else J_MAX=J_BBB;
  delta_r=10.;// per permettere confronto con background
  J_MAX=J_MAX+delta_r/RES_HOR;
  
  // J_MAX=x_size; // torno al controllo sul campo totale  
  // algoritmo 
  for (i=0; i<NUM_AZ_X_PPI;i++)
    for (j=0; j<J_MAX;j++)   
      {
	if ( vol_pol[0][i].ray[j] > DBtoBYTE(-19.) && (float)(quota[i][j])/1000. < hbbb ) // fisso una soglia minima per calcolare il back soglia=-19 DB e prendo punti sotto bright band
	  {
	    lista_bckg[np][0]=i;
	    lista_bckg[np][1]=j;
	    np=np+1;
	  } 
      } 
  convective_radius=(float *)malloc(np*sizeof(float));
  printf("numero nuclei Steiner = %i \n",np);
  fprintf(log_vpr,"numero nuclei Steiner = %i \n",np);
  if (np > N_PUNTI_MIN) {
    Z_bckgr=(double *) malloc(np*sizeof(double));
    bckgr=(float *) malloc(np*sizeof(float));
    for (i=1;i<np-1;i++) bckgr[i]=0.;
    for(i=0; i<np-1;i++){
      npoints=0;
      // estremi della finestra in range
      kmin=lista_bckg[i][1]-22;
      kmax=lista_bckg[i][1]+22;
      if (kmin<0) kmin=0;  
      if (kmax>x_size) kmax=x_size;    
      // estremi della finestra in azimut (qui devo fare un calcolo)
      if (lista_bckg[i][1]>0) d_az=ceil(11./(lista_bckg[i][1]*0.9));/// da mettere parametrizzato
      else d_az=0;
      
      jmin=lista_bckg[i][0]-d_az;
      jmax=lista_bckg[i][0]+d_az;
      ///  capire come gestire il minimo e il massimo in azimut quando va sottozero
      if (jmin<0) {
	jmin= NUM_AZ_X_PPI+jmin;
	for (j= jmin  ; j< NUM_AZ_X_PPI ; j++)
	  for (k= kmin ; k< kmax  ; k++){
	    if (BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k])>BYTEtoDB(1)) { 
	      Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]) ;
	      npoints=npoints+1;
	    }
	  }
	jmin=0; 
      }
      
      if (jmax>NUM_AZ_X_PPI) {
	jmax=jmax-NUM_AZ_X_PPI;
	for (j= 0  ; j< jmax ; j++)
	  for (k= kmin ; k< kmax  ; k++){
	    if (BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k])>BYTEtoDB(1)){ 
	      Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);  
	      npoints=npoints+1;
	    }
	  }
	jmin=0; 
      }
      
      if (jmax>NUM_AZ_X_PPI) {
	jmax=jmax-NUM_AZ_X_PPI;
	for (j= 0  ; j< jmax ; j++)
	  for (k= kmin ; k< kmax  ; k++){
	    if (BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k])>BYTEtoDB(1)){ 
	      Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);
	      npoints=npoints+1;}
	  } 
	jmax=NUM_AZ_X_PPI; 
      }
      for (j=jmin   ; j<jmax  ; j++)
	for (k=kmin  ; k<kmax   ; k++){
	  if (BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k])>BYTEtoDB(1)){ 
	    Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);
	    npoints=npoints+1;
	  }
	} 
      Z_bckgr[i]=Z_bckgr[i]/npoints;
      if (Z_bckgr[i]>0) bckgr[i]=10.*log10(Z_bckgr[i]);
      
      if (  bckgr[i] < 25.) convective_radius[i] = 1.;
      if (  bckgr[i] >= 25. && bckgr[i] <30. ) convective_radius[i] = 2.;
      if (  bckgr[i] >= 30. && bckgr[i] <35. ) convective_radius[i] = 3.;
      if (  bckgr[i] >= 35. && bckgr[i] <40. ) convective_radius[i] = 4.;
      if (  bckgr[i] > 40.)  convective_radius[i] = 5.;
    }  
  }
  return;
}

void ingrasso_nuclei(float cr,int ja,int kr) // per rimpolpare i nuclei convettivi uso l'algoritmo di Steiner:
// a tale scopo devo calcolare la Z di background per ogni cella convettiva
{ 
  int daz, dr,jmin,jmax,kmin,kmax,j,k;

 

  daz=ceil(cr/(ja*0.9));/// da mettere parametrizzato
  jmin=ja-daz;
  jmax=ja+daz;
  dr=ceil(cr/RES_HOR);
  kmin=kr-dr;
  kmax=kr+dr;  
  if (kmin<0) kmin=0;  
  if (kmax>x_size) kmax=x_size; 

  if (jmin<0) {
    jmin= NUM_AZ_X_PPI+jmin;
    for (j= jmin  ; j< NUM_AZ_X_PPI ; j++)
      for (k= kmin ; k< kmax  ; k++) conv_STEINER[j][k]=CONV_VAL;
    jmin=0; 
  }    
  if (jmax>NUM_AZ_X_PPI) {
    jmax=jmax-NUM_AZ_X_PPI;
    for (j= 0  ; j< jmax ; j++)
      for (k= kmin ; k< kmax  ; k++) conv_STEINER[j][k]=CONV_VAL;
    jmax=NUM_AZ_X_PPI; 
  }
  for (j=jmin   ; j<jmax  ; j++)
    for (k=kmin  ; k<kmax   ; k++) conv_STEINER[j][k]=CONV_VAL;

  return;
}

void merge_metodi()
{
  int j,k;
  
  for(j=0; j<NUM_AZ_X_PPI; j++){
    for(k=0; k<x_size; k++)
      {
	if (conv_STEINER[j][k] >= conv_VIZ[j][k]) {
	  conv[j][k]=conv_STEINER[j][k];	  
	}	
	else  
	  {
	    conv[j][k]=conv_VIZ[j][k]+(unsigned char)(50);  	  
	  }

      }
  }
  return;
}
#endif


int trovo0term()
{
  int ier;
  FILE *file0;
  ier=access(getenv("FILE_ZERO_TERMICO"),R_OK);
  if (!ier) { 
    file0=fopen(getenv("FILE_ZERO_TERMICO"),"r");//metti condizione di non trovato!!!
    fscanf(file0,"%f", &zeroterm );   
    fclose(file0);
  } 
  else {
    printf("non ho trovato lo zero termico \n ");
  }
  return ier;
} 
