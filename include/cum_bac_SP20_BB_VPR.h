#ifndef ARCHIVIATORE_CUM_BAC_H
#define ARCHIVIATORE_CUM_BAC_H
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
#include <sys/types.h> 
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
#include <setwork.h>
#include <func_Q3d.h>

/**
 *  
 *  @brief funzione che elabora il dato radar rimuovendo anaprop e beam blocking
 *  @details  partendo dal livello della mappa dinamica esegue il test di continuità verticale e laddove non si verifica un cambio di elevazione corregge il beam bocking per ottenere infine un campo bidimensionale adeguato alla stima di R che ricopia su tutti i livelli del volume  a partire dallo 0 fino al livello della mappa dinamica . Memorizza elevazione finale usata per il campo bidimensionale e l'output della rimozione della propagazione anomala e quota al livello scelto per la stima di R.
 *  
 *  @return 0 se ok 1 se errore
*/ 
int elabora_dato(); 


//int trovo_top();

// lettura mappa statica 

/**
 *  
 *  @brief funzione che legge la mappa statica e la mappa di elevazioni da beam blocking e le condensa in un unica mappa
 *  @details crea una mappa di elevazioni scelte in ogni pixel del volume per ottenere successivamante la stima di R
 *  @return non ritorna nulla
*/
void leggo_first_level();

  /**
 *  
 *  @brief funzione che a partire dal tempo in secondi arrotonda al NMIN-esimo minuto precedente o successivo
 *  @details 
 *  @param  time intero che rappresenta il numero di secondi
 *  @return restituisce l'arrotondamento in secondi al NMIN-esimo minuto precedente o successivo , se fallisce -1
*/  
 
time_t 	NormalizzoData(time_t  time); 

 /**
 *  
 *  @brief funzione che scrive il tempo parziale inercorso dalla chiamata precedente e il totale dall'inizio del programma
 *  @return
*/
void prendo_tempo(); 

 /**
 *  @brief   funzione scrittura matrici statistica
 *  @details scrive le statistiche di beam blocking, anaprop, cambio di elevazione in un unsigined char DIM1_ST*DIM1_ST
 * 
 *  @return non ritorna valori
 *  
*/
void ScrivoStatistica();

/**
 *  @brief funzione che verifica se il volume dati in input e' consistente con le esigenze del programma 
 *  @details  TESTA: A) risoluzione B) numero raggi per elevazione
 *  @param    tipofile  puo' essere 0(corto senza declutter) 1(corto con declutter) 2(short hail) o 3(medio)
 *  @return    0  in caso di errore   1  in caso di successo  .
 *  
*/
int test_file(char*); 
// funzioni di conversione cartesiana:associa a pixel matrice alta ris azimut e range, crea alta risoluzione e crea bassa risoluzione  
  /**
 *  
 *  @brief funzione che calcola range e azimut su un quadrante centrato in 0,0
 *  @details 
 *  @return non ritorna nulla
*/ 

void creo_matrice_conv(); 
 /**
 *  
 *  @brief funzione che crea l'output cartesiano dal polare
 *  @details cicla sui quadranti e su i e j, usando il range e l'azimut ottenuti tramite la funzione creo_matrice_conv()
 *  @return 
*/ 
void creo_cart();// conversione da polare a cartesiano alta risoluzione 
//#ifdef Z_LOWRIS 
 /**
 *  
 *  @brief funzione che  trasforma il dato da cartesiano alta risoluzione a cartesiano bassa risoluzione
 *  @details prende il massimo tra i punti considerati, il passo di ricerca è ZLR_N_ELEMENTARY_PIXEL, cioè il rapporto tra dimensioni ad alta risoluzione e le dimensioni bassa risoluzione.
 *  @return
*/   

void creo_cart_z_lowris(); 
// funzione di scrittura matrici output binario
 /**
 *  
 *  @brief funzione che scrive in un file di output una matrice di byte di dimensione size
 *  @details il formato del nome e' $dir/aaammgghhmm_$ext 
 *  @param[in] ext estensione file output (aaammgghhmm_$ext)
 *  @param[in] content contenuto del file 
 *  @param[in] dir directory dove scrivere il file 
 *  @param[in] size dimensione della matrice
 *  @param[in] matrice matrice di dati da scrivere
 *  @return 
*/ 
void scrivo_out_file_bin (char *ext,char *content,char *dir,size_t size, void  *matrice);
//#endif 

// funzione di controllo esistenza file e restituisce puntatore a file+messaggio ev errore o info su status apertura, utile? 

//funzioni associate a logging , vedere se togliere 
  /**
 *  
 *  @brief funzione che restituisce l'orario 
 *  @return asctime(tempo) ritorna l'orario in formato asctime
*/
char *PrendiOra();//togliere? 

// funzioni per QUALITA' 
//#ifdef QUALITY
//lettura dem e quote fascio 
 /**
 *  
 *  @brief funzione che legge il dem
 *  @details legge il dem e lo memorizza nella matrice dem
 *  @return
*/
void leggo_dem();  
  /**
 *  
 *  @brief funzione che legge la quota del centro fascio e del limite inferiore del fascio da file
 *  @details legge la quota del centro fascio e del limite inferiore del fascio da file e li memorizza nei vettori hray_inf e hray
 *  @return 
*/  
      
void leggo_hray(); 

/**
 *  @brief   funzione che calcola l'att enuazione totale 
 *  @details Ricevuto in ingresso il dato di Z in byte e l'attenuazione complessiva sul raggio fino al punto in considerazione, calcola l'attenuazione totale 
 *  @param[in]  DBZbyte valore in byte della Z nel pixel
 *  @param[in]  PIA attenuazione totale fino al punto
 *  @return att_tot attenuazione incrementata del contributo nel pixel corrente
 *  
*/
double attenuation(unsigned char DBZbyte, double  PIA); //ingresso:dbz in quel punto e attenuazione fin lì

/**
 *  @brief funzione che caratterizza i volumi polari tramite la qualita'
 *  @details utilizza i parametri di ouput dell'elaborazione e la PIA calcolata qui per calcolare un valore finale di qualita' del dato. Inoltre calcola il top dell'eco in base a soglia su ogni pixel.
 *  @param  sito identificativo del radar (spc o gat)
 *  @return non ritorna nulla
 *  
*/
void caratterizzo_volume();

/**
 *  
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
int func_vpr(); // crea vpr istantaneo 
  /**
 *  
 *  @brief funzione che corregge per il profilo verticale
 *  @details ciclando su tutti i bins della cartesiana polare scelta per la stima della pioggia,
 * - trovo la quota del centro del pixel
 * - correggo se: hbin>hliq , valore maggiore di soglia correzione (0 dBZ)
 * - se sefinita CLASS e pixel convettivo non correggo 
 * @param[in] sito radar su cui calcolare profilo
 * @return 0 se ok 1 se fallisce
*/
int corr_vpr(char*); // correzione vpr 
 /**
 *  
 *  @brief funzione che  calcola il no di quarti d'ora che intercorrono dall'ultimo profilo calcolato (combinato)  memorizzato in 'LAST_VPR'  
 *  @param[in] nomefile nome del file LAST_VPR dove c'e' la data cui si riferisce l'ultimo profilo prodotto in n0 di secondi a partire da istante di riferimento
 *  @return gap1 ritorna il  no di quarti d'ora che intercorrono dall'ultimo profilo calcolato
*/ 
long int profile_gap(); // calcolo gap in qualti d'ora 
 /**
 *  
 *  @brief funzione che calcola quanto il profilo è 'caldo' restituendo la variabile heating
 *  @details calcola il 'riscaldamento' del profilo (numero di combinazioni successive) e stampa il file che contiene questo valore oltre al file contenente l'ultima data se il calcolo del profilo è andato ok. Se il numero è superiore a WARM passa direttamente al valore MEMORY che rappresenta la 'memoria del profilo'.
 *   heating=heating-gap; se il profilo non è stato aggiornato 
 *   heating=heating-gap+2; se il profilo è stato aggiornato 
 *   heating=MEMORY;   se heating raggiunge WARM resta costante finchè non inizia raffreddamento 
 *  @return heating , numero di combinazioni di riscaldamento
*/ 
int profile_heating(); //calcola riscaldamento in quarti d'ora 
 /**
 *  
 *  @brief funzione che combina il profilo verticale corrente con quello precedente tramite il metodo di Germann
 *  @details  oltre a lanciare il calcolo del profilo istantaneo provvede alla combinazione del profilo calcolato con il precedente calcolato entro  un limite massimo di distanza temporale pari a 10 quarti d'ora.  restituisce un codice integer pari a 0 se ok 1 se fallisce 
 *  @param[in] sito radar corrente 
 *  @return 0 se combinazione ok 1 se fallisce
*/ 

int combina_profili(char *sito);//combina profili
  /**
 *  
 *  @brief funzione che calcola quanto il profilo è 'caldo' restituendo la variabile heating
 *  @details calcola il 'riscaldamento' del profilo (numero di combinazioni successive) e stampa il file che contiene questo valore oltre al file contenente l'ultima data se il calcolo del profilo è andato ok. Se il numero è superiore a WARM passa direttamente al valore MEMORY che rappresenta la 'memoria del profilo'
 *  @return heating , numero di combinazioni di riscaldamento
*/ 

int stampa_vpr();// stampa profilo combinato 
 /**
 *  @brief   funzione che restituisce la temperatura al suolo
 *  @details  apre file temperature , legge lon lat e t, calcola differenze rispetto coordinate radar, se diff < soglia media il dato, stampa il nr di dati usati per la media e ritorna la temperatura
 *  @param[out] t_gr temperatura al suolo
 *  @return ierr codice di uscita (0=ok 1=fallito)
 *  
*/ 
int get_t_ground(float *t_gr);//fornisce temperatura al suolo, da lettura file esterno 

/**
 *  @brief   funzione che analizza il profilo 
 *  @details analizza il profilo usando : la temperatura al suolo, la quota del massimo, e una funzione di interpolazione 
 *  @param[out]  vpr_liq valore del profilo al livello liquido
 *  @param[out]  snow matrice che indica se è presente neve o no secondo l'analisi fatta
 *  @param[out]  hliq quota del livello liquido
 *  @param[in]   sito radar in esame (spc o gat)
 *  @return ier_ana valore che indica se tutto è andato a buon fine (0) o no (1)
 *  
*/ 
int analyse_VPR(float *vpr_liq,int *snow,float *hliq, char *sito);
 
/** 
 *  @brief funzione che trova la quota del massimo valore del profilo
 *  @details trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca come massimo di almeno 5 dBZ più alto in quota
 *  @param[in] sito radar su cui calcolare profilo
 *  @return 0 se ok 1 se fallisce
 */ 
int trovo_hvprmax(int *hmax); // trova il massimo del profilo 

 /**
 *  @brief   funzione che esegue interpolazione del profilo
 *  @details interpola profilo usando una funzione gaussiana+lineare  y= B*exp(-((x-E)/G)^2)+C+Fx 
 *  @param[out] a[] vettore dei parametri della funzione
 *  @param[out] ma  dimensione vettore parametri
 *  @return ier_int codice di uscita (0=ok 1=fallito)
 *  
*/ 
int interpola_VPR(float a[], int ma); // interpola profilo 

/**
 *  @brief   funzione che  calcola derivate della gaussiana + lineare rispetto ai parametri e valore (*y) in un punto x 
 *  @details    y(x,a) is the sum of a gaussian and a linear function with amplitude B=a[1], center E=a[2] and width G=a[3] and a linear function with coefficients C(shift)=a[4] and F(slope)=a[5]  
 *  a is the parameters vector and dyda is the vector of derivatives respect to the different parameters
 *  @param[in] x vettore delle quote
 *  @param[in] a vettore dei parametri
 *  @param[in] y vettore dei valori della funzione 
 *  @param[out]  dyda derivate
 *  @return 0 codice di uscita 0
 *  
*/  
void lineargauss(float x, float a[], float *y, float dyda[],int na);// restituisce il valore della combinazione tra gaussiana e lineare

int n_close_(); //? 
//#endif
 
//funzioni per classificazione 

//#ifdef CLASS 
 
/**
 *  @brief funzione che legge il file contenente la quota dello 0 termico 
 *  @return ier codice di errore apertura file
 */
int trovo0term();//trova lo zero termico 
 
/**
 *  
 *  @brief funzione  che classifica secondo il metodo VIZ
 *  @details calcolo per ogni pixel polare l'integrale verticale esclusa la fascia della bright band
 * @return
*/ 
void classifico_VIZ();// classifica tramite Vertical Integrated Z 

/**
 *  
 *  @brief funzione  che classifica secondo STEINER
 *  @details segna come convettivi i punti che hanno valore superiore a 40 dBZ e differenza col background elevata, quindi ingrandisce i nuclei di un raggio variabile
 * @return
*/ 
void classifico_STEINER();//  classifica con metodi di Steiner
 /**
 *  
 *  @brief funzione  che calcola il background 
 *  @details la classificazione di Steiner non ha bisogno di ricampionamento cilindrco perciò  uso direttamente la matrice polare
 * @return 
*/  
void calcolo_background(); // calcola valore di background per individuare pixel convettivo
 /**
 *  
 *  @brief funzione  che ingrandisce i nuclei di Steiner
 *  @details ingrandisce i nuclei di Steiner di un valore pari al raggio convettivo
 *  @param[in] cr raggio convettivo
 *  @param[in] ja indice di azimut
 *  @param[in] kr indice di range
*/ 
void ingrasso_nuclei(float cr,int j,int k);// ingrassa nuclei convettivi 

 /**
 *  
 *  @brief funzione  che interseca i punti convettivi delle due classificazioni Viz e Steiner e sottrae quelli con  picco stratiforme
 *  @return non ritorna valori
*/ 
void merge_metodi();// fa il merge dei metodi

// funzione che calcola la qualita' per R
float func_q_R(unsigned char cl, unsigned char bb, float dst, float dr, float dt, float dh, float dhst, float PIA, float dZ, float sdevZ);

     
#endif //  ARCHIVIATORE_CUM_BAC_H
