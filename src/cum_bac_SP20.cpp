/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief codice originale di elaborazione dei volumi di riflettivita' radar usato per impulso medio
 *  @details questo codice contiene il main e alcune funzioni usate nell'elaborazione della riflettivita' radar 
 *  e rappresenta il programma originale usato solo per l'elaborazione dell'impulso medio
*/

/*--------------------------------------------------------------------------------
 | CUMULO SU BACINI
 | 
 | Alberoni Pier Paolo  S.M.R. 
 |--------------------------------------------------------------------------------
 |cc cum_bac_SP20.c  -lm  -DANAPROP  -DTIME -DZ_AVERAGE
 |                 -DZ_LOWRIS -DSIDME -DSHORT -I$HOME/include -lSP20_utility -L$HOME/lib       
 | MODIFICHE
 |
 |
 | 31-07-1996 - Aggiunta la costruzione e la scrittura di una matrice di 
 |              riflettivita' alla risoluzione di 1Kmx1Km in output
 |              l'opzione e' attivabile tramite la dichiarazione di una variabile 
 |		in fase di compilazione
 |
 | 06-08-1996 - Scrittura file di log
 |
 | 15-10-1996 - trovato errore nell'estrazione precipitazione sui bacini
 |		la matrice era stata trasposta
 |
 | 02-12-1996 - Aggiunta la routine write_xdr che scrive il vettore delle precipitazioni
 |   		medie sui bacini in formato xdr.
 |
 | 07-02-1997 - Ristretto il casi su cui viene fatto il test di anaprop:
 |              solo se (bin_low > fondo_scala && bin_high >= fondo_scala)
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
 |              questo permette di attivare la rimozione del. clutter e della propagazione anomala,
 |              del calcolo della cumulata a livello di bacino e permette di gestire i file a 
 |              impulso corto e impulso medio.
 |              Vedi sotto per una descrizione completa delle variabili di compilazione
 !
 |--------------------------------------------------------------------------------
 | Questo programma legge i file polari DATAMAT MDB_2.0
 | scopo del programma e' eliminare il clutter tramite una maschera e 
 | implementare un algoritmo di rimozione della propagazione anomala
 | come presentato al RADME 96.
 | Il programma provvede inoltre a convertire il PPI, ripulito, a minor
 | elevazione in formato cartesiano ,estrarre se presente la precipitazione
 | sui bacini regionali e a scrivere una matrice di precipoitazione media 
 | su un grigliato di 5x5 Km.
 | I file polari prima di essere utilizzati viene verificata la consistenza
 | con i seguenti criteri :
 |  - deve essere presente la rimozione del clutter tramite filtro doppler
 |  - devono essere presenti almeno NUM_MIN_BEAM per ognuna delle prime
 |    4 elevazioni
 |  - il file deve essere piu' recente dell'ultimo file processato
 |    ( questa condizione e' stata aggiunta per facilitare l'utilizzo 
 |       operativo a S.P.C.)
 |--------------------------------------------------------------------------------
 | 
 | Vengono utilizzate le seguenti variabili di ambiente :
 |
 | LISTA_FILE          - nome del file contenete i file dati da utilizzare
 | LAST_FILE           - nome del file contenente l'ultima data processata
 | ANAP_STAT_FILE      - nome del file contenente la statistica sull'anaprop
 | FIRST_LEVEL_FILE    - nome del file contenente la mappa del primo livello 
 |                       sopra il clutter
 | OUTPUT_Z_DIR        - directory di output per i file del grigliato 5X5 Km 
 |                       mediati in riflettivita'
 | OUTPUT_RAIN_DIR     - directory di output per i file del grigliato 5X5 Km
 |                       mediati in precipitazione
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
 | NOME_SIDME          - File di output per procedura sidme
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
 |               nella porcessazione
 |
 | WRITE_DBP_REORDER -  permette di scrivere in output il volume polare
 |                      dopo la fase di rilettura e di riordino ad azimut
 |                      fissi ed elevazioni fisse
 |
 | DECLUTTER -  permette di scrivere in output il volume polare
 |              dopo la fase di rimozione del clutter senza attivare la 
 |              rimozione della propagazione anomala
 |
 | Z_AVERAGE -  abilita l'elaborazione e la scrittura dei dati sulla matrice
 |	        5X5 mediata in riflettivita'
 |
 | R_AVERAGE -  abilita l'elaborazione e la scrittura dei dati sulla matrice
 |	        5X5 mediata in precipitazione
 |
 | Z_LOWRIS  -  abilita la costruzione e la scrittura di una matrice di
 |              riflettivita' a bassa risoluzione (1km)
 |
 | SIDME     -  abilita la scrittura del file per procedura sidme, deve essere definita
 |              la variabile d'ambiente NOME_SIDME
 |
 | BACINI    -  abilita il calcolo e las crittura in output della cumulata a livello di bacini
 |
 | SHORT     -  configura il programma all'esecuzione per dati provenienti da impulso corto (125 km range; 250 m ris)
 |
 | MEDIUM    -  configura il programma all'esecuzione per dati provenienti da impulso medio (250 km range; 1000 m ris)
 |
 |--------------------------------------------------------------------------------
 |
 | Queste varialibi di compilazione vengono momentaneamente a perdere di validità 17-09-2003
 | 
 | Deve essere inoltre definita una delle seguenti variabili che individua
 | l'ambiente di esecuzione del programma e utilizza i file di include 
 | Datamat corretti:
 |
 | BOLOGNA - per compilare il programma a Bologna
 |
 | SPC     - per compilare il programma a San Pietro Capofiume
 |
 |--------------------------------------------------------------------------------
 |
 | In fase di link deve essere specificata la libreria matematica
 |
 |--------------------------------------------------------------------------------

----------------------------------------------------------------------------*/

#include <stdio.h>
#include <rpc/rpc.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <string.h>
#include <memory.h>
#include <unistd.h>
extern "C" {
#include <func_SP20read.h>
}
#include <radar_parameter.h>

/*#define NEL 10*/
#define MAX_BIN 512
#define MAX_DIM 512
#define FATT_MOLT_EL ((double) 360./(double)4096.)
#define FATT_MOLT_AZ ((double) 360./(double)4096.)

#ifdef SHORT
#define RANGE_KM                 125
#define SIZE_CELLA_KM            5
#define CART_DIM_ZLR             256
#define ZLR_N_ELEMENTARY_PIXEL   4
#define ZLR_OFFSET               0
#define NMIN 5
#define MAX_TIME_DIFF 3
int elev_array[NEL]={6,17,27,37,43,52,62,71,79,88};
#endif

#ifdef MEDIUM
#define RANGE_KM                 250
#define SIZE_CELLA_KM            5
#define CART_DIM_ZLR             512
#define ZLR_N_ELEMENTARY_PIXEL   1
#define ZLR_OFFSET               CART_DIM_ZLR/2
int elev_array[NEL]={6,15,24,34,57,52,62,71,79,88};
#define NMIN 1
#define MAX_TIME_DIFF 1
#endif

#define CART_DIM 2*RANGE_KM/SIZE_CELLA_KM

//   #include <datamat_file.h> 

/*----------------------------------------------------------------------------*/
/*	DICHIARATIVE GLOBALI						      */
/*----------------------------------------------------------------------------*/
T_MDB_ap_beam_header  old_beam_header;  // queste due dichiarazioni le metto qui perchè sono sfaticato
T_MDB_data_header   old_data_header;    

//extern char *sys_errlist[];
//extern int errno;
char errori[256];


#define MAX_DIF 30
#define MIN_VALUE -10
#define MAX_DIF_NEXT 15
#define MIN_VALUE_NEXT 0
#define MISSING 0
#define NUM_MIN_BEAM 200
#define PIOVE     1
#define NON_PIOVE 0
#define INIZIO_EVENTO PIOVE
#define FINE_EVENTO   NON_PIOVE
#define DIFF_TIME 1*60*60         /*--- tempo utile per chiudere evento ---*/
#define RAIN_MIN  0.5             /*--- soglia minima presenza rain ----*/
#define STEP_STAT_ANAP_RANGE  40
#define STEP_STAT_ANAP_AZ     25
#define N_MIN_BIN           500 /*--- numero minimo di celle presenti in un 
                                       settore per la statistica            ---*/

#define SHORT_DEC         0
#define SHORT_FULL_VOLUME 1
#define SHORT_HAIL        2
#define MEDIUM_PULSE      3 

/*-----------------------------------------------------------
  Variabili globali
------------------------------------------------------------*/


struct VOL_POL vol_pol[NEL][NUM_AZ_X_PPI]; 

unsigned char first_level[400][MAX_BIN];

float azimut[MAX_BIN][MAX_BIN];
float range[MAX_BIN][MAX_BIN];
unsigned char cart[MAX_BIN*2][MAX_BIN*2];
int stat_anap[16][13];
int stat_anap_tot[16][13];

#ifdef Z_AVERAGE
double cart_5x5_z[CART_DIM][CART_DIM];
#endif
#ifdef R_AVERAGE
double cart_5x5_r[CART_DIM][CART_DIM];
#endif
struct tm *time_dbp;

float size_cell[]={62.5,125.,250.,500.,1000.,2000.};
int nbeam_elev[NEL];

#ifdef Z_LOWRIS
unsigned char z_out[CART_DIM_ZLR][CART_DIM_ZLR];
#endif

/*----------------------------------------------------------------------------*/
/*	FUNCTION PROTOTYPE						      */
/*----------------------------------------------------------------------------*/
float DBZ(unsigned char z);
int remove_anap();
int write_dbp(const char* nome);
void leggo_first_level();
void creo_matrice_conv();
void creo_cart();
float mp_func(unsigned char byte);
void creo_cart_5x5(int imin,int imax,int jmin,int jmax,int dx,int dy);
void scrivo_5x5 (int imin,int imax,int jmin,int jmax,int dx,int dy);
float Z_func(unsigned char byte);
float Z_to_mp_func(double Z);
time_t 	NormalizzoData(time_t time);
int test_file();

void bacini();
void prendo_tempo();
void AggiornoHistoryFile(char* data, int flag);
void ScrivoBacini(const char* nome,char* data, float* array_bacini, int quanti);
void ScrivoStatistica();
#ifdef Z_LOWRIS
void creo_cart_z_lowris();
void scrivo_z_lowris();
#endif
char *PrendiOra();
void ScrivoLog(int i, char* stringa);
void   write_xdr(float* bacini, time_t Time);

/*sostituisco 
  unsigned int old_data_header.norm.maq.acq_date(secondi dal 1970) con    volume->getStartEpochs()?    
  unsigned int old_data_header.maq.resolution(indice di cell_size) con  .. non esiste il corrispondente ma c'è scan->getRangeScale() 
  unsigned int old_data_header.declutter_rsp(declutter o no)        dovrebbe essere sempre sì               
 
  unsigned char vol_pol[el][az].ray[MAX_DIM] dati del raggio (vettore) 	PolarScanData* data = scan->getQuantityData(quantityName);data->readTranslatedData(matrix);matrix.elem(ray, b) con b che va da 0 a data->getNumBins()
  char vol_pol[el][az].flag ??????????????????????
  short vol_pol[el][az].b_header.alfa azimuth in gradi/0.9 (0->399) scan->getAzimuthAngles();
  short vol_pol[el][az].b_header.teta elevazione in gradi*10? (0-->189) =scan->getElevationAngles(); si ottiene coppia (start, stop)
  short vol_pol[el][az].b_header.tipo_gran tipo grandezza ..mah  scan->getQuantityDataIndex(PRODUCT_QUANTITY_DBZH)
  short vol_pol[el][az].b_header.max_bin n0 bin per determinato raggiPolarScanData* data = scan->getQuantityData(PRODUCT_QUANTITY_DBZH);
		data->getQuantity()<<endl;	
		int bins = data->getNumBins();

  int nbeam_elev[] numero raggi per elevazione  con scan->getNumRays()

*/

/* ==================== */
int main (int argc, char **argv)
/* ==================== */
{
  char *nome_file , *nome_file_volume;
  int ier,ier_test;
  int tipo_dati_richiesti = INDEX_Z; 

  /*---------------------------------------
    definisco le variabili per estrarre 
    pezzi della matrice cartesiana
    (ad esempio 5x5)
    ---------------------------------------*/
  int imin,imax;
  int jmin,jmax;
  int dx,dy;

if (argv[1]==NULL) {
    printf("MAIN:argomento non passato \n");
    exit(1);
  }
 nome_file_volume=argv[1];
 nome_file=nome_file_volume;

  ScrivoLog(0,nome_file);
  ScrivoLog(1,nome_file);
  /*----------------------------------------
    | testo presenza e apro il file di lista |ELABORO UNO PER UNO
    ----------------------------------------*/
  /* ScrivoLog(2,nome_file); */
  /* if(access(getenv("LISTA_FILE"),R_OK) == -1) return; */
  /* file2 = fopen(getenv("LISTA_FILE"),"r"); */
  /* if(file2 == NULL) */
  /*   { */
  /*     ScrivoLog(3,nome_file); */
  /*     return; */
  /*   } */
  
  /* /\*----------------------- */
  /*   | inizio ciclo sui dati | */
  /*   -----------------------*\/ */
  /* while(1) */
  /*   { */
      /* if(fscanf(file2,"%s",nome_file) != 1) break; */
#ifdef TIME
      prendo_tempo();
#endif
      
      /* pulisco tutte le matrici */
      memset(vol_pol,0,sizeof(vol_pol));
      memset(cart,0,sizeof(cart));
#ifdef Z_AVERAGE
      memset(cart_5x5_z,0,sizeof(cart_5x5_z));
#endif
      memset(z_out,0,sizeof(z_out));
      memset(nbeam_elev,0,sizeof(nbeam_elev));
 
      ScrivoLog(4,nome_file);
      printf("processo file dati: %s\n",nome_file);
      ier = read_dbp_SP20(nome_file,vol_pol,&old_data_header,
			  tipo_dati_richiesti,nbeam_elev);
         printf("fatta lettura volume polare");
      ScrivoLog(5,nome_file);

      /*      stampa di controllo 
      for (i_el=0; i_el<NEL; i_el++)
	for (j_az=0; j_az<NUM_AZ_X_PPI; j_az++){
	  printf(" %2d - %3d -- %4d - %4d %4d - %1d",i_el,j_az,
		 vol_pol[i_el][j_az].b_header.teta,
		 vol_pol[i_el][j_az].b_header.alfa,vol_pol[i_el][j_az].flag,
		 vol_pol[i_el][j_az].b_header.max_bin);
	  for(k_ray=495; k_ray<500; k_ray++)
	    printf(" %3d",vol_pol[i_el][j_az].ray[k_ray]);
	  printf("\n");
	}
      */

#ifdef TIME
      prendo_tempo();
#endif
      ier_test=test_file();
      printf("ier -- test  %d  %d\n",ier,ier_test);
      //      PrintOldHeader(&old_data_header);

      if(ier == OK  && ier_test)
	{
#ifdef WRITE_DBP_REORDER
	  
	  /*-------------------------------------------------
	    | eventuale scrittura del volume polare ripulito e|
	    -------------------------------------------------*/ 
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
#ifdef ANAPROP
	      ier = remove_anap();  
#endif

#ifdef TIME
	      prendo_tempo();
#endif
	      
#ifdef WRITE_DBP
	      /*-------------------------------------------------
		| eventuale scrittura del volume polare ripulito e|
		-------------------------------------------------*/ 
#ifdef DECLUTTER
	      strcat(nome_file,"_decl");
#else
	      strcat(nome_file,"_anap");
#endif
	      //exit (1);
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
	      
#ifdef BACINI
	      /*------------------------------------
		| estraggo i dati sui singoli bacini |
		------------------------------------*/
	      ScrivoLog(10,nome_file);
	      bacini();
#endif
	      
#ifdef TIME
	      prendo_tempo();
#endif

	      /*--------------------------------
		| creo matrice con passo di 5 Km |
		--------------------------------*/
	      ScrivoLog(11,nome_file);
	      
	      imin = -RANGE_KM/SIZE_CELLA_KM;
	      imax = -imin-1;
	      jmin = -RANGE_KM/SIZE_CELLA_KM;
	      jmax = -jmin-1;
	      dx = (int) (SIZE_CELLA_KM*1000./size_cell[old_data_header.norm.maq.resolution]);
	      dy = (int) (SIZE_CELLA_KM*1000./size_cell[old_data_header.norm.maq.resolution]);
	      creo_cart_5x5(imin,imax,jmin,jmax,dx,dy);
	      
	      printf(" finita cart_5x5 \n");
#ifdef TIME
	      prendo_tempo();
#endif

	      /*--------------------------------------------
		| scrivo matrice con passo di 5 Km in uscita |
		--------------------------------------------*/
	      ScrivoLog(12,nome_file);
	      scrivo_5x5(imin,imax,jmin,jmax,dx,dy);
	      printf(" finita scrivo_5x5 \n");
	      
#ifdef TIME
	      prendo_tempo();
#endif
	      
#ifdef Z_LOWRIS
	      ScrivoLog(13,nome_file);
	      creo_cart_z_lowris();
	      printf(" finita z_lowris\n");
	      ScrivoLog(14,nome_file);
	      scrivo_z_lowris();
	      printf(" finita scrivo_z_lowris\n");
#endif
#endif
	    }
	  else
	    {
	      sprintf(errori,"Errore Normalizzo data ");
	      ScrivoLog(16,errori);	
	      exit(1);
	    }
	}
    /* }  COMMENTO FINE LOOP*/ 
  
  /*-----------------------------
    | fine loop sui volumi polari |
    -----------------------------*/
  /* fclose(file2);    */ 
  
  ScrivoLog(15,nome_file);
  return 0; 
}                                //fine del main

 





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
int test_file()
{
  FILE *f_aus;
  time_t last_time;
  int k;
  int n_elev,resolution;
  int file_type;

  sscanf(getenv("TIPO_FILE"),"%d",&file_type);
  printf("tipo file %1d\n",file_type);
  switch (file_type)
    {
    case SHORT_DEC:
      //printf("sono qua\n");  
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
    }                                                       // end switch

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
      printf("Trovati Pochi Beam Elevazione %2d - num.: %3d",k,nbeam_elev[k]);
      sprintf(errori,"Trovati Pochi Beam Elevazione %2d - num.: %3d",k,nbeam_elev[k]);

      ScrivoLog(16,errori);
      exit(1);
    }
  }                                                             //end for
      sprintf(errori,"Primi test passati");
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
          /*----------------------------
           |  aggiorno la data nel file |
            ----------------------------*/
	   else{
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
}                                                       //end funzione test_file()




float DBZ(unsigned char z)
{
  return (z*80./255.-20.);
}

/*===============================================*/
int remove_anap()
/*===============================================*/
{
  int i,l,k;
  float bin_low,bin_high;
  float fondo_scala;
  int el_inf,el_up;
  unsigned char flag_anap;

/*---------------------
 | calcolo fondo scala |
  ---------------------*/
  flag_anap = 1;
  fondo_scala = DBZ(flag_anap);

  leggo_first_level();
#ifdef DECLUTTER
  for(i=0; i<400; i++)
  {
      for(k=1; k<vol_pol[0][i].b_header.max_bin; k++)  
      {
	el_inf = first_level[i][k];
        for(l=0; l<el_inf; l++)
          vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
      }
   }
   return 0;
#endif

   /*---------------------------------------
     | azzero matrici per statistica anaprop |
     ---------------------------------------*/
  memset(stat_anap_tot,0,sizeof(stat_anap_tot));
  memset(stat_anap,0,sizeof(stat_anap_tot));

  for(i=0; i<400; i++)
  {
      flag_anap = 0;
      for(k=1; k<vol_pol[0][i].b_header.max_bin; k++)   
      {
        stat_anap_tot[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
	el_inf = first_level[i][k];
	el_up = el_inf +1;
        bin_low  = DBZ(vol_pol[el_inf][i].ray[k]);
        bin_high = DBZ(vol_pol[el_up][i].ray[k]);
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
	  |	if(bin_low != -20. && bin_high != -20 )   |
	  ------------------------------------------------*/
	/*printf(" %6.1f %6.1f %6.1f %1d %2d %2d",
	bin_low,bin_high,bin_low-bin_high,flag_anap,el_inf,el_up);*/

	if(bin_low > fondo_scala && bin_high >= fondo_scala )
	{
          if(flag_anap)
	  {
	    if(bin_low-bin_high >= MAX_DIF_NEXT || bin_high <= MIN_VALUE_NEXT) 
	    {
	      for(l=0; l<el_up; l++)
	        vol_pol[l][i].ray[k]=1;
	      flag_anap = 1;
              stat_anap[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++; 
	    }
	    else 
	    {
	      for(l=0; l<el_inf; l++)
	        vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
	      flag_anap = 0;
	    }
	  }
	  else    /* flag_anap */
	  {
	    if(bin_low-bin_high >= MAX_DIF || bin_high <= MIN_VALUE) 
	    {
	      for(l=0; l<el_up; l++)
	        vol_pol[l][i].ray[k]=1;
	      flag_anap = 1;
              stat_anap[i/STEP_STAT_ANAP_AZ][k/STEP_STAT_ANAP_RANGE]++;
	    }
	    else 
	    {
	      for(l=0; l<el_inf; l++)
	        vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
	      flag_anap = 0;
	    }
	  }
	}
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
	}
	else if (bin_low == fondo_scala || bin_high < fondo_scala) /* && bin_high != -20. ) */
	{
	  for(l=0; l<el_inf; l++)
	  {
	    vol_pol[l][i].ray[k]=vol_pol[el_inf][i].ray[k];
	    if(!vol_pol[l][i].flag)
	    {
              vol_pol[l][i].flag = 1;
              vol_pol[l][i].b_header.alfa =(short)(i*.9/FATT_MOLT_AZ);
              vol_pol[l][i].b_header.teta = elev_array[l];
              vol_pol[l][i].b_header.tipo_gran=vol_pol[el_inf][i].b_header.tipo_gran;
              vol_pol[l][i].b_header.max_bin=vol_pol[el_inf][i].b_header.max_bin;
	    }
 	  }
	}
    }			
   
  } 		
  ScrivoStatistica();
  return 0;
}                         


/*===============================================*/
void ScrivoStatistica()
/*===============================================*/
{
  FILE *f_stat;
  int az,ran;
  unsigned char statistica[16][13];
  struct tm *tempo;
  time_t Time;
  char date[20];

  memset(statistica,255,16*13);

  Time = NormalizzoData(old_data_header.norm.maq.acq_date);
  tempo = gmtime(&Time);
  sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
    tempo->tm_mday,tempo->tm_hour, tempo->tm_min);

  for(az=0; az<16; az++)
    for(ran=0; ran<13; ran++)
      if(stat_anap_tot[az][ran] >= N_MIN_BIN)
        statistica[az][ran] = 
             (unsigned char)((stat_anap[az][ran]*100)/stat_anap_tot[az][ran]);
  
  f_stat = fopen(getenv("ANAP_STAT_FILE"),"a");
  fwrite(date,13,1,f_stat);
  fwrite(statistica,16*13,1,f_stat);
  fclose(f_stat);
  return;
}


/*===============================================*/
int write_dbp(const char* nome)
/*===============================================*/
{
  int i,j,size_h;
  FILE *file1;
  
  size_h = sizeof(vol_pol[0][0].b_header);
  file1=fopen(nome,"w");
  if(file1 == NULL) return 1;
  
  if(fwrite(&old_data_header,sizeof(old_data_header),1,file1)== 0) return 2;
  
  for (i=0; i<NEL; i++)
    for (j=0; j<400; j++)
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
	FILE *file1;
	int i,j;
	file1=fopen(getenv("FIRST_LEVEL_FILE"),"r");

	if(file1 == NULL){ 
	  sprintf(errori,"Errore Apertura mappa statica ");
	  ScrivoLog(16,errori);	
return ;
}
	for(i=0; i<400; i++)
	  fread(&first_level[i][0],MAX_BIN,1,file1);
	fclose(file1);
  /*-------------------------------
    patch per espandere il clutter
    -------------------------------*/
#ifdef MEDIUM
  {
    unsigned char first_level_tmp[400][MAX_BIN];
    int k;
    memcpy(first_level_tmp,first_level,sizeof(first_level));
    for (i=400; i<800; i++)
      {
	for (j=0; j<440; j++)
	  {
	    for (k=i-1; k<i+2; k++)
	      if(first_level[i%400][j] < first_level_tmp[k%400][j])
		first_level[i%400][j]=first_level_tmp[k%400][j];
	  }
      }
  }
#endif

}


/*===============================================*/
void creo_cart()
/*===============================================*/
{

	int i,j,quad,x,y,irange,az,iaz,az_min,az_max;
	static int flag = 1;

	if(flag)
	{
	  creo_matrice_conv();
	  flag = 0;
	}

	for(i=0; i<MAX_BIN*2; i++)
	  for(j=0; j<MAX_BIN*2; j++)
	    cart[i][j] = MISSING;


	for(quad=0; quad<4; quad++)
	  for(i=0; i<MAX_BIN; i++)
	    for(j=0; j<MAX_BIN; j++)
	    {
	      irange = (int)range[i][j];
	      if(range[i][j] - irange >= .5) irange++;
	      if(irange < MAX_BIN)
	      {
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
	          az_min = az_min + 400;
	          az_max = az_max + 400;
	        }
	        for(iaz = az_min; iaz<az_max; iaz++)
	          if(cart[x][y]<=vol_pol[0][iaz%400].ray[irange])
	            cart[x][y] = vol_pol[0][iaz%400].ray[irange];
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
void	creo_cart_5x5(int imin,int imax,int jmin,int jmax,int dx,int dy)
/*===============================================*/
{
	int i,j,x,y;
	double sum_z,sum_r,ntot;

	for(i=imin; i<=imax; i++)
	  for(j=jmin; j<=jmax; j++)
	  {
	    sum_r = 0;
	    sum_z = 0;
	    ntot = 0;
	    for(x = 0; x < dx; x++)
	      for(y = 0; y < dy; y++)
	        if(cart[MAX_BIN+i*dx+x][MAX_BIN+j*dy+y] != MISSING)
	        {
#ifdef Z_AVERAGE
	          sum_z = sum_z + Z_func(cart[MAX_BIN+i*dx+x][MAX_BIN+j*dy+y]); 
#endif
#ifdef R_AVERAGE
	          sum_r = sum_r + mp_func(cart[MAX_BIN+i*dx+x][MAX_BIN+j*dy+y]);
#endif
	          ntot++;
	        }
	    if(ntot == 0)
	    {
#ifdef Z_AVERAGE
	      cart_5x5_z[i-imin][j-jmin] = -0.01;
#endif
#ifdef R_AVERAGE
	      cart_5x5_r[i-imin][j-jmin] = -0.01;
#endif
	    }
	    else
	    {
#ifdef Z_AVERAGE
	      cart_5x5_z[i-imin][j-jmin] = Z_to_mp_func(sum_z/ntot);
#endif
#ifdef R_AVERAGE
	      cart_5x5_r[i-imin][j-jmin] = (sum_r/ntot);
#endif
	    }
	  }
}


/*===============================================*/
float Z_func(unsigned char byte)
/*===============================================*/
{
   static char flag=1;
   static float Z[256];
   int i;
   if(flag)
   {
     for(i=0;i<256; i++)
       Z[i]=pow(10.,i*8./255.-2);
     flag = 0;
   }
    return Z[byte];
}


/*===============================================*/
float Z_to_mp_func(double Z)
/*===============================================*/
{
  return (pow(Z/200.,1./1.6));
}

/*===============================================*/
float mp_func(unsigned char byte)
/*===============================================*/
{
   static char flag=1;
   static float MP[256];
   int i;
   if(flag)
   {
     for(i=0;i<256; i++)
       MP[i]=pow(pow(10.,i*8./255.-2)/200.,1./1.6);
     flag = 0;
   }
   return MP[byte];
/*---	return (pow(pow(10.,byte*8./255.-2)/200.,1./1.6)); --*/
}


/*===============================================*/
void scrivo_5x5 (int imin,int imax,int jmin,int jmax,int dx,int dy)
/*===============================================*/
{
	struct tm *tempo;
	time_t time;
	/*----------------------------------------------------------------------------*/
	/*	apertura file dati di output					      */
	/*----------------------------------------------------------------------------*/

	time = NormalizzoData(old_data_header.norm.maq.acq_date);
	printf( "%ld \n", time);
	printf("%d\n",old_data_header.norm.maq.acq_date);

	tempo = gmtime(&time);
#ifdef Z_AVERAGE
	sprintf(nome_file,"%s/%02d%02d%02d%02d%02d.5X5",getenv("OUTPUT_Z_DIR"),
	  tempo->tm_year, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

	printf(" scrivo il file \n%s\n",nome_file);
	output = fopen(nome_file,"w");
	if(output == NULL )
	{
	  sprintf(errori,"Errore Apertura File Dati Output 5X5 %s",nome_file);
	  ScrivoLog(16,errori);	
	  exit(1);
	}

	fprintf(output,"%04d%02d%02d%02d%02d 6\n",
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);
#ifdef MEDIUM
	fprintf(output,"%5d%5d%5d%5d%5d%5d%10.2f%10.2f\n",
		MAX_BIN*2,MAX_BIN*2,imin,imax,jmin,jmax,
		dx*size_cell[old_data_header.norm.maq.resolution],
	  dy*size_cell[old_data_header.norm.maq.resolution]);
#endif
	for (x=imin ; x<=jmax ; x++)
	  {
	    for (y=jmin ; y<=jmax ; y++)
	      fprintf(output,"%6d",(int)(cart_5x5_z[y-jmin][x-imin]*100.));
	    fprintf(output,"\n");
	  }
	fclose(output);
#endif
#ifdef R_AVERAGE
	sprintf(nome_file,"%s/%02d%02d%02d%02d%02d.5X5",getenv("OUTPUT_RAIN_DIR"),
	  tempo->tm_year, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

	output = fopen(nome_file,"w");
	if(output == NULL )
	{
	  sprintf(errori,"Errore Apertura File Dati Output 5X5 %s",nome_file);
	  ScrivoLog(16,errori);	
	  exit(1);
	}

	fprintf(output,"%04d%02d%02d%02d%02d 6\n",
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

	for (x=imin ; x<jmax ; x++)
	  for (y=jmin ; y<jmax ; y++)
	    fprintf(output,"%6d",(int)(cart_5x5_z[y-jmin][x-imin]*100.));

	fclose(output);
#endif

#ifdef SIDME
/*--------------------------
 | scrittura file per sidme |
  --------------------------*/
	name_sidme=fopen(getenv("NOME_SIDME"),"a");
       	if (name_sidme != NULL)
	{
	  fprintf(name_sidme,"%02d%02d%02d%02d%02d.5X5\n",
		  tempo->tm_year, tempo->tm_mon+1, tempo->tm_mday,
		  tempo->tm_hour, tempo->tm_min);
	}
	fclose(name_sidme);

#endif
        return  ;
}                        //end funzione scrivo_5x5 ()


/*===============================================*/
time_t 	NormalizzoData(time_t time)
/*===============================================*/
{
	int itime;

	itime = time/(NMIN*60);
	/*
	  printf(" esco da Normalizzo %d %d %d \n",time,itime,(time - itime*NMIN*60));
	  printf("%s\n",ctime(&time));
	  printf("%d\n",(NMIN-MAX_TIME_DIFF)*60);
	*/
	if(time - itime*NMIN*60 <=MAX_TIME_DIFF*60) return (itime*NMIN*60); //time è in secondi dal 1970
	if(time - itime*NMIN*60 >(NMIN-MAX_TIME_DIFF)*60) return ((itime+1)*NMIN*60);
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
   char date[80],dateold[80];
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
	       printf("antonio");
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
           bacini[ind] +=   mp_func(cart[12+i][12+j]);
 	   if(bacini[ind] <0. )
           {
	     sprintf(errori,"ERRORE GRAVE Bacino con Precipitzaione Negativa");
	     ScrivoLog(16,errori);	
	     sprintf(errori,"Coordinate Punto %4d %4d  Precipitazione %f",i,j,
		mp_func(cart[12+i][12+j]));
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
		 /*-----------------------------------------------
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
}                   //end funzione bacini()



/*===============================================*/
void ScrivoBacini(const char* nome,char* data, float* array_bacini, int quanti)
/*===============================================*/
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
void AggiornoHistoryFile(char* data, int flag)
/*===============================================*/
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
    //        printf("%s\n",sys_errlist[errno]); 
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
}                   //end funzione AggiornoHistoryFile(data,flag)



/*===============================================*/
void	creo_cart_z_lowris()
/*===============================================*/
{
	int i,j,x,y;
	unsigned char z;

        memset(z_out,0,sizeof(z_out));
	for(i=0; i<CART_DIM_ZLR; i++)
	  for(j=0; j<CART_DIM_ZLR; j++)
	  {
	    z = 0;
	    for(x = 0; x < ZLR_N_ELEMENTARY_PIXEL; x++)
	      for(y = 0; y < ZLR_N_ELEMENTARY_PIXEL; y++)
	        if(cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] != MISSING)
	          if(cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] > z) 
		    z= cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
	    z_out[i][j]=z;
	  }
}


/*===============================================*/
void scrivo_z_lowris ()
/*===============================================*/
{

	char nome_file [150];
	FILE *output;
	struct tm *tempo;
	time_t time;
	/*----------------------------------------------------------------------------*/
	/*	apertura file dati di output					      */
	/*----------------------------------------------------------------------------*/

#ifdef SHORT
	time = NormalizzoData(old_data_header.norm.maq.acq_date);
	tempo = gmtime(&time);
#endif
#ifdef MEDIUM
	//	tempo = gmtime(&old_data_header.norm.maq.acq_date);
	time = NormalizzoData(old_data_header.norm.maq.acq_date);
	tempo = gmtime(&time);
#endif
	sprintf(nome_file,"%s/%04d%02d%02d%02d%02d.ZLR",
	  getenv("OUTPUT_Z_LOWRIS_DIR"),
	  tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
	  tempo->tm_hour, tempo->tm_min);

	output = fopen(nome_file,"w");
	if(output == NULL )
	{
	  sprintf(errori,"Errore Apertura File Dati Output 1X1 %s",nome_file);	
          ScrivoLog(16,errori);	
	  exit(1);
	}

	fwrite(z_out,sizeof(z_out),1,output);
	fclose(output);
        return  ;
}


/*===============================================*/
void ScrivoLog(int i, char* stringa)
/*===============================================*/
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



	     fprintf(log,"\n");
	     fprintf(log,"-----------------------------------------------------------------\n");
	     fprintf(log,"Variabili d'Ambiente\n");
	     fprintf(log,"LISTA_FILE = %s\n",getenv("LISTA_FILE"));
	     fprintf(log,"LAST_FILE = %s\n",getenv("LAST_FILE"));
	     fprintf(log,"ANAP_STAT_FILE = %s\n",getenv("ANAP_STAT_FILE"));
	     fprintf(log,"FIRST_LEVEL_FILE = %s\n",getenv("FIRST_LEVEL_FILE"));
	     fprintf(log,"OUTPUT_Z_DIR = %s\n",getenv("OUTPUT_Z_DIR"));
	     fprintf(log,"OUTPUT_RAIN_DIR = %s\n",getenv("OUTPUT_RAIN_DIR"));
	     fprintf(log,"MATRICE_BACINI = %s\n",getenv("MATRICE_BACINI"));
	     fprintf(log,"BACINI_AUS_FILE = %s\n",getenv("BACINI_AUS_FILE"));
	     fprintf(log,"BACINI_DIR = %s\n",getenv("BACINI_DIR"));
	     fprintf(log,"BACINI_HISTORY_FILE = %s\n",getenv("BACINI_HISTORY_FILE"));
	     fprintf(log,"OUTPUT_Z_LOWRIS_DIR = %s\n",getenv("OUTPUT_Z_LOWRIS_DIR"));
	     fprintf(log,"NOME_SIDME = %s\n",getenv("NOME_SIDME"));
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
  }

}                        //end funzione ScrivoLog(i,stringa)


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
void   write_xdr(float* bacini, time_t Time)
/*===============================================*/
{
  FILE *f_output;
  XDR xdr;

  f_output=fopen(getenv("BACINI_BOLOGNA"),"a");
  xdrstdio_create(&xdr ,f_output,XDR_ENCODE);
  xdr_long(&xdr,&Time);
  //  if(xdr_vector(&xdr,bacini,500,4,xdr_float) == 1)
  xdr_destroy(&xdr);
  fclose(f_output);
  return;
}

