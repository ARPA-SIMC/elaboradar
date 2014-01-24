
/*----------------------------------------------------------------------------*/
  /*    INCLUDE file                                      */
/*----------------------------------------------------------------------------*/

#include "cum_bac.h"
#include "logging.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <func_Z_R.h>
#ifdef __cplusplus
}
#endif

#include <Q_components.h>


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




#ifdef QUALITY
#include <qual_par.h>
#endif
#ifdef CLASS
#include <par_class.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
// nr
#include <nrutil.h>
#include <nr.h>
#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------------*/
/*      FINE SEZIONE    INCLUDE                   */
/*----------------------------------------------------------------------------*/
// DEFINIZIONI PREPROCESSORE NB: VPR E CLASS IMPLICANO QUALITY ,

//Definizioni geometriche
#define AMPLITUDE 0.9 /* esternalizzo?*/ // ampiezza fascio radar
#define DTOR  M_PI/180. /* esternalizzo?*/ //fattore conversione gradi-radianti
#define CONV_RAD 360./4096.*DTOR  // fattore conversione unità angolare radar-radianti

//Dimensioni massime range polare
#define MAX_BIN 512
#define MAX_DIM 512

// Fattori moltiplicativi per passare da unità azimut / unità elevazione a gradi
#define FATT_MOLT_EL ((double) 360./(double)4096.)
#define FATT_MOLT_AZ ((double) 360./(double)4096.)

// parametri ereditati da programma beam blocking:numero elevazioni da programma beam blocking ; le matrici ivi definite considerano questo
#define NSCAN 6

/*----------------------------------------------------------------------------*/
/*  DICHIARATIVE GLOBALI                              */
/*----------------------------------------------------------------------------*/

// Soglie algoritmi
#define OVERBLOCKING 51 /* minimo BB non accettato*/
#define THRES_ATT 0 /* minimo valore di Z in dBZ per calcolare att rate */
#define SOGLIA_TOP 20 // soglia per trovare top
#define MISSING 0 /*valore mancante*/

/*extern char *sys_errlist[];  ANNA 17-2-2006 sostituito con strerror che sta in string.h */
extern int errno;

/* comstart lineargauss
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


#ifdef QUALITY
/*
   comstart caratterizzo_volume
   idx calcola qualita' volume polare
   calcola qualita' volume polare
   NB il calcolo è fatto considerando q=0 al di sotto della mappa dinamica.
   per ora drrs=dist nche nel caso di Gattatico, mentre dtrs è letto da file
   si puo' scegliere tra qualita' rispetto a Z e rispetto a R, in realtà per ora sono uguali.

   double PIA;  path integrated attenuation
   float dh: dimensione verticale bin calcolata tramite approcio geo-ottico
   float drrs: distanza radiosondaggio,
   float dtrs: tempo dal radiosondaggio
   float dist: distanza dal radar
   float quo:  quota bin
   float el:   elevazione
   float rst:  raggio equivalente in condizioni standard
   float dZ:   correzione vpr
   float sdevZ:  stand. dev. correzione vpr
   unsigned char bb: beam blocking
   unsigned char cl: indice clutter da anaprop

//--- nb. non ho il valore di bb sotto_bb_first_level
comend
*/
void CUM_BAC::caratterizzo_volume()
{
    int i,l,k;
    double PIA;
    float dh=1.,dhst=1.,drrs=1.,dist=1.,elevaz;
    unsigned char bb=0,cl=0 ;

    //----------ciclo su NSCAN(=6), cioè sul numero di elevazioni (nominali) per le quali ho calcolato il beam blocking
    /* a questo punto servono: bb, cl,  PIA, dtrs e drrs radiosond, quota, hsup e hinf beam-----------------*/

    //for (l=0; l<NSCAN; l++)/*ciclo elevazioni*/// NSCAN(=6) questo lascia molti dubbi sul fatto che il profilo verticale alle acquisizioni 48, 19 etc..  sia realmente con tutti i dati! DEVO SOSTITUIRE CON nel E FARE CHECK.

    for (l=0; l<NEL; l++)/*ciclo elevazioni*/// VERIFICARE CHE VADA TUTTO OK
    {
        for (i=0; i<NUM_AZ_X_PPI; i++)/*ciclo azimuth*/
        {
            //-----elevazione reale letta da file* fattore di conversione 360/4096
            elevaz=(float)(vol_pol[l][i].teta_true)*CONV_RAD;//--- elev reale

            //--assegno PIA=0 lungo il raggio NB: il ciclo nn va cambiato in ordine di indici!
            PIA=0.;

            for (k=0; k<vol_pol[0][i].b_header.max_bin; k++)/*ciclo range*/
            {
                //---------distanza in m dal radar (250*k+125 x il corto..)
                dist= k*size_cell[old_data_header.norm.maq.resolution]+size_cell[old_data_header.norm.maq.resolution]/2.;/*distanza radar */

                //-----distanza dal radiosondaggio (per GAT si finge che sia colocato ..), perchè? (verificare che serva )
                drrs=dist;
                /* if (!(strcmp(sito,"GAT")) ) {  */
                /*     drrs=dist; */
                /* } */
                /* if (!(strcmp(sito,"SPC")) ) {  */
                /*     drrs=dist; */
                /* } */


                //assegno la PIA (path integrated attenuation) nel punto e POI la incremento  (è funzione dell'attenuazione precedente e del valore nel punto)
                if (l == elev_fin[i][k]) att_cart[i][k]=DBtoBYTE(PIA);
                PIA=attenuation(vol_pol[l][i].ray[k],PIA);

                //------calcolo il dhst ciè l'altezza dal bin in condizioni standard utilizzando la funzione quota_f e le elevazioni reali
                dhst =quota_f(elevaz+0.45*DTOR,k)-quota_f(elevaz-0.45*DTOR,k);


                //----qui si fa un po' di mischione: finchè ho il dato dal programma di beam blocking uso il dh con propagazione da radiosondaggio, alle elevazioni superiori assegno dh=dhst  e calcolo quota come se fosse prop. standard, però uso le elevazioni nominali

                if (l<NSCAN-1   ) {
                    dh=hray_inf[k][l+1]-hray_inf[k][l]; /* differenza tra limite sup e inf lobo centrale secondo appoccio geo-ott*/
                }
                else {
                    dh=dhst; /* non ho le altezze oltre nscan-1 pero' suppongo che a tali elevazioni la prop. si possa considerare standard*/
                    hray[k][l]=quota_f(elevaz,k);//non lo assegno

                }

                if (l-elev_fin[i][k] <0) {
                    cl=ANAP_YES;
                    bb=BBMAX;
                }
                if (l == elev_fin[i][k] ) {
                    cl=dato_corrotto[i][k];  /*cl al livello della mappa dinamica*/
                    bb=beam_blocking[i][k];  /*bb al livello della mappa dinamica *///sarebbe da ricontrollare perchè con la copia sopra non è più così
                }
                if (l-elev_fin[i][k] >0 ) {
                    cl=0;       /*per come viene scelta la mappa dinamica si suppone che al livello superiore cl=0 e bb=0*/
                    bb=0;   // sarebbe if (l-bb_first_level[i][k] >0  bb=0;  sopra all'elevazione per cui bb<soglia il bb sia =0 dato che sono contigue o più però condiz. inclusa
                }

                //------dato che non ho il valore di beam blocking sotto i livelli che ricevo in ingresso ada progrmma beam blocking e
                //--------dato che sotto elev_fin rimuovo i dati come fosse anaprop ( in realtà c'è da considerare che qui ho pure bb>50%)
                //--------------assegno qualità zero sotto il livello di elev_fin (si può discutere...), potrei usare first_level_static confrontare e in caso sia sotto porre cl=1
                if (l-elev_fin[i][k] <0) {
                    qual[l][i][k]=0;
                    cl=2;
                }
                //--------bisogna ragionare di nuovo su definizione di qualità con clutter se si copia il dato sopra.----------
                else {

                    //-----------calcolo la qualità----------
                    qual[l][i][k]=(unsigned char)(func_q_Z(cl,bb,dist,drrs,dtrs,dh,dhst,PIA)*100);
                }

                if(qual[l][i][k]==0) qual[l][i][k]=1;//????a che serve???
#ifdef VPR
                /* sezione PREPARAZIONE DATI VPR*/
                if(cl==0 && bb<BBMAX_VPR )   /*pongo le condizioni per individuare l'area visibile per calcolo VPR, riduco il bb ammesso (BBMAX_VPR=20)*/ //riveder.....?????
                    flag_vpr[l][i][k]=1;
#endif
                //------------trovo il top per soglia
                if (BYTEtoDB(vol_pol[l][i].ray[k]) > SOGLIA_TOP )
                    top[i][k]=(unsigned char)((quota_f(elevaz,k))/100.); //top in ettometri

            }

        }
    }

    return;
}
#endif


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
        att_tot=att_tot+2.*att_rate*0.001*size_cell[old_data_header.norm.maq.resolution];
        if (att_tot>BYTEtoDB(254)) att_tot=BYTEtoDB(254);
    }
    return att_tot;
}


void CUM_BAC::creo_cart()
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
                if(irange < MAX_BIN)        {
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
                            topxy[x][y]=top[iaz%NUM_AZ_X_PPI][irange];
#ifdef QUALITY
                            qual_Z_cart[x][y]=qual[elev_fin[iaz%NUM_AZ_X_PPI][irange]][iaz%NUM_AZ_X_PPI][irange];
                            quota_cart[x][y]=quota[iaz%NUM_AZ_X_PPI][irange];
                            dato_corr_xy[x][y]=dato_corrotto[iaz%NUM_AZ_X_PPI][irange];
                            beam_blocking_xy[x][y]=beam_blocking[iaz%NUM_AZ_X_PPI][irange];
                            elev_fin_xy[x][y]=elev_fin[iaz%NUM_AZ_X_PPI][irange];
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

void CUM_BAC::creo_matrice_conv()
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

void CUM_BAC::creo_cart_z_lowris()
{
    int i,j,x,y,cont;
    unsigned char z,q,nv,c1x1,traw,dc1x1,el1x1,bl1x1,Zr1x1;
    unsigned short q1x1;
    float zm;

    //tolta qui inizializzazione di z_out che era duplicata (già fatta all'inizio del main)
    // ciclo sui punti della nuova matrice. per il primo prenderò il massimo tra i primi sedici etc..
    for(i=0; i<CART_DIM_ZLR; i++)
        for(j=0; j<CART_DIM_ZLR; j++)
        {
            //reinizializzo tutte le variabili calcolate dentro la funzione .
            z = 0;
            q = 0;
            zm = 0.;
            dc1x1=0;
            el1x1=0;
            q1x1=0;
            c1x1=0;
            bl1x1=0;
            cont=0;
            traw=0;
            for(x = 0; x < ZLR_N_ELEMENTARY_PIXEL; x++)
                for(y = 0; y < ZLR_N_ELEMENTARY_PIXEL; y++)
                    //ciclo a passi di 4 in x e y nella matrice a massima risoluzione, cercando il valore massimo di z tra i primi sedici e attribuendolo al primo punto della matrice a bassa risoluzione e poi i tra i secondi sedici e attribuendolo al secondo punto etc...
                {
                    if(cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] != MISSING)
                        if(cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET] > z){
                            z= cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
                            traw=topxy[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
#ifdef QUALITY
                            q=qual_Z_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
                            q1x1=quota_cart[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
                            dc1x1=dato_corr_xy[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
                            el1x1=elev_fin_xy[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];
                            bl1x1=beam_blocking_xy[i*ZLR_N_ELEMENTARY_PIXEL+x+ZLR_OFFSET][j*ZLR_N_ELEMENTARY_PIXEL+y+ZLR_OFFSET];

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
#ifdef QUALITY
                    qual_Z_1x1[i][j]=q;
                    quota_1x1[i][j]=128+(unsigned char)(q1x1/100);
                    dato_corr_1x1[i][j]=dc1x1;
                    elev_fin_1x1[i][j]=el1x1;
                    beam_blocking_1x1[i][j]=bl1x1;
#endif
                    top_1x1[i][j]=traw;

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

void CUM_BAC::scrivo_out_file_bin (const char *ext,char *content,char *dir,size_t size, void  *matrice)
{
    char nome_file [512];
    FILE *output;
    struct tm *tempo;
    time_t time;
    /*----------------------------------------------------------------------------*/
    /*    apertura file dati di output                                          */
    /*----------------------------------------------------------------------------*/
    //definisco stringa data in modo predefinito
#ifdef SHORT
    time = NormalizzoData(old_data_header.norm.maq.acq_date);
    tempo = gmtime(&time);
#endif
#ifdef MEDIUM
    tempo = gmtime(&old_data_header.norm.maq.acq_date);
    time = NormalizzoData(old_data_header.norm.maq.acq_date);
    tempo = gmtime(&time);
#endif
    snprintf(nome_file, 512, "%s/%04d%02d%02d%02d%02d%s",dir,
            tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
            tempo->tm_hour, tempo->tm_min, ext);

    output=controllo_apertura(nome_file,content,"w");
    printf("aperto file %s dimensione matrice %zd\n",nome_file,size);

    fwrite(matrice,size,1,output);
    fclose(output);
    return;
}

/*=======================================================================================*/



#ifdef VPR
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
cappi1100[1200*MAX_BIN]: vettore volumi buoni tra 1000 e 1200 m
long int area_vpr[NMAXLAYER]; area totale usata per calcolo vpr

*/

int CUM_BAC::func_vpr(long int *cv, long int *ct, float vpr1[], long int area_vpr[], const char *sito)
{
    int l,i,iA,k,ilay,il,ilast,iaz_min,iaz_max,icounter,naz,nra;
    int strat;
    long int dist,dist_plain,vert_ext,vol_rain,somma;
    float noval,area,elevaz,quota_true_st;
    float grad;

    //----------------inizializzazioni----------------
    somma=0;
    stdev = -1.;
    noval=NODATAVPR;
    icounter=0;
    vert_ext=0;

    *cv=0;
    *ct=0;


    /*SECTION A: cv and ct retrieval and calculation of current area- weighted vpr*/


    //------------riconoscimento sito per definizione limiti azimut---------

    if (!(strcmp(sito,"SPC"))){  /* limiti azimuth per calcolo VPR sarebbe bello fare un CASE ma con stringhe non si puo'!*/
        iaz_min=IAZ_MIN_SPC;
        iaz_max=IAZ_MAX_SPC;
    }
    if (!(strcmp(sito,"GAT"))){
        iaz_min=IAZ_MIN_GAT;
        iaz_max=IAZ_MAX_GAT;
    }
    if (strcmp(sito,"SPC") && strcmp(sito,"GAT")) {
        fprintf(log_vpr,"errore, sito sconosciuto o non definito");
        exit (1);
    }

    //

    for (int l=0; l<NEL; l++)//ciclo elevazioni
    {
        for (int k=0; k<MAX_BIN; k++)/*ciclo range*/
        {
            //-------------calcolo distanza-----------
            dist=k*(long int)(size_cell[old_data_header.norm.maq.resolution])+(int)(size_cell[old_data_header.norm.maq.resolution])/2.;

            //-----ciclo settore azimut???
            for (iA=iaz_min; iA<iaz_max; iA++)//ciclo sulle unità di azimut  (0,9°)    ---------
            {
                //----------ottengo il valore dell'indice in unità di azimut (0,9°) ---------
                i=(iA+NUM_AZ_X_PPI)%NUM_AZ_X_PPI;

                //--------calcolo elevazione e quota---------
                elevaz=(float)(vol_pol[l][i].teta_true)*CONV_RAD;
                quota_true_st=quota_f(elevaz,k);

                //--------trovo ilay---------
                ilay=floor(quota_true_st/TCK_VPR);//  in teoria quota indipendente da azimuth  , in realtà no (se c'è vento)

                if (ilay <0 || ilay >= NMAXLAYER) {
                    //fprintf(log_vpr,"ilay %d errore\n",ilay);
                    break;
                }


                vol_rain=0;
                // dist=(long int)(dist*cos((float)(vol_pol[elev_fin[i][k]][i].teta_true)*CONV_RAD));

                /* //---------calcolo la distanza proiettata sul piano-------------  */

                dist_plain=(long int)(dist*cos(elevaz));
                if (dist_plain <RMIN_VPR || dist_plain > RMAX_VPR )
                    flag_vpr[l][i][k]=0;

                if(qual[l][i][k]<QMIN_VPR ) flag_vpr[l][i][k]=0;

                //AGGIUNTA PER CLASS
# ifdef CLASS
                if(conv[i][k]>= CONV_VAL){
                    flag_vpr[l][i][k]=0;
                }
#endif

                // ------per calcolare l'area del pixel lo considero un rettangolo dim bin x ampiezzamediafascio x flag vpr/1000 per evitare problemi di memoria?
                area=size_cell[old_data_header.norm.maq.resolution]*dist_plain*AMPLITUDE*DTOR*flag_vpr[l][i][k]/1000.; // divido per  mille per evitare nr troppo esagerato

                // ------incremento il volume totale di area

                *cv=*cv+(long int)(area);


                //---------------------condizione per incrementare VPR contributo: valore sopra 13dbz, qualità sopra 20 flag>0 (no clutter e dentro settore)------------------
                if (BYTEtoDB(vol_pol[l][i].ray[k])> THR_VPR &&  flag_vpr[l][i][k]>0 )
                {
                    //-------incremento il volume di pioggia = pioggia x area
                    vol_rain=(long int)(BYTE_to_mp_func(vol_pol[l][i].ray[k],aMP,bMP)*area);//peso ogni cella con la sua area

                    //-------incremento l'area precipitante totale ct,aggiungendo però,cosa che avevo messo male una THR solo per ct, cioè per il peso
                    if (BYTEtoDB(vol_pol[l][i].ray[k])> THR_PDF)
                        *ct=*ct+(long int)(area);

                    //------se l'area in quello strato è già maggiore di 0 allora incremento il volume dello strato altrimenti lo scrivo ex novo. poi vpr1 andrà diviso per l'area
                    if (area_vpr[ilay]> 0) vpr1[ilay]=vpr1[ilay]+(float)(vol_rain);
                    else vpr1[ilay]=(float)(vol_rain);

                    //------incremento l'area dello strato----------
                    area_vpr[ilay]=area_vpr[ilay]+area;

                    /* //------cappi 1100 m */
                    /* if (abs(quota_true_st-1100)<TCK_VPR/2) { */
                    /*   cappi1100[icounter]=vol_rain;     */
                    /*   icounter=icounter+1; */
                    /* }   */
                }
            }
        }
    }

    fprintf(log_vpr,"calcolati ct e cv ct= %li cv= %li\n",*ct,*cv);

    /*SECTION B: vpr quality checks and re-normalisation of vpr*/

    //--------------CONTROLLO DI QUALITA' E NORMALIZZAZIONE DEL PROFILO ISTANTANEO CALCOLATO
    //-------- se il volume supera quello minimo------
    if ((*ct) > CT_MIN*(*cv)) {

        ilast=0;
        vert_ext=0;


        //----- calcolo 'estensione verticale del profilo , negli strati dove l'area è troppo piccola assegno NODATAVPR,  NORMALIZZO il profilo, e se l'estensione è minore di VEXTMIN_VPR esco-


        for (ilay=0; ilay<NMAXLAYER; ilay++){


            fprintf(log_vpr,"  ilay %d area_vpr= %ld  ct= %ld  cv= %ld \n", ilay, area_vpr[ilay],*ct,*cv );

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
                    fprintf(log_vpr,"raggiunta cima profilo \n");
                    ilast=ilay-1;// c'era errore!!!

                    //---------- raggiunta la cima profilo faccio check immediato sull'estensione verticale
                    if (vert_ext<VEXTMIN_VPR ){
                        fprintf(log_vpr,"estensione profilo verticale troppo bassa\n");
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
        fprintf(log_vpr,"volume precipitante troppo piccolo\n");
        *ct=0;
        ilast=0;
        for  (il=0; il<NMAXLAYER; il++) vpr1[il]=NODATAVPR; //--devo riassegnare o mi rimane 'sporco' forse si potrebbe usare una ver diversa
        return(1);
    }


    //------calcolo il gradiente del profilo come media del gradiente negli ultimi 3 strati per assegnare la parte 'alta' (novità)

    grad=((vpr1[ilast]-vpr1[ilast-1]) + (vpr1[ilast-1]-vpr1[ilast-2])/(2.)+ (vpr1[ilast-2]-vpr1[ilast-3])/(3.) ) /3.;
    if (grad > 0.002)
        grad=-0.03 ;
    fprintf(log_vpr," %f \n", grad);

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

float CUM_BAC::comp_levels(float v0, float v1, float nodata, float peso)
{
    float result;
    /* if ((v0<nodata+1)&&(v1<nodata+1)) result=nodata; */
    /* if (v0<nodata+1) result=v1; */
    /* if (v1<nodata+1) result=v0;         */
    if ((v0>nodata) && (v1>nodata)  ) result=((1.-peso)*v0+peso*v1); /* in questa configurazione il vpr è di altezza costante  nel tempo ma un po' 'sconnesso' in alto*/

    else result=nodata;
    return(result);
}

long int CUM_BAC::profile_gap(char nomefile[])
{
    FILE *file;
    time_t last_time;
    long int gap1;

    gap1=(long int)(100 );
    file = fopen(nomefile,"r");
    if(file != NULL ){ /*contemplo la prima iterazione dopo installazione*/
        fread(&last_time,4,1,file);
        fclose(file);
        gap1=abs(old_data_header.norm.maq.acq_date-last_time)/900;
        fprintf (log_vpr,"old_data_header.norm.maq.acq_date last_time gap %d %ld %ld \n",old_data_header.norm.maq.acq_date,last_time,gap1);
    }

    return(gap1);
}


/*===============================================*/

//----------ALGORITMO
/*  combina il profilo verticale corrente con quello precedente tramite il metodo di Germann (2003)
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
*/


int CUM_BAC::combina_profili(const char *sito)
{
    long int c0,*cv,*ct;
    float vpr0[NMAXLAYER],vpr1[NMAXLAYER],vpr_dbz;
    float alfat,noval;
    int mode,ilay,gap_res,heating_res,i,foundlivmin=0,il,ier_ap,ier_cp,combinante=0; // combinante: variabile che contiene presenza vpr alternativo
    long int  area[NMAXLAYER],ar=0;
    int n=0,diff=0;
    char nomefile[150],stringa[100];
    struct tm *T_tempo;
    time_t Time,T_Time;
    FILE *file;

    mode=MOD_VPR;

    for (i=0;i<NMAXLAYER;i++) {
        area[i]=0;//inizializzo area
        area_vpr[i]=0;
        vpr0[i]=NODATAVPR;
        vpr1[i]=NODATAVPR;
    }
    noval=NODATAVPR;


    /* questo per fare ciclo sul vpr vecchio*/
    Time = NormalizzoData(old_data_header.norm.maq.acq_date);

    //--------inizializzo cv e ct-------------//
    //-----calcolo del profilo istantaneo:faccio func_vpr-----//

    cv=( long int *)malloc(sizeof(long int));
    ct=( long int *)malloc(sizeof(long int));
    if (ct == NULL) printf("malloc fallita per ct\n");
    if (cv == NULL) printf("malloc fallita per cv\n");
    *cv=0;
    *ct=0;

    ier_vpr=func_vpr(cv,ct,vpr1,area_vpr,sito); // ho fatto func_vpr, il profilo istantaneo
    fprintf(log_vpr,"fatta func vpr\n");


    /*modalità VPR combinato*/

    if(mode == 0) {

        /*----calcolo il peso c0 per la combinazione dei profili*/

        c0=2*(*cv);

        /*------calcolo la distanza temporale che separa l'ultimo profilo calcolato dall'istante attuale--*/
        /* (dentro il file LAST_VPR c'è una data che contiene la data cui si riferisce il vpr in n0 di secondi dall'istante di riferimento)*/

        gap=profile_gap(getenv("LAST_VPR"));


        /*------leggo il profilo vecchio più recente di MEMORY ----*/
        /*------nota bene: è in R ovvero  pioggia!! ----*/

        file=fopen(getenv("VPR0_FILE"),"r");
        if(file == NULL ) {
            fprintf(log_vpr,"non esiste file vpr vecchio \n",getenv("VPR0_FILE"));

            //----se file non esiste assegno gap=100----
            gap=100;
        }

        //------------se gap < MEMORY leggo vpr e area per ogni strato-----------
        //--------qui dentro c'è la funzione controllo_apertura, per la quale rimandiamo a dopo qualsiasi commento--------

        if (gap<=MEMORY){
            combinante=1;
            controllo_apertura(getenv("VPR0_FILE")," old VPR ","r");
            for(ilay=0; ilay<NMAXLAYER; ilay++){

                //-----leggo vpr e area per ogni strato----
                ier_cp=fscanf(file,"%f %li\n",&vpr0[ilay],&area[ilay]);
            }
            fprintf(log_vpr,"fatta lettura vpr\n");
            fclose(file);

        }

        //-----Se gap > MEMORY

        //a)----- tento .. sono in POST-ELABORAZIONE:----

        //-----devo andare a ricercare tra i profili 'buoni' in archivio quello con cui combinare il dato----
        //---- trattandosi di profili con data nel nome del file, costruisco il nome a partire dall'istante corrente ciclando su un numero di quarti d'ora
        //---- pari a memory finchè non trovo un profilo. se non lo trovo gap resta=100

        else {
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
                        ier_cp=fscanf(file," %i %f %li", &il, &vpr0[ilay], &ar);  //---NB il file in archivio è in dBZ e contiene anche la quota----

                        //---- converto in R il profilo vecchio--
                        if (vpr0[ilay]>0){
                            vpr_dbz=vpr0[ilay];
                            vpr0[ilay] = DBZtoR(vpr_dbz,aMP,bMP);
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
        fprintf(log_vpr,"gap %li \n",gap);

        //TOLTA: combinazione dell'istantaneo col vecchio dell'altro radar purchè sia 'caldo' (non prevista la post-combinazione)


        //-----se è andata male la ricerca dell'altro e anche il calcolo dell'istantaneo esco

        if ( !combinante && ier_vpr){

            free(cv);
            free(ct);
            return (1);
        }


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
            alfat=(float)(*ct)/(c0+(*ct));
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
        if (ier_vpr) {
            free(cv);
            free(ct);
            return (1);
        }
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
    fprintf(log_vpr," livmin %i  \n", livmin);

    if (livmin>=(NMAXLAYER-1)*TCK_VPR+TCK_VPR/2  || !foundlivmin) return (1);



    //-----scrivo il profilo e la sua area-----
    file=controllo_apertura(getenv("VPR0_FILE")," ultimo vpr ","w");
    for (ilay=0;  ilay<NMAXLAYER; ilay++)
        fprintf(file," %10.3f %li\n",vpr[ilay], area_vpr[ilay]);
    fclose(file);


    //-----libero memoria-----
    free(cv);
    free(ct);

    return(0);
}

int CUM_BAC::profile_heating()
#include <vpr_par.h>
{
    int heating;
    FILE *file;

    //---leggo ultimo file contenente riscaldamento , se non esiste impongo heating=0 (verificare comando)

    if (!access(getenv("VPR_HEATING"),F_OK)){
        file=fopen(getenv("VPR_HEATING"),"r");
        if (file == NULL )  {
            fprintf (log_vpr,"non ho il file di riscaldamento %s \n",getenv("VPR_HEATING"));
            heating=0;
        }
        else{
            fscanf(file,"%i ",&heating);
            fclose(file);
        }
    } /*contemplo la prima iterazione dopo installazione testando l'esistenza del file*/
    else heating=0;


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

        file = controllo_apertura(getenv("LAST_VPR")," ultimo vpr ","w");
        fwrite(&old_data_header.norm.maq.acq_date,4,1,file);
        fclose(file);
    }
    if (heating<0) heating=0;

    //----stampo heating su file
    file=controllo_apertura(getenv("VPR_HEATING")," riscaldamento vpr ","w");
    fprintf(file," %i \n",heating);
    fclose(file);

    //----stampo log vpr
    fprintf (log_vpr,"gap %li heating %i \n",gap,heating);

    return(heating);
}


int CUM_BAC::stampa_vpr()
{
    float vpr_dbz;
    int ilay;
    FILE *file;

    file=controllo_apertura(getenv("VPR_ARCH")," ultimo vpr in dBZ per il plot","w");
    fprintf(file," QUOTA   DBZ    AREA PRECI(KM^2/1000)\n" );
    for (ilay=0;  ilay<NMAXLAYER; ilay++){
        if (vpr[ilay]> 0.001 ) {
            vpr_dbz=RtoDBZ(vpr[ilay],aMP,bMP);
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
int CUM_BAC::corr_vpr(const char *sito)
    //* ====correzione profilo====================================*/

#include <vpr_par.h>

{


    int ilray,ilref,ilay2,i,k,ier_ana,snow,strat;
    float corr,vpr_liq,vpr_hray,hbin,hliq,elevaz,db1,db2;

    /*inizializzazione variabili */
    snow=0;
    //vpr al livello liquido liv liquido e liv max
    vpr_liq=NODATAVPR;
    hliq=NODATAVPR;
    hvprmax=INODATA;

    // analisi vpr

    ier_max=trovo_hvprmax(&hvprmax);
    ier_ana=analyse_VPR(&vpr_liq,&snow,&hliq,sito);
    fprintf (log_vpr,"ier_analisi %i \n",ier_ana) ;

    /* se analisi dice che non è il caso di correggere non correggo (NB in questo caso non riempio la matrice di neve)*/
    if (ier_ana) return 1;

    fprintf (log_vpr,"altezza bright band %i \n",hvprmax) ;
    fprintf (log_vpr,"CORREGGO VPR \n") ;



    //correzione vpr
    for (i=0; i<NUM_AZ_X_PPI; i++){
        for (k=0; k<vol_pol[0][i].b_header.max_bin; k++){
            corr=0.;
            /* trovo elevazione reale e quota bin*/
            //elevaz=(float)(vol_pol[elev_fin[i][k]][i].teta_true)*CONV_RAD;
            hbin=(float)quota[i][k];

            /* se dall'analisi risulta che nevica assegno neve ovunque*/
            if (snow) neve[i][k]=1;
            strat=1;
#ifdef CLASS
            if (conv[i][k] >= CONV_VAL){
                strat=0;
            }
#endif
            //--- impongo una soglia per la correzione pari a 0 dBZ

            if (BYTEtoDB(vol_pol[0][i].ray[k])>THR_CORR && hbin > hliq  && strat){

                //---trovo lo strato del pixel, se maggiore o uguale a NMAXLAYER lo retrocedo di 2, se minore di livmn lo pongo uguale a livmin
                ilray=(hbin>=livmin)?(floor(hbin/TCK_VPR)):(floor(livmin/TCK_VPR));//discutibile :livello del fascio se minore di livmin posto=livmin
                if (ilray>= NMAXLAYER ) ilray=NMAXLAYER-2;//livello del fascio se >= NMAXLAYER posto =NMAXLAYER-2

                //---trovo ilay2 strato con cui mediare per calcolare il vpr a una quota intermedia tra 2 livelli, se l'altezza del bin è sopra metà strato prendo quello sopra altrimenti quello sotto
                if ((int)hbin%TCK_VPR > TCK_VPR/2) ilay2=ilray+1;
                else ilay2=(int)(fabs(ilray-1));
                if (ilay2< livmin) ilay2=livmin;

                //trovo ilref: livello di riferimento per ricostruire il valore vpr al suolo nel caso di neve.
                // in caso di profilo di pioggia mi riporto sempre al valore del livello liquido e questo può essere un punto critico.. vedere come modificarlo.

                ilref=(dem[i][k]>livmin)?(floor(dem[i][k]/TCK_VPR)):(floor(livmin/TCK_VPR));//livello di riferimento; se livello dem>livmin = livello dem altrimenti livmin


                if (vpr[ilref] > 0 && vpr[ilray] > 0 ){ /*devo avere dati validi nel VPR alle quote considerate!*/
                    //-- calcolo il valore del profilo alla quota di interesse
                    vpr_hray=vpr[ilray]+((vpr[ilray]-vpr[ilay2])/(ilray*TCK_VPR-TCK_VPR/2-ilay2*TCK_VPR))*(hbin-ilray*TCK_VPR-TCK_VPR/2); /*per rendere la correzione continua non a gradini */
                    //--identifico le aree dove nevica stando alla quota teorica dello zero termico

                    if (dem[i][k]> hvprmax+HALF_BB-TCK_VPR || snow){ /*classifico neve*/
                        neve[i][k]=1;

                    }

                    //--se nevica la correzione consiste solo nel riportare il valore del vpr al suolo: PROPOSTA: qui si potrebbe generare una mappa di intensità di neve ma deve essere rivisto tutto


                    //if(snow) //A rimosso, faccio una cosa diversa
                    if(neve[i][k]){

                        //faccio la regressione lineare dei punti del profilo sopra il punto del dem
                        //calcolo il valore al livello del dem e lo sostituisco a vpr[ilref] nella correzione
                        // faccio linearizzazione in maniera becera:
                        //vpr[ilref]=(vpr[ilref+7]-vpr[ilref+2])/(5)*(ilref-(ilref+2))+vpr[ilref+2];

                        //passaggio=BYTEtoR(vol_pol,aMP_SNOW,bMP_SNOW)

                        //volpol[0][i].ray[k]=RtoBYTE(passaggio)

                        corr=RtoDBZ(vpr[ilref],aMP,bMP)-RtoDBZ(vpr_hray,aMP,bMP);

                        vol_pol[0][i].ray[k]=DBtoBYTE(RtoDBZ( BYTE_to_mp_func(vol_pol[0][i].ray[k],aMP_SNOW,bMP_SNOW),aMP_class,bMP_class )) ;

                    }
                    else
                        // -- altrimenti correggo comunque a livello liquido :
                        corr=RtoDBZ(vpr_liq,aMP_class,bMP_class)-RtoDBZ(vpr_hray,aMP_class,bMP_class);/*riporto comunque al valore liquido anche se sono sopra la bright band*/

                    // --  controllo qualità su valore correzione
                    if (corr>MAX_CORR) corr=MAX_CORR; /*soglia sulla massima correzione*/
                    if (hbin<hvprmax && corr>0.) corr=0; /*evito effetti incrementi non giustificati*/

                    //controllo qualità su valore corretto e correzione
                    if ( BYTEtoDB(vol_pol[0][i].ray[k])+corr > BYTEtoDB(255) ) // se dato corretto va fuori scala assegno valore massimo
                        vol_pol[0][i].ray[k]=255;
                    else
                        if ( BYTEtoDB(vol_pol[0][i].ray[k])+corr <BYTEtoDB(1)) // se dato corretto va a fodoscala assegno valore di fondo scala
                            vol_pol[0][i].ray[k]=1;
                        else
                            vol_pol[0][i].ray[k]=DBtoBYTE(BYTEtoDB(vol_pol[0][i].ray[k])+corr);  // correggo

                    corr_polar[i][k]=(unsigned char)(corr)+128;


                    //inserisco un ponghino per rifare la neve con aMP e bMP modificati // DA SCOMMENTARE SE DECIDO DI FARLO

                    //if (neve[i][k]) vol_pol[0][i].ray[k]=DBtoBYTE(RtoDBZ( BYTE_to_mp_func(vol_pol[0][i].ray[k],aMP_SNOW,bMP_SNOW),aMP_class,bMP_class )) ;


                }
            }
        }
    }
    return(0);
}

int CUM_BAC::trovo_hvprmax(int *hmax)
{
    int i,imax,istart,foundlivmax;
    float vprmax,h0start,peak,soglia;


    if (! ier_t ){

        printf("trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca\n");
        h0start=t_ground/9.8*1000 ;
        istart=h0start/TCK_VPR -2;
        if (istart< livmin/TCK_VPR) istart=livmin/TCK_VPR;
        printf("t_ground h0start istart %f %f %i  \n",t_ground,h0start,istart);
    }
    else {
        printf("trovo hvprmax  a partire da livmin\n");
        istart=livmin/TCK_VPR+1;
    }


    /* trovo hvprmax e il suo livello a partire dal livello istart */

    //--inizializzazione
    foundlivmax=0;
    peak=NODATAVPR;
    *hmax=INODATA;
    vprmax=NODATAVPR;
    imax=INODATA;
    soglia=DBZtoR(THR_VPR,200,1.6); // CAMBIATO, ERRORE, PRIMA ERA RtoDBZ!!!!VERIFICARE CHE IL NUMERO PARAMETRI FUNZIONE SIA CORRETTO

    //--se vpr al livello corrente e 4 layer sopra> soglia, calcolo picco
    if (vpr[istart] >soglia && vpr[istart+4] > soglia){
        peak=10*log10(vpr[istart]/vpr[istart+4]);//inizializzo il picco
        printf("peak1 = %f\n",peak);
    }
    //----se picco > MINIMO il punto è ok
    if(peak> MIN_PEAK_VPR){
        imax=istart;
        vprmax=vpr[imax];
        printf("il primo punto soddisfa le condizioni di picco \n");
    }
    for (i=istart+1;i<NMAXLAYER-4;i++) //la ricerca è un po' diversa dall'originale.. trovo il picco + alto con valore  rispetto a 4 sopra > soglia
    {
        if (vpr[i] <soglia || vpr[i+4] < soglia) break;
        peak=10*log10(vpr[i]/vpr[i+4]);
        if (vpr[i]>vpr[i-1]  && peak> MIN_PEAK_VPR ) // se vpr(i) maggiore del massimo e picco sufficientemente alto
        {
            imax=i;
            vprmax=vpr[imax];
        }

    }

    if ( imax  > INODATA ){
        foundlivmax=1;
        peak=10*log10(vpr[imax]/vpr[imax+4]);
        *hmax=imax*TCK_VPR+TCK_VPR/2;
        printf("trovato ilaymax %i %i \n",*hmax,imax);
        printf (" picco in dbR %f  \n",peak);
    }

    printf("exit status trovo_hvprmax %i \n",foundlivmax);
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
int CUM_BAC::analyse_VPR(float *vpr_liq,int *snow,float *hliq, const char *sito)
    /*=======analisi profilo============ */
{
    int i,ier=1,ier_ana=0,ier_ap,liv0;
    int tipo_profilo;
    int npar=5;
    float *a,*dyda,v1000,v1500,v600sottobb,vliq,vhliquid,vprmax; //*togliere gli ultimi tre*/;
    char date[]="000000000000";
    struct tm *tempo;
    time_t Time;
    int ndata=10;
    FILE *file;


    // ------------inizializzazione delle variabili ----------

    //strcpy(date,"000000000000");

    tipo_profilo=-1;
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

    //ier_max=trovo_hvprmax(&hvprmax);


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
                fprintf(log_vpr," temperatura alta e non ho trovato hvprmax, esco \n");
                free_vector(a,1,npar);
                free_vector(dyda,1,npar);
                return 1;
            }
            tipo_profilo=0;
        }


        // if (t_ground >= T_MAX_SN+0.65*(float)(livmin+TCK_VPR/2)/100  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100. )
        if (t_ground >= T_MAX_SN  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100. )
        {

            if (  ier_max ) {
                fprintf(log_vpr," temperatura da scioglimento e massimo in quota\n");
                tipo_profilo=2;
            }
            else{
                fprintf(log_vpr," temperatura da scioglimento ma superiore a temperatura max neve e non ho trovato hvprmax, esco \n");
                free_vector(a,1,npar);
                free_vector(dyda,1,npar);
                return 1;
            }
            // solo una scritta per descrivere cos'è accaduto
            liv0=livmin+HALF_BB;
            if (hvprmax > liv0) fprintf(log_vpr," il livello %i è sotto la Bright band, ma T bassa  interpolo\n",livmin);
            else fprintf(log_vpr," il livello %i potrebbe essere dentro la Bright Band, interpolo\n",livmin);

        }

        //if (t_ground >= T_MIN_ML  && t_ground < T_MAX_SN+0.65*(float)(livmin+TCK_VPR/2)/100.)
        if (t_ground < T_MAX_SN)
        {
            if ( ier_max ){
                fprintf(log_vpr," temperatura da neve o scioglimento e massimo in quota\n");
                tipo_profilo=2;
            }
            else {
                fprintf(log_vpr," temperatura da neve o scioglimento e massimo non trovato,neve , non interpolo\n");
                tipo_profilo=3;
                hvprmax=0;
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
                        //lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);
                        if (*hliq<0)
                            *hliq=0;  /*con casi di bright band bassa.. cerco di correggere il più possibile*/
                        *vpr_liq=vpr[(hvprmax+1000)/TCK_VPR]*2.15;
                    }
                    else {
                        *hliq=(a[2]-2.1*a[3])*1000.;
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
    fprintf(log_vpr,"TIPO_PROFILO= %i vpr_liq %f hliq %f \n", tipo_profilo, *vpr_liq,*hliq );


    /* parte di stampa test vpr*/

    /* nome data */
    //definisco stringa data in modo predefinito
    Time = NormalizzoData(old_data_header.norm.maq.acq_date);
    tempo = gmtime(&Time);
    sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
            tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
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

    fprintf(test_vpr,"%s %i %i %f %f %f  %f %f %f %f %f %f %f %f  %f %f %f  %f \n",date,hvprmax,tipo_profilo,stdev,chisqfin,*hliq,vliq,vhliquid,v600sottobb,v1000+6,v1500+6,vprmax,rmsefin,a[1],a[2],a[3],a[4],a[5]);
    fclose(test_vpr);

    /*fine parte di stampa test vpr*/

    free_vector(a,1,npar);
    free_vector(dyda,1,npar);



    /* fine parte di stampa test vpr*/

    //---SCRIVO ALTEZZA MASSIMO PER CLASSIFICAZIONE AL RUN SUCCESSIVO

    ier_ap=access(getenv("VPR_HMAX"),R_OK);
    file = fopen(getenv("VPR_HMAX"),"w");
    fprintf(file,"%d", hvprmax);
    fprintf(log_vpr,"fatta scrittura hmax vpr = %d \n",hvprmax);
    fclose(file);

    return (ier_ana);
}

int CUM_BAC::get_t_ground(float *t_gr)
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

        sito=getenv("SITO");

        if (!(strcmp(sito,"SPC"))) { // spostata qui
            radar_lat=SPC_LAT;
            radar_lon=SPC_LON; }
        if (!(strcmp(sito,"GAT"))) {
            radar_lat=GAT_LAT;
            radar_lon=GAT_LON; }

        while (1) {
            if(fscanf(file_t,"%f %f %f \n",&lon,&lat,&t) == EOF) break;
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
int CUM_BAC::interpola_VPR(float a[], int ma)
# include <vpr_par.h>
{
    float *x, *y,*sig,alamda,y1=0,*dyda,B,E,C,G,F,xint,qdist,*abest;
    float chisq=100.;
    float chisqold=0.0;
    float chisqin=0.0;
    int i,in1,in2,in3,in4,*ia,ifit,ii,ndati_nok,ndati_ok,k,ier_int;
    //int ma=5;
    int ndata=10;
    float **covar;
    float **alpha;
    char file_vprint[200];
    FILE *file;

    fprintf(log_vpr,"sono in interpola_vpr\n");
    ier_int=0;

    in1=(hvprmax-TCK_VPR/2)/TCK_VPR; //indice del massimo
    in2=(hvprmax+HALF_BB)/TCK_VPR; //indice del massimo + 500 m
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

        alamda=-0.01;

        for (i=1;i<=ma;i++) ia[i]=1;
        qdist=0;
        ii=1;
        ndati_nok=0;

        for (i=1; i<=ndata; i++)
        {
            sig[ii]=0.5;
            x[ii]= ((hvprmax-1000.)>livmin)? (i*TCK_VPR+(hvprmax-800)-TCK_VPR)/1000. : (livmin+(i-1)*TCK_VPR)/1000.;
            y[ii]= ((hvprmax-1000.)>livmin)? vpr[i+((hvprmax-800)-TCK_VPR)/TCK_VPR] : vpr[i-1+livmin/TCK_VPR];
            // x[ii]= ((hvprmax-800.)>livmin)? (i*TCK_VPR+(hvprmax-600)-TCK_VPR)/1000. : (livmin+(i-1)*TCK_VPR)/1000.;
            //y[ii]= ((hvprmax-800.)>livmin)? vpr[i+((hvprmax-600)-TCK_VPR)/TCK_VPR] : vpr[i-1+livmin/TCK_VPR];
            lineargauss(x[ii], a, &y1, dyda, ndata);
            qdist=(y1-y[ii])*(y1-y[ii]);
            if (sqrt(qdist) < DIST_MAX)
            {
                ii+=1;
                chisqin=qdist+chisqin;
            }
            else
                ndati_nok=ndati_nok+1;

            if    ( ndati_nok > 2  )
            {
                fprintf (log_vpr,"  first guess troppo lontano dai punti , interpolazione fallisce \n");
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


            if (chisq <chisqfin)
            {
                chisqfin=chisq;
                for (i=1;i<=ma;i++) abest[i]=a[i];
            }
        }
    }

    for (i=1; i<=ndata-ndati_nok; i++)
    {
        lineargauss(x[i], abest, &y1, dyda, ndata);
        rmsefin=rmsefin+  (y[i]-y1)*(y[i]-y1) ;
    }
    rmsefin=sqrt(rmsefin/(float)((ndata-ndati_nok)*(ndata-ndati_nok)));
    fprintf(log_vpr,"RMSEFIN %f \n ", rmsefin );


    if (chisqfin>CHISQ_MAX)
    {
        ier_int=1;
    }
    else {
        for (i=1;i<=ma;i++) a[i]=abest[i];
        sprintf(file_vprint,"%s_int",getenv("VPR_ARCH"));
        file=controllo_apertura(file_vprint," vpr interpolato ","w");
        for (i=1; i<=NMAXLAYER; i++)
        {
            xint=(i*TCK_VPR-TCK_VPR/2)/1000.;
            lineargauss(xint, a, &y1, dyda, ndata);
            // stampa del profilo interpolato
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

int CUM_BAC::testfit(float a[], float chisq, float chisqin)
{
    if (a[1]<0. || a[1] >15.) return 1;
    if (a[2] >10.) return 1;
    if (a[3]<0.2 || a[3] > 0.6 ) return 1; //da analisi set dati
    if (a[4]<0. ) return 1;
    if (a[5]>0 ) return 1;
    if (chisq>chisqin ) return 1;
    return 0;
}

#endif

#ifdef CLASS
void CUM_BAC::classifica_rain()
{
    float a;// raggio terra, non so perchè lo rendo variabile
    float range[MAX_BIN];
    float zz[MAX_BIN][NEL];
    float xx[MAX_BIN][NEL];
    int  i_xx[MAX_BIN][NEL],i_zz[MAX_BIN][NEL],i_xx_min[MAX_BIN][NEL],i_xx_max[MAX_BIN][NEL],i_zz_min[MAX_BIN][NEL],i_zz_max[MAX_BIN][NEL];
    int  im[MAX_BIN][NEL], ix[MAX_BIN][NEL], jm[MAX_BIN][NEL], jx[MAX_BIN][NEL];
    int i,j,kx,kz,k,iel,iaz,ibin,imin,imax,jmin,jmax,RHI_ind[NEL][MAX_BIN],jbb;
    int wimin,wjmin,wimax,wjmax;
    int hmax=-9999, ier_ap,ier_0term=0;

    //float w_size[2]={3.,1.5}; //dimensione della matrice pesi
    float w_size[2]={3.,0.3}; //dimensione della matrice pesi
    float **rhi_cart,**rhi_weight,RHI_beam[NEL][MAX_BIN],*w_x,*w_z,**w_tot,**beamXweight[MAX_BIN]; // da inizializzare in fase di programma
    float range_min,range_max,xmin,zmin,xmax,zmax;
    int w_x_size,w_z_size,w_x_size_2,w_z_size_2;
    FILE *file;

    //definisco e così inizializzo resol e a
    resol[0]=RES_HOR_CIL; // uguale a dimensione cella volume polare .. va parametrizzato
    resol[1]=RES_VERT_CIL;
    a=REARTH;
    // inizializzazione variabili



    /* ;---------------------------------- */
    /* ;          FASE 0 :                  */
    /* ;---------------------------------- */
    // DEFINISCO QUOTE DELLA BASE E DEL TOP DELLA BRIGHT BAND USANDO IL DATO quota del picco  DEL PRECEDENTE RUN O, SE NON PRESENTE LA QUOTA DELLO ZERO DA MODELLO

    // Lettura quota massimo da VPR  calcolo base e top bright band
    fprintf(log_class,"data= %s \n",date);
    // calcolo il gap
    gap=profile_gap(getenv("LAST_VPR"));
    //-- se gap < memory leggo hmax da VPR
    if (gap<=MEMORY){
        ier_ap=access(getenv("VPR_HMAX"),R_OK);
        if (!ier_ap) {
            file = fopen(getenv("VPR_HMAX"),"r");
            //file=controllo_apertura(getenv("VPR_HMAX")," altezza hmax ultimo vpr ","r");
            fscanf(file,"%i", &hmax);
            fprintf(log_class,"fatta lettura hmax vpr = %i \n",hmax);
            printf("fatta lettura hmax vpr = %i \n",hmax);
            //---suppongo una semiampiezza massima della bright band di 600 m e definisco htopbb e hbasebb come hmassimo +600 m (che da clima ci sta) e hmassimo -600 m
            hbbb=(hmax-600.)/1000.;
            htbb=(hmax+600.)/1000.;
            fclose(file);
        }
    }
    //-- se gap > memory o se non ho trovato il file
    if (hmax <0 ){
        // Lettura 0 termico da modello, e calcolo base e top bright band
        fprintf(log_class,"leggo 0termico per class da file %s \n",getenv("FILE_ZERO_TERMICO"));
        //  leggo informazioni di temperatura da modello*/
        ier_0term=trovo0term();
        if( !ier_0term) {
            //-- considerato che lo shift medio tra il picco e lo zero è tra 200 e 300 m, che il modello può avere un errore, definisco cautelativamente htbb come quota zero + 400 m e hbbb come  quota zero -700 m  .
            htbb=zeroterm/1000. + 0.4; // se non ho trovato il vpr allora uso un range più ristretto, potrebbe essere caso convettivo
            hbbb=zeroterm/1000. - 1.0;

        }
        else{
            printf("attenzione, non ho trovat zero termico ne da vpr ne da radiosondaggio \n");
            htbb=0; // discutibile così faccio tutto con VIZ
            hbbb=0;
        }
    }

    // se hbasebb è <0 metto 0
    if (hbbb<0) hbbb=0;

    fprintf(log_class,"calcolati livelli sopra e sotto bright band hbbb=%f  htbb=%f \n",hbbb,htbb);
    /* ---------------------------------- */
    /*           FASE 1 */
    /* ---------------------------------- */
    /*    Costruzione matrice generale per indicizzare il caricamento dei dati */
    /*    da coordinate radar a coordinate (X,Z) */
    /*  Calcolo distanza per ogni punto sul raggio */
    /*  xx contiene la distanza sulla superfice terrestre in funzione del range radar e dell'elevazione */
    /*  Metodo 1 - -Calcolo le coordinate di ogni punto del RHI mediante utilizzo di cicli */

    // estremi x e z (si procede per rhi)
    range_min=0.5*size_cell[old_data_header.norm.maq.resolution]/1000.;
    range_max=(MAX_BIN-0.5)*size_cell[old_data_header.norm.maq.resolution]/1000.;

    xmin=floor(range_min*cos(elev_array[NEL-1]*CONV_RAD)); // distanza orizzontale minima dal radar
    zmin=pow(pow(range_min,2.)+pow(4./3*a,2.)+2.*range_min*4./3.*a*sin(elev_array[0]*CONV_RAD),.5) -4./3.*a+h_radar; // quota  minima in prop standard
    xmax=floor(range_max*cos(elev_array[0]*CONV_RAD)); // distanza orizzontale massima dal radar
    zmax=pow(pow(range_max,2.)+pow(4./3*a,2.)+2.*range_max*4./3.*a*sin(elev_array[NEL-1]*CONV_RAD),.5) -4./3.*a+h_radar;//quota massima

    x_size=(xmax-xmin)/resol[0]; //dimensione orizzontale
    z_size=(zmax-zmin)/resol[1]; //dimensione verticale
    fprintf(log_class,"calcolati range_min e range_max , dimensione orizzontale e dimensione verticale range_min=%f  range_max=%f x_size=%d z_size=%d \n",range_min,range_max,x_size,z_size);
    if (x_size > MAX_BIN) x_size=MAX_BIN;

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

        for (k=0; k<NEL; k++){
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
                cil[i][k][j]= -20. ;
        }
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

        for (i=0;i<NEL;i++){
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

        for (iel=0;iel<NEL;iel++){
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
                jbb=ceil((htbb)/RES_HOR_CIL);
                if (j == jbb ) cappi[iaz][i] = DBtoBYTE(cil[iaz][i][j]);

            }
        }
    }


    //output = fopen("CAPPI","w");
    //fwrite(cappi,sizeof(cappi),1,output);
    //fclose(output);

    //-------------------------------------------------------------------------------------------------------------------------
    // faccio la classificazione col metodo Vertical Integrated Reflectivity
    classifico_VIZ();
    //classificazione con STEINER
    //  if (hmax > 2000.) {// per evitare contaminazioni della bright band, si puo' tunare
    // if (hbbb > 500.) {// per evitare contaminazioni della bright band, si puo' tunare
    calcolo_background();
    classifico_STEINER();
    //  }
    merge_metodi();
    return ;
}


void CUM_BAC::classifico_VIZ()
{
    int i,j,k,kbbb=0,ktbb=0,kmax=0;
    float cil_Z,base,Zr;
    double *Zabb[NUM_AZ_X_PPI],*Zbbb[NUM_AZ_X_PPI], ext_abb,ext_bbb;
    float LIM_VERT= 8.;//questo l'ho messo io
    FILE *file_Zabb;


    kbbb=floor(hbbb/resol[1]);   //08/01/2013...MODIFICA, inserito questo dato
    ktbb=ceil(htbb/resol[1]);
    kmax=ceil(LIM_VERT/resol[1]);
    // kmax=ceil(z_size/resol[1]);
    if (t_ground < T_MAX_ML) kmax=0;/////se t suolo dentro t melting layer pongo kmax=00 e in tal modo non classifico
    if (ktbb>z_size) ktbb=z_size;
    printf ("kmax= %i \n kbbb= %i \n ktbb= %i \n  z_size= %i \n ",kmax,kbbb,ktbb,z_size);

    //inizializzazione vettori e matrici
    for (i=0; i<NUM_AZ_X_PPI; i++){
        Zabb[i]= (double *) malloc(x_size*sizeof(double ));
        Zbbb[i]= (double *) malloc(x_size*sizeof(double ));
        conv_VIZ[i]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;
        conv_STEINER[i]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;
        conv[i]=(unsigned char *) malloc(x_size*sizeof(unsigned char )) ;



        for (j=0; j<x_size; j++){ // cambiato da x_size
            Zabb[i][j]=0.;
            Zbbb[i][j]=0.;
            conv_VIZ[i][j]=MISSING;
            conv_STEINER[i][j]=MISSING;
            conv[i][j]=MISSING;
            stratiform[i][j]=MISSING;
        }
    }

    //inizio l'integrazione
    for(i=0; i<NUM_AZ_X_PPI; i++){
        for(j=0; j<x_size; j++)
        {
            ext_abb=0.;
            ext_bbb=0.;

            //modifica 08/01/2013 .. afggiungo questo ..per fare l'integrazione anche con i dati sotto la bright band
            for(k=0; k<kbbb; k++)
            {
                if (cil[i][j][k] > -19.){   // 08/01/2013..modifica, prendo fin dove ho un segnale
                    base=(cil[i][j][k])/10.;
                    cil_Z=pow(10.,base);
                    Zbbb[i][j] = Zbbb[i][j] + resol[1]*cil_Z;
                    ext_bbb=resol[1]+ext_bbb;
                }
            }
            for(k=kbbb; k<ktbb; k++)
            {
                if (k < 4 ){
                    if (cil[i][j][k]>10. &&  cil[i][j][k+4]> 5.){
                        if (cil[i][j][k] - cil[i][j][k+4] > 5.)
                            stratiform[i][j]=1;
                    }
                }
                else if (cil[i][j][k]>10. &&  cil[i][j][k+4]> 5. &&  cil[i][j][k-4] > 5.){
                    if (cil[i][j][k] - cil[i][j][k+4] > 5.&&   cil[i][j][k]- cil[i][j][k-4] > 5. )
                        stratiform[i][j]=1;

                }


                if (cil[i][j][k] - cil[i][j][k+4] > 5.)
                    stratiform[i][j]=1;

                for(k=ktbb; k<kmax; k++)
                {
                    if (cil[i][j][k] > -19.){    // 08/01/2013..modifica, prendo fin dove ho un segnale
                        base=(cil[i][j][k])/10.;
                        cil_Z=pow(10.,base);
                        Zabb[i][j] = Zabb[i][j] + resol[1]*cil_Z;
                        ext_abb=resol[1]+ext_abb;
                        //if (j == 400) printf(" %i  %i %f  %f  %f \n", i, j, cil[i][j][k],Zabb[i][j],ext ) ;
                    }

                }


                //solo se l'estensione verticale del segnale sopra il top della bright band è maggiore di 0.8 Km classifico

                if (ext_bbb +  ext_abb>0.8) {
                    //if ( ext_abb>0.8) {
                    if ((Zabb[i][j] +Zbbb[i][j])/(ext_bbb+ext_abb) > THR_VIZ){
                        //if ((Zabb[i][j] /ext_abb) > THR_VIZ){
                        conv_VIZ[i][j]=CONV_VAL;
                        lista_conv[ncv][0]= i;
                        lista_conv[ncv][1]= j;
                        ncv=ncv+1;
                    }
                }

            }
        }
    }
    printf("numero nuclei VIZ = %li \n",ncv);

    return;
}


void  CUM_BAC::classifico_STEINER()
{
    static bool warned = false;
    int i,j,k;

    unsigned char BYTE;
    float diff_bckgr,cr;

    for(i=0; i<np; i++){
        j=lista_bckg[i][0]; //az=lista_bckg[i][0]
        k=lista_bckg[i][1]; //ra=lista_bckg[i][1]

        // calcolo valore BYTE nel punto
        if (j == -999)
        {
            if (!warned)
            {
                fprintf(stderr, "elev_fin[%d][%d]\n", j, k);
                fprintf(stderr, " = %d\n", (int)elev_fin[j][k]);
                fprintf(stderr, "vol_pol[%d][%d].ray[%d]\n", (int)elev_fin[j][k], j, k);
                //fprintf(stderr, " = %d\n", (int)vol_pol[elev_fin[j][k]][j].ray[k]);
                warned = true;
            }
            BYTE=56;
#warning This code used to work like this, by accident, in C: this is a workaround to maintain functional equivalence during porting to c++
        } else
            BYTE=vol_pol[elev_fin[j][k]][j].ray[k];
        // calcolo diff col background
        diff_bckgr=BYTEtoDB(BYTE)-bckgr[i];
        /* if (k < 160 ){ */
        /*   fprintf(log_class," bckgr[i] = %f diff_bckgr %f \n ",bckgr[i],diff_bckgr ); */
        /* } */
        // test su differenza con bckground , se soddisfatto e simultaneamente il VIZ non ha dato class convettiva (?)
        if ((BYTEtoDB(BYTE)>40.)||
                (bckgr[i]< 0 && diff_bckgr > 10) ||
                (bckgr[i]< 42.43 &&  bckgr[i]>0 &&  diff_bckgr > 10. - bckgr[i]*bckgr[i]/180. )||
                (bckgr[i]> 42.43 &&  diff_bckgr >0)  ) {
            // assegno il punto  nucleo di Steiner
            conv_STEINER[j][k]=CONV_VAL;

            // ingrasso il nucleo
            cr=convective_radius[i];
            printf (" %f cr \n", cr);
            ingrasso_nuclei(cr,j,k);
            ncs=ncs+1;

        }
    }

    return;
}

void CUM_BAC::calcolo_background() // sui punti precipitanti calcolo bckgr . nb LA CLASSIFICAZIONE DI STEINER NON HA BISOGNO DI RICAMPIONAMENTO CILINDRICO PERCIÒ uso direttamente la matrice polare
    // definisco una lista di punti utili per l'analisi, cioè con valore non nullo e sotto la bright band o sopra (NB. per l'analisi considero i punti usati per la ZLR che si suppongono non affetti da clutter e beam blocking<50%. questo ha tutta una serie di implicazioni..... tra cui la proiezione, il fatto che nel confronto entrano dei 'buchi'...cioè il confronto è fatto su una matrice pseudo-orizzontale bucata e quote variabili tuttavia, considerato che i gradienti verticali in caso convettivo e fuori dalla bright band non sono altissimi spero funzioni ( andro' a verificare che l'ipotesi sia soddisfatta)
    // per trovare il background uso uno pseudo quadrato anzichè cerchio, mantenendo come semi-lato il raggio di steiner (11KM); devo perciò  definire una semi-ampiezza delle finestre in azimut e range corrispondente al raggio di 11 km
    // quella in azimut  dipende dal range

{
    int i,j,k,jmin,jmax,J_MAX,J_BBB,J_ABB,kmin,kmax,delta_naz=0,delta_nr;
    float z,zmax,range,delta_r,delta_az,az_min,az_max;
    long int npoints=0;




    //traccio una lista dei punti che hanno valore non nullo e sotto base bright band (lista_bckg) contenente iaz e irange e conto i punti
    for (i=0; i<NUM_AZ_X_PPI;i++) {
        for (j=0; j<MAX_BIN;j++)  // propongo max_bin visto che risoluzione è la stessa
        {
            //if ( vol_pol[0][i].ray[j] > 1 &&  (float)(quota[i][j])/1000. < hbbb ) //verifico che il dato usato per la ZLR cioè la Z al lowest level sia > soglia e la sua quota sia sotto bright band o sopra bright band

            {     if ( vol_pol[0][i].ray[j] > 1 )
                lista_bckg[np][0]=i;  //IAZIMUT
                lista_bckg[np][1]=j;  //IRANGE
                np=np+1;
            }
        }

    }

    // per il calcolo della finestra range su cui calcolare il background divido il raggio di Steiner (11km) per la dimensione della cella
    delta_nr=(int)(STEINER_RADIUS*1000./size_cell[old_data_header.norm.maq.resolution]);//definisco ampiezza semi-finestra range corrispondente al raggio di steiner (11km), unità matrice polare
    printf("delta_nrange per analisi Steiner = %i \n",delta_nr);



    if (np > 1) {
        //inizializzo vettori
        convective_radius=(float *)malloc(np*sizeof(float));
        Z_bckgr=(double *) malloc(np*sizeof(double));
        bckgr=(float *) malloc(np*sizeof(float));

        for (i=0;i<np;i++){ //M: tolto np-1 messo np
            bckgr[i]=0.;
            Z_bckgr[i]=0;
            convective_radius[i]=0;
        }

        for(i=0; i<np;i++){       // M:tolto np -1 messo np
            npoints=0;

            // estremi della finestra in range
            kmin=lista_bckg[i][1]-delta_nr;
            kmax=lista_bckg[i][1]+delta_nr;

            if (kmin>0) {

                if (kmax>MAX_BIN) kmax=MAX_BIN;

                //definisco ampiezza semi finestra nazimut  corrispondente al raggio di steiner (11km)  (11/distanzacentrocella)(ampiezzaangoloscansione)
                delta_naz=ceil(11./((lista_bckg[i][1]*size_cell[old_data_header.norm.maq.resolution]/1000.+size_cell[old_data_header.norm.maq.resolution]/2000.)/(AMPLITUDE*DTOR)));
                if (delta_naz>NUM_AZ_X_PPI/2)
                    delta_naz=NUM_AZ_X_PPI/2;

                jmin=lista_bckg[i][0]-delta_naz;
                jmax=lista_bckg[i][0]+delta_naz;


                if (jmin<0) {
                    jmin=NUM_AZ_X_PPI-jmin%NUM_AZ_X_PPI;
                    for (j= jmin  ; j< NUM_AZ_X_PPI ; j++) {
                        for (k= kmin ; k< kmax  ; k++){
                            //        if ( vol_pol[elev_fin[j][k]][j].ray[k]> 1 &&  (float)(quota[j][k])/1000. < hbbb ) {  // aggiungo condizione quota
                            if ( vol_pol[elev_fin[j][k]][j].ray[k]> 1  ){
                                Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]) ;
                                bckgr[i]=bckgr[i]+BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k]);
                                npoints=npoints+1;
                            }
                        }
                    }
                    jmin=0;
                }

                if (jmax>NUM_AZ_X_PPI) {
                    jmax=jmax%NUM_AZ_X_PPI;
                    for (j= 0  ; j< jmax ; j++) {
                        for (k= kmin ; k< kmax  ; k++){
                            // if (vol_pol[elev_fin[j][k]][j].ray[k]> 1 &&  (float)(quota[j][k])/1000. < hbbb ) {
                            if ( vol_pol[elev_fin[j][k]][j].ray[k]> 1  ) {
                                Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);
                                bckgr[i]=bckgr[i]+BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k]);
                                npoints=npoints+1;
                            }
                        }
                        }
                        jmax=NUM_AZ_X_PPI;
                }

                for (j=jmin   ; j<jmax  ; j++) {
                    for (k=kmin  ; k<kmax   ; k++){
                        // if (vol_pol[elev_fin[j][k]][j].ray[k]> 1 &&  (float)(quota[j][k])/1000. < hbbb ) {
                        if ( vol_pol[elev_fin[j][k]][j].ray[k]> 1  ) {
                            Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);
                            bckgr[i]=bckgr[i]+BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k]);
                            npoints=npoints+1;
                        }
                    }
                }
            }
            else{
                for (j=0   ; j<NUM_AZ_X_PPI/2  ; j++){
                    for (k=0  ; k<kmax   ; k++){
                        // if (vol_pol[elev_fin[j][k]][j].ray[k]> 1 &&  (float)(quota[j][k])/1000. < hbbb ) {
                        if ( vol_pol[elev_fin[j][k]][j].ray[k]> 1  ) {
                            Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);
                            bckgr[i]=bckgr[i]+BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k]);
                            npoints=npoints+1;
                        }
                    }
                }
                for (j= NUM_AZ_X_PPI/2  ; j<NUM_AZ_X_PPI  ; j++) {
                    for (k=0  ; k<-kmin   ; k++){
                        // if (vol_pol[elev_fin[j][k]][j].ray[k]> 1 &&  (float)(quota[j][k])/1000. < hbbb ) {
                        if ( vol_pol[elev_fin[j][k]][j].ray[k]> 1  ) {
                            Z_bckgr[i]=Z_bckgr[i]+ BYTEtoZ(vol_pol[elev_fin[j][k]][j].ray[k]);
                            bckgr[i]=bckgr[i]+BYTEtoDB(vol_pol[elev_fin[j][k]][j].ray[k]);
                            npoints=npoints+1;
                        }
                    }
                }
            }
            if (npoints > 0){
                Z_bckgr[i]=Z_bckgr[i]/npoints;
                //bckgr[i]=bckgr[i]/npoints; //no
                if (Z_bckgr[i]>0) bckgr[i]=10*(log10(Z_bckgr[i]));
            }
            //il valore del raggio convettivo varia a seconda del background, da 1 a 5 km
            if (  bckgr[i] < 25.) convective_radius[i] = 1.;
            if (  bckgr[i] >= 25. && bckgr[i] <30. ) convective_radius[i] = 2.;
            if (  bckgr[i] >= 30. && bckgr[i] <35. ) convective_radius[i] = 3.;
            if (  bckgr[i] >= 35. && bckgr[i] <40. ) convective_radius[i] = 4.;
            if (  bckgr[i] > 40.)  convective_radius[i] = 5.;

        }

    }

    return;
}

void CUM_BAC::ingrasso_nuclei(float cr,int ja,int kr)
{
    int daz, dr,jmin,jmax,kmin,kmax,j,k;


    dr=(int)(cr*1000./size_cell[old_data_header.norm.maq.resolution]);//definisco ampiezza semi-finestra range corrispondente al raggio di steiner (11km), unità matrice polare

    kmin=kr-dr;
    kmax=kr+dr;

    daz=ceil(cr/((kr*size_cell[old_data_header.norm.maq.resolution]/1000.+size_cell[old_data_header.norm.maq.resolution]/2000.)/(AMPLITUDE*DTOR)));
    jmin=ja-daz;
    jmax=ja+daz;

    printf ("dr cr kmin kmax  %d %f %d %d %d %d \n", dr,cr, kmin,kmax,jmin,jmax);

    if (kmin>0) {
        if (kmax>x_size) kmax=x_size;

        if (jmin<0) {
            jmin=NUM_AZ_X_PPI-jmin%NUM_AZ_X_PPI;
            for (j=jmin; j< NUM_AZ_X_PPI ; j++) {
                for (k=kmin ; k<kmax  ; k++) {
                    conv_STEINER[j][k]=CONV_VAL;
                }
            }
            printf ("jmin %d \n", jmin);
            jmin=0;

        }

        if (jmax>NUM_AZ_X_PPI) {
            jmax=jmax%NUM_AZ_X_PPI;
            for (j=0; j<jmax ; j++) {
                for (k=kmin; k<kmax  ; k++) {
                    conv_STEINER[j][k]=CONV_VAL;
                }
            }
            printf ("jmax %d \n", jmax);
            jmax=NUM_AZ_X_PPI;
        }
        for (j=jmin; j<jmax ; j++) {
            for (k=kmin; k<kmax  ; k++) {
                conv_STEINER[j][k]=CONV_VAL;
            }
        }
    }
    else
    {
        for (j=0   ; j<NUM_AZ_X_PPI/2  ; j++)
            for (k=0  ; k<kmax   ; k++){
                conv_STEINER[j][k]=CONV_VAL;
            }
        for (j= NUM_AZ_X_PPI/2  ; j<NUM_AZ_X_PPI  ; j++)
            for (k=0  ; k<-kmin   ; k++){
                conv_STEINER[j][k]=CONV_VAL;
            }

    }

    return;
}

void CUM_BAC::merge_metodi()
{
    int j,k;

    for(j=0; j<NUM_AZ_X_PPI; j++){
        for(k=0; k<x_size; k++)
        {

            if ( conv_STEINER[j][k] == conv_VIZ[j][k] &&  conv_STEINER[j][k]> 0 && stratiform[j][k]<1 ){
                conv[j][k]=conv_VIZ[j][k];
            }

        }
    }
    return;
}

int CUM_BAC::trovo0term()
{

    int ier;
    FILE *file0;
    ier=access(getenv("FILE_ZERO_TERMICO"),R_OK);
    if (!ier) {
        file0=fopen(getenv("FILE_ZERO_TERMICO"),"r");//metti condizione di non trovato!!!
        fscanf(file0,"%f", &zeroterm );
        fprintf(log_class," zero termico da modello= %f \n ",zeroterm);
        fclose(file0);
    }
    else {
        fprintf(log_class,"non ho trovato il file dello zero termico \n ");
    }
    return ier;
}
#endif

bool CUM_BAC::esegui_tutto(const char* nome_file, int file_type, const char* sito)
{
    // Legge e controlla il volume dal file SP20
    if (!read_sp20_volume(nome_file, sito, file_type))
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

    if (NormalizzoData(old_data_header.norm.maq.acq_date) == -1)
        return true;

    /*------------------------------------------
      | rimozione propagazione anomala e clutter |
      ------------------------------------------*/
    LOG_INFO("%s %s -- Cancellazione Clutter e Propagazione Anomala", PrendiOra(), nome_file);

    // --- ricavo il mese x definizione first_level e  aMP bMP ---------
    //definisco stringa data in modo predefinito
    Time = NormalizzoData(old_data_header.norm.maq.acq_date);
    tempo = gmtime(&Time);
    month=tempo->tm_mon+1;

    // scrivo la variabile char date con la data in formato aaaammgghhmm
    sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
            tempo->tm_mday,tempo->tm_hour, tempo->tm_min);

    // ------setto ambiente statico (le first level e il nome del dem )--------------
    setstat(sito,month, nome_dem, nome_fl);
    LOG_INFO("nome dem %s nome fl %s", nome_dem, nome_fl);

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


    //--------------se def anaprop : rimozione propagazione anomala e correzione beam blocking-----------------//
#ifdef ANAPROP
    LOG_INFO("inizio rimozione anaprop e beam blocking");
    int ier = elabora_dato();
#endif


    //--------------se definita la qualita procedo con il calcolo qualita e del VPR (perchè prendo solo i punti con qual > soglia?)-----------------//

#ifdef QUALITY

    //-------------calcolo qualita' e trovo il top
    printf ("calcolo Q3D \n") ;
    caratterizzo_volume();

    /* //---------trovo il top (a X dbZ) */
    /* printf ("trovo top \n") ; */
    /* ier=trovo_top(); */


    //--------------se definito VPR procedo con ricerca t_ground che mi serve per classificazione per cui la metto prima-----------------//
#ifdef VPR
    log_vpr=fopen(getenv("LOG_VPR"),"a+");
    ier_t=get_t_ground(&t_ground);
    fclose(log_vpr);
#endif


    //--------------se definita CLASS procedo con  classificazione -----------------//
#ifdef CLASS
    log_vpr=fopen(getenv("LOG_VPR"),"a+");
    log_class=fopen(getenv("LOG_CLASS"),"a+");
    classifica_rain();
    fclose(log_vpr);
    fclose(log_class);
#endif


    //--------------se definito VPR procedo con calcolo VPR -----------------//

#ifdef VPR
    log_vpr=fopen(getenv("LOG_VPR"),"a+");
    test_vpr=fopen(getenv("TEST_VPR"),"a+");

    fprintf(log_vpr,"processo file dati: %s\n",nome_file);
    printf ("calcolo VPR \n") ;

    //VPR  // ------------inizializzo hvprmax ---------------

    hvprmax=INODATA;

    //VPR  // ------------chiamo combina profili con parametri sito, sito alternativo ---------------

    //  ier_comb=combina_profili(sito,argv[4]);
    ier_comb=combina_profili(sito);
    printf ("exit status calcolo VPR istantaneo: (1--fallito 0--ok)  %i \n",ier_vpr) ; // debug
    printf ("exit status combinaprofili: (1--fallito 0--ok) %i \n",ier_comb) ; // debug


    //VPR  // ------------chiamo profile_heating che calcola riscaldamento profilo ---------------

    heating=profile_heating();
    printf ("heating %i \n", heating);
    fprintf (log_vpr,"ier_vpr %i ier_comb %i \n",ier_vpr,ier_comb);

    //VPR  // ------------se combina profili ok e profilo caldo correggo --------------
    if (!ier_comb && heating >= WARM){

        ier=corr_vpr(sito);
        printf ("exit status correggo vpr: (1--fallito 0--ok) %i \n",ier) ; // debug


        //VPR // ------------se la correzione è andata bene e il profilo è 'fresco' stampo profilo con data-------

        if ( ! ier && ! ier_vpr)
            ier_stampa_vpr=stampa_vpr();
    }

    fclose(log_vpr);

#endif
#endif

#ifdef CLASS
    log_class=fopen(getenv("LOG_CLASS"),"a+");
    for (int i=0; i<NUM_AZ_X_PPI; i++){
        for (int k=0; k<vol_pol[0][i].b_header.max_bin; k++){

            if (conv[i][k] > 0){

                vol_pol[0][i].ray[k]=DBtoBYTE(RtoDBZ( BYTE_to_mp_func(vol_pol[0][i].ray[k],aMP_conv,bMP_conv),aMP_class,bMP_class )) ;

            }
        }

    }
    fclose(log_class);
#endif


#ifdef TIME
    prendo_tempo();
#endif


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
#ifndef WRITE_DBP

    /*--------------------------------------------------
      | conversione di coordinate da polare a cartesiana |
      --------------------------------------------------*/
    LOG_INFO("%s -- Creazione Matrice Cartesiana",PrendiOra());
    creo_cart();

    //------------------- STA PRENDOTEMPO!!!!!!!!!!!!!!!
#ifdef TIME
    prendo_tempo();
#endif



    //-------------------Se definita Z_LOWRIS creo matrice 1X1  ZLR  stampo e stampo coeff MP (serve?)------------------

#ifdef Z_LOWRIS

    LOG_INFO("%s -- Estrazione Precipitazione 1X1", PrendiOra());
    creo_cart_z_lowris();

    //-------------------scritture output -----------------------
    LOG_INFO("%s -- Scrittura File Precipitazione 1X1 %s\n", PrendiOra(), nome_file);
    scrivo_out_file_bin(".ZLR","dati output 1X1",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(z_out),z_out);


    char nome_file_output[512];
    sprintf(nome_file_output,"%s/MP_coeff",getenv("OUTPUT_Z_LOWRIS_DIR"));
    FILE* output=controllo_apertura(nome_file_output,"file coeff MP","w");
    fwrite(MP_coeff,sizeof(MP_coeff),1,output);
    fclose(output);
    printf(" dopo scrivo_z_lowris\n");



#ifdef TIME
    prendo_tempo();
#endif





    //-------------------scritture output -----------------------
#ifdef QUALITY

    //------------------ output qual  per operativo:qualità in archivio e elevazioni e anap in dir a stoccaggio a scadenza
    // temporanee
    // in archivio
    scrivo_out_file_bin(".qual_ZLR","file qualita' Z",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(qual_Z_1x1),qual_Z_1x1);

    //------------------ stampe extra per studio
#ifdef STAMPE_EXTRA
    //scrivo_out_file_bin(".corrpt","file anap",getenv("DIR_QUALITY"),sizeof(dato_corrotto),dato_corrotto);
    //scrivo_out_file_bin(".pia","file PIA",getenv("DIR_QUALITY"),sizeof(att_cart),att_cart);
    //scrivo_out_file_bin(".bloc","file bloc",getenv("DIR_QUALITY"),sizeof(beam_blocking),beam_blocking);
    //scrivo_out_file_bin(".quota","file quota",getenv("DIR_QUALITY"),sizeof(quota),quota);
    scrivo_out_file_bin(".elev","file elevazioni",getenv("DIR_QUALITY"),sizeof(elev_fin),elev_fin);
    scrivo_out_file_bin(".bloc_ZLR","file bloc",getenv("DIR_QUALITY"),sizeof(beam_blocking_1x1),beam_blocking_1x1);
    scrivo_out_file_bin(".anap_ZLR","file anap",getenv("DIR_QUALITY"),sizeof(dato_corr_1x1),dato_corr_1x1); //flag di propagazione anomala
    scrivo_out_file_bin(".quota_ZLR","file qel1uota",getenv("DIR_QUALITY"),sizeof(quota_1x1),quota_1x1);// m/100 +128
    scrivo_out_file_bin(".elev_ZLR","file elev",getenv("DIR_QUALITY"),sizeof(elev_fin_1x1),elev_fin_1x1);
    scrivo_out_file_bin(".top20_ZLR","file top20",getenv("DIR_QUALITY"),sizeof(top_1x1),top_1x1);

#endif
    //------------------ stampe correzioni da profili verticali in formato ZLR
#ifdef VPR
    scrivo_out_file_bin(".corr_ZLR","file correzione VPR",getenv("DIR_QUALITY"),sizeof(corr_1x1),corr_1x1);
    //scrivo_out_file_bin(".neve","punti di neve",getenv("DIR_QUALITY"),sizeof(neve),neve);
    // scrivo_out_file_bin(".neve_ZLR","file presunta neve ",getenv("DIR_QUALITY"),sizeof(neve_1x1),neve_1x1);
#endif

    //------------------se definita CLASS  stampo punti convettivi
#ifdef CLASS
    // in archivio?
    // scrivo_out_file_bin(".conv_ZLR","punti convettivi",getenv("OUTPUT_Z_LOWRIS_DIR"),sizeof(conv_1x1),conv_1x1);
    scrivo_out_file_bin(".conv_ZLR","punti convettivi",getenv("DIR_QUALITY"),sizeof(conv_1x1),conv_1x1);

#endif

#endif


#endif// ifdef Z_lOWRIS
#endif // ifndef WRITE_DBP

    return true;
}

  /*non cancellare*/
  /* int trovo_top() */
  /* { */
  /*   int ier,i,k,l; */
  /*   FILE *file0; */

  /*   for(i=0; i<NUM_AZ_X_PPI; i++){ */
  /*     for (k=0; k<vol_pol[0][i].b_header.max_bin; k++){ */
  /*       for (l=first_level_static[i][k]; l<NEL; l++)     */
  /*    {  */
  /*      if (BYTEtoDB(vol_pol[l][i].ray[k]) > 10. ) */
  /*        top[i][k]=(unsigned char)(quota_f((float)(vol_pol[l][i].teta_true)*CONV_RAD+0.45*DTOR,k)/100.);{ //top in ettometri */
  /*        if (l >= NEL -1 ) top[i][k]=0; */
  /*      } */
  /*    } */
  /*     } */
  /*   } */
  /*   return ier; */
  /* } */

/**
 * Write a banner with program information, compile flags, environment and so on.
 */
static void startup_banner()
{
    LOG_CATEGORY("radar.banner");

    /// Log initial program state
    LOG_INFO("%s -- Lancio Programma", PrendiOra());
    LOG_INFO("-----------------------------------------------------------------");
    LOG_INFO("Flag di Compilazione: "
#ifdef BOLOGNA
            " BOLOGNA"
#endif
#ifdef SPC
            " SPC"
#endif
#ifdef WRITE_DBP
            " WRITE_DBP"
#endif
#ifdef TIME
            " TIME"
#endif
#ifdef WRITE_DBP_REORDER
            " WRITE_DBP_REORDER"
#endif
#ifdef DECLUTTER
            " DECLUTTER"
#endif
#ifdef Z_AVERAGE
            " Z_AVERAGE"
#endif
#ifdef R_AVERAGE
            " R_AVERAGE"
#endif
#ifdef Z_LOWRIS
            " Z_LOWRIS"
#endif
#ifdef ANAPROP
            " ANAPROP"
#endif
#ifdef SHORT
            " SHORT"
#endif
#ifdef MEDIUM
            " MEDIUM"
#endif
#ifdef STATIC
            " STATIC"
#endif
#ifdef BEAMBLOCKING
            " BEAMBLOCKING"
#endif
#ifdef QUALITY
            " QUALITY"
#endif
            ".");

    LOG_INFO("-----------------------------------------------------------------");
    LOG_INFO("Variabili d'Ambiente:");
    LOG_INFO("LISTA_FILE = %s", getenv("LISTA_FILE"));
    LOG_INFO("LAST_FILE = %s", getenv("LAST_FILE"));
    LOG_INFO("ANAP_STAT_FILE = %s", getenv("ANAP_STAT_FILE"));
    LOG_INFO("BLOC_STAT_FILE = %s", getenv("BLOC_STAT_FILE"));
    LOG_INFO("ELEV_STAT_FILE = %s", getenv("ELEV_STAT_FILE"));
    LOG_INFO("DIR_OUT_PP_BLOC = %s", getenv("DIR_OUT_PP_BLOC"));
    LOG_INFO("FILE_DEM_SPC = %s", getenv("FILE_DEM_SPC"));
    LOG_INFO("FILE_DEM_GAT = %s", getenv("FILE_DEM_GAT"));
    LOG_INFO("DIR_QUALITY = %s", getenv("DIR_QUALITY"));
    LOG_INFO("FIRST_LEVEL_FILE = %s", getenv("FIRST_LEVEL_FILE"));
    LOG_INFO("OUTPUT_Z_DIR = %s", getenv("OUTPUT_Z_DIR"));
    LOG_INFO("OUTPUT_RAIN_DIR = %s", getenv("OUTPUT_RAIN_DIR"));
    LOG_INFO("OUTPUT_Z_LOWRIS_DIR = %s", getenv("OUTPUT_Z_LOWRIS_DIR"));
    LOG_INFO("-----------------------------------------------------------------");
}

/// Report a command line error and quit
void commandline_error(const char* progname, const char* msg)
{
    fprintf(stderr, "%s\n\n", msg);
    fprintf(stderr, "Usage: %s volume_file file_type site_name\n", progname);
    exit(1);
}

/* ================================ */
int main (int argc, char **argv)
    /* ================================ */
{
    char *nome_file;
    int i,k,l;//indici azimut ,range, elevazione
    int ier,ier_test, ier_main=0;//uscite errore generico (lettura volume e succ.anap) , di test_file, del main
    FILE *output;
    char *sito;//GAT O SPC
    float a,b,c;

    // Initialize logging
    Logging logging;

    LOG_CATEGORY("radar.main");

    //------- verifica n0 argomenti ------

    if (argc < 4)
        commandline_error(argv[0], "some argument is missing.");

    nome_file = argv[1];

    int file_type; // -- x definire n_elev e reolution e se è =1 esco
    if (sscanf(argv[2], "%d", &file_type) != 1)
        commandline_error(argv[0], "file_type must be a number");

    sito = argv[3];  //---- assegnazioni legate a argv[3]-- sito,  ambiente di lavoro, elevazioni
    setwork(argv[3]);  //-------setto ambiente lavoro (se var amb lavoro non settate le setta in automatico) ------

    startup_banner();

#ifdef TIME
    prendo_tempo();
#endif

    CUM_BAC *cb = new CUM_BAC;

    // Set feature flags
#ifdef QUALITY
    cb->do_quality = true;
#endif
#ifdef BEAMBLOCKING
    cb->do_beamblocking = true;
#endif
#ifdef DECLUTTER
    cb->do_declutter = true;
#endif

    if (cb->esegui_tutto(nome_file, file_type, sito))
        ier_main = 0;
    else
        ier_main = 1;
    delete cb;

    // è stato tolto il loop sui volumi

    LOG_INFO("End of processing, result: %d", ier_main);

    return ier_main;
}

