
/*
 omstart Q_components.c 
 idx libreria componenti qualita'
 contiene funzioni di qualita' relative a ciascun  fattore di errore 
 rispetto a calcolaqual.c sostituisco il file  bb_par.h con qual_par.h.

 q= 1.-(1.-qd)*(1.-qc)

 qd= qualita dato
 qc=qualita' correzione

 se errore non corretto qc=0

 omend
*/

#include <Q_components.h>
#include <elaboradar/algo/anaprop.h>
#include <stdlib.h>
#include <math.h> //file di sistema definisce costanti matematiche e richiama funzioni mat
#include <qual_par.h> //file creato da me, parametri qualita'


/*--------------------------------------*/
  /* funzione qualita' da  BB */
/*--------------------------------------*/
/*
 omstart qBB 
 idx funzione qualita' relativa a beam blocking

 qd=1.-(pow((float)(bb)/(float)(BBMAX),1/1.5))
 qc=fcBB*fcDt*fcDr*fcErrpt
 
 float qd;  qual. dato non corretto
 float qc;  qual.  correzione
 float fcBB;  comp. qual. corr. da BB
 float fcDt; comp. qual. corr. da distanza temporale radiosondaggio 
 float fcDr;  comp. qual. corr. da distanza spaziale radiosondaggio 
 float fcErrpt;  comp. qual. corr. da errore puntamento antenna
 omend
*/

   float qBB(unsigned char bbc,float dr,float dt)

   //unsigned char  bbc;
   //float          dr,dt;   /* distanza spaziale radiosondaggio*/

{

float qd; /* qual. dato non corretto*/
float qc; /* qual.  correzione*/
float fcBB; /* comp. qual. corr. da BB*/
float fcDt; /* comp. qual. corr. da distanza temporale radiosondaggio */
float fcDr; /* comp. qual. corr. da distanza spaziale radiosondaggio */
float fcErrpt; /* comp. qual. corr. da errore puntamento antenna*/

   if (bbc >  BBMAX) bbc=BBMAX;
   qd=1.-(pow((float)(bbc)/(float)(BBMAX),1/1.5)); 
  
   fcBB=1.-(pow((float)(bbc)/(float)(BBMAX),1/1.5));
   if (dt >  DTMAX) dt=DTMAX;
   fcDt=exp(-(double)(dt)/DTLIM);
   fcDr=exp(-(double)(dr)/DRLIM);
   fcErrpt=1.-pow(ERRPT,1/1.5);
   
   qc=fcBB*fcDt*fcDr*fcErrpt;
   
   return (1.-(1.-qd)*(1.-qc));
   // return (qd); // caso non corretto; 
}

/*--------------------------------------*/
  /* funzione qualita' da  Clutter */
/*--------------------------------------*/

/*
 omstart qCl
 idx funzione qualita' relativa a clutter
 funzione qualita' relativa a clutter
 dipende dal valore di ritorno del controllo anaprop

   0 = dato ok    qCl=1
   1 = anaprop    qCl=0.5
   2=  no data    qCl=0
   3=  no control qCl=0.8
 omend
*/

float qCl( unsigned char  clc)
{
    using namespace elaboradar::algo;
    switch(clc)
    {
        case ANAP_OK:
            return (1.);     //  ok
        case ANAP_YES:
            return (0.5);    // anaprop
        case ANAP_NODAT:
            return (0.);     // no  data
        case ANAP_NOCONTROL:
            return (0.8);   //  no control
    }
   
 }
/*--------------------------------------*/
  /* funzione qualita' da  Distanza (per ora ricavata da Koistinen e Puhakka con parametri Robby) */
/*--------------------------------------*/
/*
 omstart qDist
 idx funzione qualita' relativa a distanza
 qd=(exp(-(float)(BETA*dr/1000.)))
 omend
*/

   float qDist( float dr)
   //  float dr;

{
   
  return (exp(-(float)(BETA*dr/1000.)));

}

/*--------------------------------------*/
  /* funzione qualita' da errore volumetrico in anaprop */
/*--------------------------------------*/
/*
 omstart qVol
 idx funzione qualita' da errore volumetrico
 funzione qualita' da errore volumetrico, approcio geometrico-ottico
 qd=1.-sqrt(pow((1.-pow(dh/dhst,1/1.5)),2.))
 omend
*/

float qVol( float dh, float dhst)

//     float dh,dhst;

{


 
  return (1.-sqrt(pow((1.-pow(dh/dhst,1/1.5)),2.)));

}

/*--------------------------------------*/
  /* funzione qualita' da attenuazione */
/*--------------------------------------*/
/*
 omstart qAtt
 idx funzione qualita' da attenuazione
 funzione qualita' da attenuazione
 qd=pow(10,-PIA/15.) 
 omend
*/

  float qAtt(float PIA)

  //float PIA;

{
  
  return (pow(10,-PIA/15.));

}
/*--------------------------------------*/
  /* funzione qualita' da errore vpr*/
/*--------------------------------------*/
/*
 omstart qVpr
 idx funzione qualita' da errore Vpr
 funzione qualita' da errore Vpr

 qd=pow(10.,-(abs(dZ)+sdevZ)/15.0)
 qc= 1.-pow(10.,-(abs(dZ)-2.0*sdevZ));
    if (qc < 0) qc=0.0;
 
 float  dZ    correzione su Z
 float  sdevZ deviazione standard della correzione su Z
 omend
*/


/*
 *  @brief funzione componente qualita' da VPR
 *  @param dZ correzione
 *  @param stdev standard deviation della correzione
*/
  float qVpr(float  dZ,float sdevZ)

  //float dZ,sdevZ;

{
  float qd,qc,eps;
  eps=0.001;
  qd= pow(10.,-(abs(dZ)+sdevZ)/15.0)  ;
  if (abs(dZ) < eps) qd=1.0;

  qc= 1.-pow(10.,-(abs(dZ)-2.0*sdevZ));
  if (qc < 0) qc=0.0;
    
  return (1.- (1.-qd)*(1.-qc));
}

