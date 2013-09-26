
#include <func_Q3d.h>

#include <stdio.h> 
#include <math.h> //file di sistema definisce costanti matematiche e richiama funzioni mat
#include <qual_par.h> //file creato da me, parametri qualita' e definizione funzioni qualita'

float func_q_Z(cl,bb,dst,dr,dt,dh,dhst,PIA)

unsigned char  bb;
unsigned char  cl;
float          PIA,dr,dt,dh,dhst,dst;

{
float q;

 q=qCl(cl)*qBB(bb,dr,dt)*qDist(dst)*qVol(dh,dhst)*qAtt(PIA);

      return q;
}
/**
*  @brief funzione che calcola la qualita' per R
*  @param[in] cl clutter da anaprop
*  @param[in] bb  beam blocking in %
*  @param[in] dst distanza da radar (metri)
*  @param[in] dr distanza da radiosondaggio (metri)
*  @param[in] dt intervallo tracorso da radiosondaggio (h)
*  @param[in] dh altezza del bin reale (m)
*  @param[in] dhst altezza del bin standard (m)
*  @param[in] PIA path integrated attenuation (dB)
*  @param[in] dZ correzione vpr (dB)
*  @param[in] sdevZ standard deviation di Z (dB)
*  @return valore qualita' finale
*/
/*
 omstart func_q_R
  idx calcola qualita' rispetto a prodotto precipitazione
  calcola qualita' rispetto a prodotto precipitazione
  q=qCl(cl)*qBB(bb,dr,dt)*qDist(dr)*qVol(dr,dh)*qAtt(PIA)*qVpr(dZ,sdevZ)
 omend
*/
float func_q_R(cl,bb,dst,dr,dt,dh,dhst,PIA,dZ,sdevZ)

unsigned char  bb;
unsigned char  cl;
float          PIA,dr,dt,dh,dst,dZ,dhst,sdevZ;
{
float q;

 q=qCl(cl)*qBB(bb,dr,dt)*qDist(dst)*qVol(dh,dhst)*qAtt(PIA)*qVpr(dZ,sdevZ);
      return q;
}

