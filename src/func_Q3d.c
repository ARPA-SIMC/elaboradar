/*
 omstart func_Q3d.c 
 idx funzioni qualita' radar 
 funzioni qualita' radar rispetto a Z e rispetto a R
 i parametri per calcolo qualita' sono in qual_par.h
 omend
*/
#include <stdio.h> 
#include <math.h> //file di sistema definisce costanti matematiche e richiama funzioni mat
#include <qual_par.h> //file creato da me, parametri qualita' e definizione funzioni qualita'
/*
 omstart func_q_Z
  idx calcola qualita' rispetto a prodotto riflettivita'
  calcola qualita' rispetto a prodotto riflettivita'
  q=qCl(cl)*qBB(bb,dr,dt)*qDist(dr)*qVol(dr,dh)*qAtt(PIA)
 omend
*/

float func_q_Z(cl,bb,dr,dt,dh,dhst,PIA)

unsigned char  bb;
unsigned char  cl;
float          PIA,dr,dt,dh,dhst;

{
float q;

 q=qCl(cl)*qBB(bb,dr,dt)*qDist(dr)*qVol(dr,dh,dhst)*qAtt(PIA);
 //printf("qCl(cl) %f qBB(bb,dr,dt) %f qDist(dr) %f qVol(dr,dh,dhst) %f qAtt(PIA) %f \n",qCl(cl),qBB(bb,dr,dt),qDist(dr),qVol(dr,dh,dhst),qAtt(PIA));
      return q;
}

/*
 omstart func_q_R
  idx calcola qualita' rispetto a prodotto precipitazione
  calcola qualita' rispetto a prodotto precipitazione
  q=qCl(cl)*qBB(bb,dr,dt)*qDist(dr)*qVol(dr,dh)*qAtt(PIA)*qVpr(dZ,sdevZ)
 omend
*/
float func_q_R(cl,bb,dr,dt,dh,dhst,PIA,dZ,sdevZ)

unsigned char  bb;
unsigned char  cl;
float          PIA,dr,dt,dh,dZ,sdevZ;
{
float q;

 q=qCl(cl)*qBB(bb,dr,dt)*qDist(dr)*qVol(dr,dh,dhst)*qAtt(PIA)*qVpr(dZ,sdevZ);
      return q;
}

