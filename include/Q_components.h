#ifndef ARCHIVIATORE_QCOMPONENTS_H
#define ARCHIVIATORE_QCOMPONENTS_H
/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief funzioni componenti di qualita' radar semplici
 *  @details  funzioni componenti di qualita' radar semplici, una per ogni fattore di contaminazione
 i parametri per calcolo qualita' sono in qual_par.h
*/
/**
 *  @brief funzione componente qualita' beam blocking
 *  @param[in] bbc beam blocking
 *  @param[in] dr distanza da radiosondaggio (m)
 *  @param[in] dt tempo intercorso da radiosondaggio (h)
 *  @return valore componente qualita' beam blocking
 * 
*/

float qBB(unsigned char bbc,float dr,float dt);

/**
 *  @brief funzione componente qualita' clutter
 *  @details  funzione qualita' relativa a clutter che dipende dal valore di ritorno del controllo anaprop:  0 = dato ok    qCl=1, 1 = anaprop    qCl=0.5,  2=  no data    qCl=0,  3=  no control qCl=0.8
 *  @param[in] clc parametro anaprop: ANAP_OK=0 , ANAP_YES=1 (presenza anaprop), ANAP_NODAT=2 (no data), ANAP_NOCONTROL=3 (no controllo) 
 *  @return valore componente qualita' clutter
 * 
*/
float qCl( unsigned char clc);

/**
 *  @brief funzione componente qualita' distanza
 *  @param[in] dr distanza in metri da radiosondaggio 
 *  @return  valore componente qualita' distanza
*/
float qDist(float dr);

/**
 *  @brief funzione componente qualita' focalizzazione fascio
 *  @param[in] dh altezza cella da propagazione secondo radiosondaggio (metri)
 *  @param[in] dhst altezza cella standard (metri)
 *  @return valore componente qualita'da focalizzazione fascio
*/
float qVol( float dh, float dhst);

/**
 *  @brief funzione componente qualita' da path integrated attenuation
 *  @param[in] PIA path integrated attenuation (dbZ)
 *  @return valore componente qualita' da attenuazione
*/
float qAtt( float PIA);

/*
 *  @brief funzione componente qualita' da VPR
 *  @param dZ correzione
 *  @param stdev standard deviation della correzione
*/
float qVpr(float dZ,float sdevZ);

#endif //ARCHIVIATORE_QCOMPONENTS_H
