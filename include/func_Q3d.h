// Guard condition: da mettere in test ad ogni .h
#ifndef ARCHIVIATORE_FUNC_Q3D_H
#define ARCHIVIATORE_FUNC_Q3D_H


/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief funzioni che combinano le componenti semplici di qualita' radar
 *  @details   funzioni qualita' radar rispetto a Z e rispetto a R che combinano le componenti semplici di qualita' radar;
 i parametri per calcolo qualita' sono in qual_par.h
*/
/**
*  @brief funzione che calcola la qualita' per Z
*  @param[in] cl clutter da anaprop
*  @param[in] bb  beam blocking in %
*  @param[in] dst distanza da radar (metri)
*  @param[in] dr distanza da radiosondaggio (metri)
*  @param[in] dt intervallo tracorso da radiosondaggio (h)
*  @param[in] dh altezza del bin reale (m)
*  @param[in] dhst altezza del bin standard (m)
*  @param[in] PIA path integrated attenuation (dB)
*  @return valore qualita' finale
*/
float func_q_Z(unsigned char cl,unsigned char bb,float dst,float dr,float dt,float dh,float dhst,float PIA);

#endif //ARCHIVIATORE_FUNC_Q3D_H
