// Guard condition: da mettere in test ad ogni .h
#ifndef ARCHIVIATORE_SETWORK_H
#define ARCHIVIATORE_SETWORK_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief settaggio ambiente lavoro nel caso non sia settato dall'esterno
 *  @details settaggio ambiente lavoro (variabili ambiente corrispondenti a files e directories di lavoro) nel caso non siano settate
*/

/**
 *  @brief funzione che setta ambiente lavoro nel caso non sia settato dall'esterno
 *  @param sito (gat o spc)
 *  @return ier_code of the last putenv call (????)
*/
int setwork(char *sito);

/**
 *  @brief funzione che ripulisce tutto l'ambiente lavoro
*/
void unsetwork();

/**
 *  @brief funzione che stampa l'ambiente lavoro con chiamata a LOG_INFO
*/
void printwork();

#ifdef __cplusplus
}
#endif

#endif // ARCHIVIATORE_SETWORK_H
