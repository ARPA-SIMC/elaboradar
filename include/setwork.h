// Guard condition: da mettere in test ad ogni .h
#ifndef ARCHIVIATORE_SETWORK_H
#define ARCHIVIATORE_SETWORK_H

/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief settaggio ambiente lavoro nel caso non sia settato dall'esterno
 *  @details settaggio ambiente lavoro (variabili ambiente corrispondenti a files e directories di lavoro) nel caso non siano settate
*/

/**
 *  @brief funzione che setta ambiente lavoro nel caso non sia settato dall'esterno
 *  @param sito (gat o spc)
*/
int setwork(char *sito);


#endif // ARCHIVIATORE_SETWORK_H
