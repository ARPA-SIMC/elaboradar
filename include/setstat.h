// Guard condition: da mettere in test ad ogni .h
#ifndef ARCHIVIATORE_SETSTAT_H
#define ARCHIVIATORE_SETSTAT_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief settaggio ambiente statico
*/

/**
 *  @brief funzione che setta le variabili d'ambiente
 *  @details  setta ambiente statico usando le getenv e  in caso in cui le variabili d'ambiente non siano definite all'esterno le definisce automaticamente
 *  @param sito sito (gat o spc)
 *  @param mese mese in numero ordinale 
 *  @param n_dem nome dem
 *  @param n_fl nome first level
*/
int setstat(const char *sito, int mese, char *n_dem, char *n_fl);

#ifdef __cplusplus
}
#endif

#endif // ARCHIVIATORE_SETSTAT_H
