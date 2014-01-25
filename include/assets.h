#ifndef ARCHIVIATORE_ASSETS_H
#define ARCHIVIATORE_ASSETS_H

#include <string>
#include <ctime>
#include <cstdio>
#include <logging.h>

/**
 * Finds resources, like data files, used by the program.
 */
class Assets
{
    // Una volta che abbiamo spostato tutta la ricerca dei file di dati in
    // questa classe, dovremmo avere un modo per sapere esattamente di cosa ha
    // bisogno il programma, e per cambiarlo in caso cambino le esigenze
    // operative.

protected:
    log4c_category_t* logging_category;
    enum {
        SITE_GAT,
        SITE_SPC,
    } conf_site;
    int conf_year;
    int conf_month;
    int conf_day;

public:
    Assets();

    /**
     * Configure asset lookup with the given details.
     *
     * sito: "GAT" or "SPC" according to the radar site.
     * time: the volume acquisition time
     */
    void configure(const char* sito, time_t acq_time);

    /**
     * Open the dem file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     *
     * TODO: cos'è il dem?
     */
    FILE* open_file_dem();

    /**
     * Open the first level file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    FILE* open_file_first_level();
};

#endif