#ifndef ARCHIVIATORE_ASSETS_H
#define ARCHIVIATORE_ASSETS_H

#include <string>
#include <ctime>
#include <cstdio>
#include <logging.h>

namespace cumbac {

struct Site;

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
    const Site* conf_site;
    time_t conf_acq_time;
    int conf_year;
    int conf_month;
    int conf_day;
    int conf_hour;
    int conf_minute;

public:
    Assets();

    /**
     * Configure asset lookup with the given details.
     *
     * sito: Site object for the radar site.
     * time: the volume acquisition time
     */
    void configure(const Site& site, time_t acq_time);

    /**
     * Configure asset lookup with the given details.
     *
     * sito: "GAT" or "SPC" according to the radar site.
     * time: the volume acquisition time
     */
    void configure(const char* site, time_t acq_time);

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

    /**
     * Read the value inside the first level dimension file
     */
    int read_file_first_level_dim();

    /**
     * Open the first level elevation BB el file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    FILE* open_file_first_level_bb_el();

    /**
     * Open the first level elevation BB bloc file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    FILE* open_file_first_level_bb_bloc();

    /**
     * Open the hray file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    FILE* open_file_hray();

    /**
     * Open the hray inf file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    FILE* open_file_hray_inf();

    /**
     * Load ground temperature
     */
    float read_t_ground() const;

    /// Read the gap in seconds between the time in $LAST_VPR and the current acquisition time
    long int read_profile_gap() const;

    /// Write the acquisition time in $LAST_VPR file
    void write_last_vpr();

protected:
    /// Compute the file name of a date/time based file in $DIR_OUT_PP_BLOC
    std::string fname_out_pp_bloc(const char* suffix) const;
};

}

#endif
