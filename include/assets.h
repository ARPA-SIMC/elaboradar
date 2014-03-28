#ifndef ARCHIVIATORE_ASSETS_H
#define ARCHIVIATORE_ASSETS_H

#include <string>
#include <ctime>
#include <cstdio>
#include <logging.h>
#include <H5Cpp.h>

namespace cumbac {

template<typename T> struct Matrix2D;

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
    mutable H5::H5File* outfile_devel_data;

public:
    Assets();
    ~Assets();

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
    void load_first_level(Matrix2D<unsigned char>& matrix);

    /**
     * Open the first level elevation BB el file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    void load_first_level_bb_el(Matrix2D<unsigned char>& matrix);

    /**
     * Open the first level elevation BB bloc file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     */
    void load_first_level_bb_bloc(Matrix2D<unsigned char>& matrix);

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
     *  fornisce temperatura al suolo, da lettura file esterno
     *  @brief   funzione che restituisce la temperatura al suolo
     *  @details  apre file temperature , legge lon lat e t, calcola differenze
     *  rispetto coordinate radar, se diff < soglia media il dato, stampa il nr
     *  di dati usati per la media e ritorna la temperatura
     *  @param[out] t_gr temperatura al suolo
     *  @return ierr codice di uscita (0=ok 1=fallito)
     */
    float read_t_ground() const;

    /**
     * Read the gap between the time in $LAST_VPR and the current acquisition time
     *
     *  @brief funzione che  calcola il no di quarti d'ora che intercorrono dall'ultimo profilo calcolato (combinato)  memorizzato in 'LAST_VPR'
     *  @param[in] nomefile nome del file LAST_VPR dove c'e' la data cui si riferisce l'ultimo profilo prodotto in n0 di secondi a partire da istante di riferimento
     *  @return gap1 ritorna il  no di quarti d'ora che intercorrono dall'ultimo profilo calcolato
     */
    long int read_profile_gap() const;

    /// Write the acquisition time in $LAST_VPR file
    void write_last_vpr();

    void write_vpr_hmax(int hvprmax);

    /**
     * Return an open HDF5 File to which we can write datasets used to debug
     * run information
     */
    H5::H5File get_devel_data_output() const;

protected:
    /// Compute the file name of a date/time based file in $DIR_OUT_PP_BLOC
    std::string fname_out_pp_bloc(const char* suffix) const;

    template<typename T>
    void load_raw(const std::string& fname, const char* desc, Matrix2D<T>& matrix);
};

}

#endif
