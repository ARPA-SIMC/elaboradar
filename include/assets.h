#ifndef ARCHIVIATORE_ASSETS_H
#define ARCHIVIATORE_ASSETS_H

#include <string>
#include <ctime>
#include <cstdio>
#include <logging.h>
#include <functional>
#include <vector>
#include <H5Cpp.h>

namespace elaboradar {

struct Config;
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
    const Config& cfg;
    const Site* conf_site;
    time_t conf_acq_time;
    int conf_year;
    int conf_month;
    int conf_day;
    int conf_hour;
    int conf_minute;
    mutable H5::H5File* outfile_devel_data;

public:
    Assets(const Config& cfg);
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
     * Save acq_time in $LAST_FILE, comparing it with the previous value.
     *
     * @returns
     *  false, if acq_time is older than the previous file processed
     *  true, if acq_time is newer, if $LAST_FILE does not exist or if
     *  $LAST_FILE is not set.
     */
    bool save_acq_time(time_t acq_time=0);

    /**
     * Open the dem file.
     *
     * The result is always a valid file: it throws an exception if something
     * goes wrong.
     *
     * TODO: cos'è il dem?
     */
    void load_dem(Matrix2D<float>& matrix);

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
     * Read the hray file, calling a callback on each parsed value.
     *
     * Returns the value of dtrs.
     */
    double read_file_hray(std::function<void (unsigned el, unsigned bin, double value)> on_sample);

    /**
     * Read the hray file, calling a callback on each parsed value.
     *
     * Returns the value of dtrs.
     */
    double read_file_hray_inf(std::function<void (unsigned el, unsigned bin, double value)> on_sample);

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

    /**
     * Read the value of $VPR_HEATING
     *
     * Returns 0 if the file does not exist or cannot be read
     */
    int read_vpr_heating() const;

    /**
     * Write a new value to $VPR_HEATING
     */
    void write_vpr_heating(int value) const;

    /**
     * Read $FILE_ZERO_TERMICO
     *
     * @returns true if the file was found and read correctly, false if zeroterm was not set
     */
    bool read_0term(float& zeroterm);

    /// Write the acquisition time in $LAST_VPR file
    void write_last_vpr();

    /// Read VPR_HMAX, returning -9999 if not found
    int read_vpr_hmax();

    void write_vpr_hmax(int hvprmax);

    bool read_vpr0(std::vector<float>& vpr0, std::vector<long int>& area);
    void write_vpr0(std::vector<float>& vpr, std::vector<long int>& area);

    /**
     * Return an open HDF5 File to which we can write datasets used to debug
     * run information
     */
    H5::H5File get_devel_data_output() const;

    /**
     * Write an image in a raw file in ${dir_env_var}, with the acquisition
     * date as file name and the given extension.
     *
     * desc is used to get better error messages.
     */
    void write_image(const Matrix2D<unsigned char>& image, const char* dir_env_var, const char* ext, const char* desc);

    template<typename T>
    void write_gdal_image(const Matrix2D<T>& image, const char* dir_env_var, const char* name, const char* format);

protected:
    /// Build a basename (without extension) for a file given the current
    /// acquisition time
    std::string fname_from_acq_time() const;

    /// Compute the file name of a date/time based file in $DIR_OUT_PP_BLOC
    std::string fname_out_pp_bloc(const char* suffix) const;

    /// Load a Matrix2D, from packed row-major binary data
    template<typename T>
    void load_raw(const std::string& fname, const char* desc, Matrix2D<T>& matrix);

    /// Load a Matrix2D, from space-separated column-major ascii floats
    void load_ascii(const std::string& fname, const char* desc, Matrix2D<float>& matrix);
};

}

#endif
