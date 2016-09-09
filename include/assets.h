/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief Gestisce risorse usate dal programma.
*/
#ifndef ARCHIVIATORE_ASSETS_H
#define ARCHIVIATORE_ASSETS_H

#include <string>
#include <ctime>
#include <cstdio>
#include <radarelab/logging.h>
#include <functional>
#include <vector>
#include <H5Cpp.h>
#include <radarelab/RadarSite.h>

namespace radarelab {
template<typename T> struct Matrix2D;

namespace algo {
class DBZ;
class VPR;
}

}

namespace elaboradar {

struct Config;

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
/**
 * \brief  Constructor
 * \param [in] cfg - objet used to pass configuration
 *
 */

    Assets(const Config& cfg);
    ~Assets();

    /**
     * Configure asset lookup with the given details.
     *
     * @param [in] site -  Site object for the radar site.
     * @param [in] acq_time - Volume acquisition time
     */
    void configure(const Site& site, time_t acq_time);

    /**
     * Configure asset lookup with the given details.
     *
     * @param [in] site - "GAT" or "SPC" according to the radar site.
     * @param [in] acq_time - Volume acquisition time
     */
    void configure(const char* site, time_t acq_time);

    /**
     * Save acq_time in $LAST_FILE, comparing it with the previous value.
     *
     * @param [in] acq_time - Volume acquisition time
     * @return false, if acq_time is older than the previous file processed
     * @return true, if acq_time is newer, if $LAST_FILE does not exist or if $LAST_FILE is not set.
     */
    bool save_acq_time(time_t acq_time=0);

    /**
     * Open the dem file.
     *
     * The result is always a valid file: it throws an exception if something goes wrong.
     *
     * @param matrix - Matrix  [float] where dem is loaded
     *
     * TODO: cos'è il dem?
     */
    void load_dem(radarelab::Matrix2D<float>& matrix);

    /**
     * Open the first level file.
     *
     * The result is always a valid file: it throws an exception if something goes wrong.
     * @param matrix - Matrix  [unsigned char] where first_elev table  is loaded
     */
    void load_first_level(radarelab::Matrix2D<unsigned char>& matrix);

    /**
     * Open the first level elevation BB el file.
     *
     * The result is always a valid file: it throws an exception if something goes wrong.
     * @param matrix - Matrix  [unsigned char] where first_elev_bb_el table  is loaded
     */
    void load_first_level_bb_el(radarelab::Matrix2D<unsigned char>& matrix);

    /**
     * Open the first level elevation BB bloc file.
     *
     * The result is always a valid file: it throws an exception if something goes wrong.
     * @param matrix - Matrix  [unsigned char] where first_elev_bb_bloc table  is loaded
     */
    void load_first_level_bb_bloc(radarelab::Matrix2D<unsigned char>& matrix);

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
     *  @return it_gr temperatura al suolo
     */
    float read_t_ground() const;
     //*  @param[out] t_gr temperatura al suolo
     //*  @return ierr codice di uscita (0=ok 1=fallito)

    /**
     * Read the gap between the time in $LAST_VPR and the current acquisition time
     *
     *  @brief funzione che  calcola il numero di quarti d'ora che intercorrono dall'ultimo profilo calcolato (combinato)  memorizzato in 'LAST_VPR'
     *  @return gap1 ritorna il numero di quarti d'ora che intercorrono dall'ultimo profilo calcolato
     */
    long int read_profile_gap() const;
     //**  param[in] nomefile nome del file LAST_VPR dove c'e' la data cui si riferisce l'ultimo profilo prodotto in n0 di secondi a partire da istante di riferimento

    /**
     * Read the value of $VPR_HEATING (counter of consecutive vpr calculated, see scientific documentation)
     *
     * @return 0 if the file does not exist or cannot be read
     * @return The value of $VPR_HEATING
     */
    int read_vpr_heating() const;

    /**
     * Write a new value to $VPR_HEATING (counter of consecutive vpr calculated, see scientific documentation)
     * @param [in] value - value to be written
     */
    void write_vpr_heating(int value) const;

    /**
     * Read $FILE_ZERO_TERMICO
     *
     * @param [out] zeroterm  - height of 0°C level [m]
     * @return true if the file was found and read correctly, false if zeroterm was not set
     */
    bool read_0term(float& zeroterm);

    /** Write the acquisition time in $LAST_VPR file
     */
    void write_last_vpr();

    /** Read in $VPR_HMAX the vpr peak's height.
     * @return value or -9999 if not found [m]
     */  
    int read_vpr_hmax();

    /** write in $VPR_HMAX the vpr peak's height.
     * @param [in] hvprmax -  value [m]
     */  
    void write_vpr_hmax(int hvprmax);

    /** Read in $VPR0_FILE the last vpr available.
     * @param [out] vpr0 - vpr profile in mmh-1 [rain intensity]
     * @param [out] area - areal coverage for each layer km^2/1000
     * @return true if succesfull
     * @return false if file does not exits 
     */  
    bool read_vpr0(radarelab::algo::VPR& vpr0);

    /** Write in $VPR0_FILE the vpr calculated.
     * @param [in] vpr - vpr profile in mmh-1 [rain intensity]
     * @param [in] area - areal coverage for each layer km^2/1000
     */  
    void write_vpr0(const radarelab::algo::VPR& vpr);

    /** Write in $OUTPUT_Z_LOWRIS_DIR/MP_coeff the MP coefficients.
     * @param [in] dbz - DBZ object with MP coefficients 
     */  
    void write_dbz_coefficients(const radarelab::algo::DBZ& dbz);

    /**
     * Return an open HDF5 File ( $DIR_QUALITY/devel-data.h5) to which we can write datasets used to debug
     * run information
     * @ return H5FILE object
     */
    H5::H5File get_devel_data_output() const;

    /**
     * Write an image in a raw file in ${dir_env_var}, with the acquisition
     * date as file name and the given extension.
     * Image matrix is tranformed in  out_image(x,image.cols()-1-y) = image(y, x);
     *
     * @param [in] image - Matrix2D to be written 
     * @param [in] dir_env_var - file path
     * @param [in] ext - file extension
     * @param [out] desc - used to get better error messages.
     */
    void write_image(const radarelab::Matrix2D<unsigned char>& image, const char* dir_env_var, const char* ext, const char* desc);

    /**
     * Write an image in a raw file in ${dir_env_var}, with the acquisition
     * date as file name and the given extension.
     * Image matrix is tranformed in  out_image(x,image.cols()-1-y) = image(y, x);
     *
     * @param [in] image - Matrix2D to be written 
     * @param [in] image_side - write only the square central part of the original image, with this size in pixels of the image side
     * @param [in] dir_env_var - file path
     * @param [in] ext - file extension
     * @param [out] desc - used to get better error messages.
     */
    void write_subimage(const radarelab::Matrix2D<unsigned char>& image, unsigned image_side, const char* dir_env_var, const char* ext, const char* desc);
    /**
     * Write an image in a raw file in ${dir_env_var}, with the acquisition
     * date as file name and the given extension.
     * Image matrix is tranformed in  out_image(x,image.cols()-1-y) = image(y, x);
     *
     * @param [in] image - Matrix2D to be written 
     * @param [in] image_side - write only the square central part of the original image, with this size in pixels of the image side
     * @param [in] algos - Algorithms applied during processing  (i.e. BB_VPR_CLASS)
     * @param [in] dir_env_var - file path
     * @param [in] ext - file extension
     * @param [out] desc - used to get better error messages.
     */
    void write_subimage(const radarelab::Matrix2D<unsigned char>& image, unsigned image_side, std::string algos, const char* dir_env_var, const char* ext, const char* desc);

    /**
     * Write a graphic image with gdal.
     *
     * @param [in] image - Matrix2D to be written 
     * @param [in] dir_env_var - file path
     * @param [in] name - file name 
     * @param [in] format  - file graphic format used.
     */
    template<typename T>
    void write_gdal_image(const radarelab::Matrix2D<T>& image, const char* dir_env_var, const char* name, const char* format);
    
    /// Build a basename (without extension) for a file given the current
    /// acquisition time
    std::string fname_from_acq_time() const;

    time_t getAcqTime();
    RadarSite getRadarSite(); 


protected:

    /// Compute the file name of a date/time based file in $DIR_OUT_PP_BLOC
    std::string fname_out_pp_bloc(const char* suffix) const;

    /// Load a Matrix2D, from packed row-major binary data
    template<typename T>
    void load_raw(const std::string& fname, const char* desc, radarelab::Matrix2D<T>& matrix);

    /// Load a Matrix2D, from space-separated column-major ascii floats
    void load_ascii(const std::string& fname, const char* desc, radarelab::Matrix2D<float>& matrix);
};

}

#endif
