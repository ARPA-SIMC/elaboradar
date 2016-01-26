/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief definisce struttura Site Contiene le informazioni di base che caratterizzano il sito radar
*/

#ifndef ARCHIVIATORE_SITE_H
#define ARCHIVIATORE_SITE_H

#include <vector>
#include <string>

namespace elaboradar {

/**
 * Radar site information
 */
/*!
 * @brief Structure gestisce le informazioni sul sito
 * 
 */
struct Site
{
/*!
 * string_name 
 * @brief Nome sito radar
 */
    std::string name;
/*!
 * radar_lat 
 * @brief latitude of radar site
 */
    float radar_lat;
/*!
 * radar_lon 
 * @brief longitude of radar site
 */
    float radar_lon;
/*!
 * vpr_iaz_min
 * @brief azimuth index of the begin of the area for vpr computation 
 */
    int vpr_iaz_min;
/*!
 * vpr_iaz_max
 * @brief  azimuth index of the end of the area for vpr computation 
 */
    int vpr_iaz_max;

/*!
 * @brief Destructor 
 */ 	
    virtual ~Site();

/*!
 * @brief Return dem file name
 * @return dem file name [char *]
 */
    virtual const char* get_dem_file_name() const = 0;
/*!
 * @brief Return first_elev file name
 * @param [in] month - month  
 * @return first_elev file name [char *]
 */
    virtual const char* get_first_level_file_name(unsigned month) const = 0;
/*!
 *  @brief return the elev array used
 *  @param [in] medium -  flag to specify medium pulse
 *  @return values of elevation used 
 */
    virtual std::vector<double> get_elev_array(bool medium=false) const = 0;
/*! @brief Return the magic number for wind to be used in clean procedure
 *  @param [in] when - time information 
 *  @return wind maginc number [unsigned char]
 */
    virtual unsigned char get_bin_wind_magic_number(time_t when) const = 0;

/*!
 * @brief Get a Site object according to a site name.
 *
 * Currently supported are "GAT" and "SPC".
 *
 * Throws an exception in case of unsupported site name.
 *
 * @param [in] name - name of radar site
 * @return - Site Object [Site]
 */
   static const Site& get(const char* name);
};

}

#endif
