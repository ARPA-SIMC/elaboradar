#ifndef ARCHIVIATORE_SITE_H
#define ARCHIVIATORE_SITE_H

#include <vector>
#include <string>

namespace cumbac {

/**
 * Radar site information
 */
struct Site
{
    std::string name;
    float radar_lat;
    float radar_lon;
    int vpr_iaz_min;
    int vpr_iaz_max;

    virtual ~Site();

    virtual const char* get_dem_file_name() const = 0;
    virtual const char* get_first_level_file_name(unsigned month) const = 0;
    virtual std::vector<double> get_elev_array(bool medium=false) const = 0;
    virtual unsigned char get_bin_wind_magic_number(time_t when) const = 0;

    /**
     * Get a Site object according to a site name.
     *
     * Currently supported are "GAT" and "SPC".
     *
     * Throws an exception in case of unsupported site name.
     */
    static const Site& get(const char* name);
};

}

#endif
