#ifndef ARCHIVIATORE_SITE_H
#define ARCHIVIATORE_SITE_H

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

    virtual ~Site();

    virtual const char* get_dem_file_name() const = 0;
    virtual const char* get_first_level_file_name(unsigned month) const = 0;
    virtual void fill_elev_array(int* elev_array, bool medium=false) const = 0;

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