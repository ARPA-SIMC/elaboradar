#include "site.h"
#include "utils.h"
#include "geo_par.h"
#include <stdexcept>
#include <cstring>

using namespace std;

namespace cumbac {

Site::~Site()
{
}

struct SiteGAT : public Site
{
    SiteGAT()
    {
        name = "GAT";
        radar_lat=GAT_LAT;
        radar_lon=GAT_LON;
    }

    virtual const char* get_dem_file_name() const
    {
        return getenv_default("FILE_DEM_GAT", "../../PP+BLOC/dati/dem_Gatta.txt");
    }

    virtual const char* get_first_level_file_name(unsigned month) const
    {
        if (1 <= month && month <= 3)
            return "../dati/FIRST_LEVEL_GAT_2006_INV";
        else if (4 <= month && month <= 9)
            return "../dati/FIRST_LEVEL_GAT_2006_PRI-EST";
        else
            return "../dati/FIRST_LEVEL_GAT_2006_AUT";
    }
} site_gat;

struct SiteSPC : public Site
{
    SiteSPC()
    {
        name = "SPC";
        radar_lat=SPC_LAT;
        radar_lon=SPC_LON;
    }

    virtual const char* get_dem_file_name() const
    {
        return getenv_default("FILE_DEM_SPC", "../../PP+BLOC/dati/dem_SanPi.txt");
    }

    virtual const char* get_first_level_file_name(unsigned month) const
    {
        if (1 <= month && month <= 3)
            return "../dati/FIRST_LEVEL_SPC_2006_INV";
        else if (4 <= month && month <= 9)
            return "../dati/FIRST_LEVEL_SPC_2006_PRI-EST";
        else
            return "../dati/FIRST_LEVEL_SPC_2006_AUT";
    }
} site_spc;

const Site& Site::get(const char* name)
{
    if (strcmp(name, "GAT") == 0)
        return site_gat;
    else if (strcmp(name, "SPC") == 0)
        return site_spc;
    else
    {
        string errmsg(name);
        throw domain_error(errmsg + " is not a valid radar site name");
    }
}

}
