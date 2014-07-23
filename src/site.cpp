#include "site.h"
#include "utils.h"
#include "geo_par.h"
#include "vpr_par.h"
#include <stdexcept>
#include <cstring>

using namespace std;

namespace {
vector<double> make_elev_array(const int* data, unsigned count)
{
    vector<double> res;
    res.reserve(count);
    for (unsigned i = 0; i < count; ++i)
        res.push_back(data[i] * 360. / 4096.);
    return res;
}
}

namespace elaboradar {

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
        vpr_iaz_min=IAZ_MIN_GAT;
        vpr_iaz_max=IAZ_MAX_GAT;
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

    virtual std::vector<double> get_elev_array(bool medium=false) const
    {
        if (medium)
        {
            static const int elev_data[]={6,16,26,37,47,57,80,109,148,205,284, 300, 305, 310, 315 };
            return make_elev_array(elev_data, sizeof(elev_data) / sizeof(int));
        } else {
            static const int elev_data[]={6,16,26,37,47,57,80,109,148,205,284, 300, 305, 310, 315 };
            return make_elev_array(elev_data, sizeof(elev_data) / sizeof(int));
        }
    }

    virtual unsigned char get_bin_wind_magic_number(time_t when) const
    {
        return 135;
    }
} site_gat;



struct SiteSPC : public Site
{
    SiteSPC()
    {
        name = "SPC";
        radar_lat=SPC_LAT;
        radar_lon=SPC_LON;
        vpr_iaz_min=IAZ_MIN_SPC;
        vpr_iaz_max=IAZ_MAX_SPC;
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

    virtual std::vector<double> get_elev_array(bool medium=false) const
    {
        if (medium)
        {
            static const int elev_data[]={6,16,26,36,47,57,80,108,148,205,284,300,305,310,315};
            return make_elev_array(elev_data, sizeof(elev_data) / sizeof(int));
        } else {
            static const int elev_data[]={6,16,26,36,47,57,80,108,148,205,284,300,305,310,315};
            return make_elev_array(elev_data, sizeof(elev_data) / sizeof(int));
        }
    }

    virtual unsigned char get_bin_wind_magic_number(time_t when) const
    {
        // After DBP2_250920131130_BOLOGNA
        if (when >= 1380108600)
          return 135;
        else
          return 131;
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
