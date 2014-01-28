#include "assets.h"
#include "utils.h"
#include "geo_par.h"
#include "vpr_par.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <stdexcept>

using namespace std;

Assets::Assets()
    : logging_category(log4c_category_get("radar.assets"))
{
}

void Assets::configure(const char* sito, time_t acq_time)
{
    if (strcmp(sito, "GAT") == 0)
        conf_site = SITE_GAT;
    else if (strcmp(sito, "SPC") == 0)
        conf_site = SITE_SPC;
    else
    {
        string errmsg(sito);
        throw domain_error(errmsg + " is not a valid radar site name");
    }

    struct tm* tempo = gmtime(&acq_time);
    conf_year = tempo->tm_year + 1900;
    conf_month = tempo->tm_mon + 1;
    conf_day = tempo->tm_mday;
    conf_hour = tempo->tm_hour;
    conf_minute = tempo->tm_min;
}

FILE* Assets::open_file_dem()
{
    const char* fname;
    switch (conf_site)
    {
        case SITE_SPC: fname = getenv_default("FILE_DEM_SPC", "../../PP+BLOC/dati/dem_SanPi.txt"); break;
        case SITE_GAT: fname = getenv_default("FILE_DEM_GAT", "../../PP+BLOC/dati/dem_Gatta.txt"); break;
    }
    LOG_INFO("Opening dem file %s", fname);
    return fopen_checked(fname, "rt", "file dem");
}

FILE* Assets::open_file_first_level()
{
    const char* fname = getenv("FIRST_LEVEL_FILE");
    if (!fname)
    {
        switch (conf_site)
        {
            case SITE_SPC:
                if (1 <= conf_month && conf_month <= 3)
                    fname = "../dati/FIRST_LEVEL_SPC_2006_INV";
                else if (4 <= conf_month && conf_month <= 9)
                    fname = "../dati/FIRST_LEVEL_SPC_2006_PRI-EST";
                else
                    fname = "../dati/FIRST_LEVEL_SPC_2006_AUT";
                break;
            case SITE_GAT:
                if (1 <= conf_month && conf_month <= 3)
                    fname = "../dati/FIRST_LEVEL_GAT_2006_INV";
                else if (4 <= conf_month && conf_month <= 9)
                    fname = "../dati/FIRST_LEVEL_GAT_2006_PRI-EST";
                else
                    fname = "../dati/FIRST_LEVEL_GAT_2006_AUT";
                break;
        }
    }
    LOG_INFO("Opening mappa statica %s", fname);
    return fopen_checked(fname, "rb", "mappa statica");
}

int Assets::read_file_first_level_dim()
{
    const char* fname = getenv("FIRST_LEVEL_DIM_FILE");
    if (!fname) throw runtime_error("FIRST_LEVEL_DIM_FILE is not set");
    FILE* fd = fopen_checked(fname, "rt", "dimensioni mappa statica");

    // We can directly read its contents: it's just one number in ASCII
    int dim;
    int res = fscanf(fd, "%i", &dim);
    if (res == EOF)
    {
        string errmsg("Error reading ");
        errmsg += fname;
        errmsg += ": ";
        errmsg += strerror(errno);
        fclose(fd);
        throw runtime_error(errmsg);
    }

    if (res != 1)
    {
        string errmsg("Error reading ");
        errmsg += fname;
        errmsg += ": the file does not seem to contain a number";
        errmsg += strerror(errno);
        fclose(fd);
        throw runtime_error(errmsg);
    }

    fclose(fd);

    return dim;
}

FILE* Assets::open_file_first_level_bb_el()
{
    string fname = fname_out_pp_bloc("mat_el.bin");
    LOG_INFO("Opening elev BB %s", fname.c_str());
    return fopen_checked(fname.c_str(), "rb", "elev BB");
}

FILE* Assets::open_file_first_level_bb_bloc()
{
    string fname = fname_out_pp_bloc("mat_bloc.bin");
    LOG_INFO("Opening elev BB %s", fname.c_str());
    return fopen_checked(fname.c_str(), "rb", "elev BB");
}

FILE* Assets::open_file_hray()
{
    string fname = fname_out_pp_bloc("h_ray.txt");
    LOG_INFO("Opening hray %s", fname.c_str());
    return fopen_checked(fname.c_str(), "rb", "hray");
}

FILE* Assets::open_file_hray_inf()
{
    string fname = fname_out_pp_bloc("h_rayinf.txt");
    LOG_INFO("Opening hray inf %s", fname.c_str());
    return fopen_checked(fname.c_str(), "rb", "hray inf");
}

std::string Assets::fname_out_pp_bloc(const char* suffix) const
{
    const char* dir = getenv("DIR_OUT_PP_BLOC");
    if (!dir) throw runtime_error("DIR_OUT_PP_BLOC is not set");

    char fname[1024];
    sprintf(fname, "%s/%04d%02d%02d%02d%02d%s", dir,
            conf_year, conf_month, conf_day, conf_hour, conf_minute, suffix);

    return fname;
}

float Assets::read_t_ground()
{
    LOG_CATEGORY("radar.vpr");
    const char* fname = getenv("FILE_T");
    if (!fname)
    {
       LOG_ERROR("FILE_T is not set");
       return NODATAVPR;
    }

    FILE* file_t = fopen(fname, "rt");
    if (!file_t)
    {
        LOG_ERROR("Cannot open FILE_T=%s: %s", fname, strerror(errno));
        return NODATAVPR;
    }

    float media_t = 0;
    int icount = 0;
    float radar_lat, radar_lon, lon, lat, t;

    switch (conf_site)
    {
        case SITE_GAT:
            radar_lat=GAT_LAT;
            radar_lon=GAT_LON;
            break;
        case SITE_SPC:
            radar_lat=SPC_LAT;
            radar_lon=SPC_LON;
            break;
    }

    while (1) {
        if(fscanf(file_t,"%f %f %f \n",&lon,&lat,&t) == EOF) break;
        if (fabs(radar_lat-lat)<=maxdlat && fabs(radar_lon-lon)<=maxdlon) {
            ++icount;
            media_t += t - 273.15;
        }
    }

    fclose(file_t);

    if (icount == 0)
    {
        LOG_ERROR("Temperature data not found in FILE_T=%s", fname);
        return NODATAVPR;
    }

    media_t /= (float)icount;
    LOG_INFO("ho %i stazioni dati affidabili e la t media è %f\n", icount, media_t);
    return media_t;
}
