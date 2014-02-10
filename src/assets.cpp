#include "assets.h"
#include "utils.h"
#include "geo_par.h"
#include "vpr_par.h"
#include "site.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <stdint.h>
#include <stdexcept>

using namespace std;

namespace cumbac {

Assets::Assets()
    : logging_category(log4c_category_get("radar.assets"))
{
}

void Assets::configure(const char* site, time_t acq_time)
{
    configure(Site::get(site), acq_time);
}

void Assets::configure(const Site& site, time_t acq_time)
{
    conf_site = &site;
    conf_acq_time = acq_time;
    struct tm* tempo = gmtime(&acq_time);
    conf_year = tempo->tm_year + 1900;
    conf_month = tempo->tm_mon + 1;
    conf_day = tempo->tm_mday;
    conf_hour = tempo->tm_hour;
    conf_minute = tempo->tm_min;
}

FILE* Assets::open_file_dem()
{
    const char* fname = conf_site->get_dem_file_name();
    LOG_INFO("Opening dem file %s", fname);
    return fopen_checked(fname, "rt", "file dem");
}

FILE* Assets::open_file_first_level()
{
    const char* fname = getenv("FIRST_LEVEL_FILE");
    if (!fname)
        fname = conf_site->get_first_level_file_name(conf_month);
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
       LOG_ERROR("$FILE_T is not set");
       return NODATAVPR;
    }

    FILE* file_t = fopen(fname, "rt");
    if (!file_t)
    {
        LOG_ERROR("Cannot open $FILE_T=%s: %s", fname, strerror(errno));
        return NODATAVPR;
    }

    float media_t = 0;
    int icount = 0;
    float lon, lat, t;

    while (1) {
        if(fscanf(file_t,"%f %f %f \n",&lon,&lat,&t) == EOF) break;
        if (fabs(conf_site->radar_lat-lat)<=maxdlat && fabs(conf_site->radar_lon-lon)<=maxdlon) {
            ++icount;
            media_t += t - 273.15;
        }
    }

    fclose(file_t);

    if (icount == 0)
    {
        LOG_ERROR("Temperature data not found in $FILE_T=%s", fname);
        return NODATAVPR;
    }

    media_t /= (float)icount;
    LOG_INFO("ho %i stazioni dati affidabili e la t media è %f\n", icount, media_t);
    return media_t;
}

long int Assets::read_profile_gap()
{
    LOG_CATEGORY("radar.vpr");
    const char* fname = getenv("LAST_VPR");
    if (!fname)
    {
        LOG_ERROR("$LAST_VPR is not set");
        return 100;
    }

    FILE *file = fopen(fname, "rb");
    if (!file)
    {
        LOG_ERROR("Cannot open $LAST_VPR=%s: %s", fname, strerror(errno));
        return 100;
    }

    // FIXME: time_t può essere 64 bit, qui viene sempre troncato.
    // FIXME: l'ideale sarebbe, in questo caso, usare fprintf/fscanf invece di
    // FIXME: fread/fwrite
    uint32_t last_time;
    fread(&last_time, 4, 1, file);
    fclose(file);

    long int gap1 = abs(conf_acq_time - last_time)/900;
    LOG_INFO("old_data_header.norm.maq.acq_date last_time gap %d %ld %ld", conf_acq_time, last_time, gap1);

    return gap1;
}

void Assets::write_last_vpr()
{
    LOG_CATEGORY("radar.vpr");
    const char* fname = getenv("LAST_VPR");
    if (!fname) throw runtime_error("$LAST_VPR is not set");
    FILE* out = fopen_checked(fname, "wb", "ultimo VPR");
    uint32_t val = conf_acq_time;
    fwrite(&val, 4, 1, out);
    fclose(out);
}

}
