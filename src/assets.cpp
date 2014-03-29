#include "assets.h"
#include "utils.h"
#include "geo_par.h"
#include "vpr_par.h"
#include "matrix.h"
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
    : logging_category(log4c_category_get("radar.assets")), outfile_devel_data(0)
{
}

Assets::~Assets()
{
    if (outfile_devel_data)
        delete outfile_devel_data;
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

bool Assets::save_acq_time(time_t acq_time)
{
    // If LAST_FILE is not set, return true
    const char* last_file = getenv("LAST_FILE");
    if (last_file == NULL)
    {
        LOG_INFO("$LAST_FILE not set");
        return true;
    }

    bool res = true;
    uint32_t last_time;

    FILE* fp = fopen(last_file, "r");

    // If the file does not exist, return true
    if (fp == NULL)
    {
        LOG_INFO("$LAST_FILE=%s does not exist", last_file);
        last_time = 0;
        goto check;
    }

    // If the file is empty, return true
    if (fread(&last_time, 4, 1, fp) != 1)
    {
        LOG_INFO("$LAST_FILE=%s cannot be read", last_file);
        last_time = 0;
        goto check;
    }

check:
    {
        int diff = acq_time - last_time;
        LOG_INFO("%s: new acq_time is old %c %d", last_file, diff < 0 ? '-' : '+', abs(diff));
    }

    if (acq_time <= last_time)
        res = false;

close:
    if (fp) fclose(fp);

update:
    if ((fp = fopen(last_file, "w")) == NULL)
    {
        LOG_WARN("cannot write to %s: %s", last_file, strerror(errno));
        throw std::runtime_error("cannot (re)create $LAST_FILE");
    }

    // Fit acq_time in 4 bytes (FIXME: so far so good, until 2036)
    last_time = acq_time;
    if (fwrite(&last_time, 4, 1, fp) != 1)
    {
        LOG_WARN("cannot write to %s: %s", last_file, strerror(errno));
        throw std::runtime_error("cannot write to $LAST_FILE");
    }
    fclose(fp);

    return res;
}

FILE* Assets::open_file_dem()
{
    const char* fname = conf_site->get_dem_file_name();
    LOG_INFO("Opening dem file %s", fname);
    return fopen_checked(fname, "rt", "file dem");
}

void Assets::load_first_level(Matrix2D<unsigned char>& matrix)
{
    const char* fname = getenv("FIRST_LEVEL_FILE");
    if (!fname)
        fname = conf_site->get_first_level_file_name(conf_month);
    load_raw(fname, "mappa statica", matrix);
}

void Assets::load_first_level_bb_el(Matrix2D<unsigned char>& matrix)
{
    load_raw(fname_out_pp_bloc("mat_el.bin"), "elev BB", matrix);
}

void Assets::load_first_level_bb_bloc(Matrix2D<unsigned char>& matrix)
{
    load_raw(fname_out_pp_bloc("mat_bloc.bin"), "elev BB", matrix);
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

float Assets::read_t_ground() const
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

long int Assets::read_profile_gap() const
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
    //LOG_CATEGORY("radar.vpr");
    const char* fname = getenv("LAST_VPR");
    if (!fname) throw runtime_error("$LAST_VPR is not set");
    FILE* out = fopen_checked(fname, "wb", "ultimo VPR");
    uint32_t val = conf_acq_time;
    fwrite(&val, 4, 1, out);
    fclose(out);
}

void Assets::write_vpr_hmax(int hvprmax)
{
    const char* fname = getenv("VPR_HMAX");
    if (!fname) throw runtime_error("$VPR_HMAX is not set");
    FILE* out = fopen_checked(fname, "wt", "hmax VPR");
    fprintf(out, "%d", hvprmax);
    fclose(out);
}

H5::H5File Assets::get_devel_data_output() const
{
    if (!outfile_devel_data)
    {
        const char* qdir = getenv("DIR_QUALITY");
        if (!qdir) throw runtime_error("$DIR_QUALITY is not set");
        string fname(qdir);
        fname += "/devel-data.h5";
        outfile_devel_data = new H5::H5File(fname, H5F_ACC_TRUNC);
    }
    return *outfile_devel_data;
}

template<class T>
void Assets::load_raw(const std::string& fname, const char* desc, Matrix2D<T>& matrix)
{
    LOG_CATEGORY("radar.io");
    LOG_INFO("Opening %s %s", desc, fname.c_str());
    FILE* in = fopen_checked(fname.c_str(), "rb", desc);

    // Read the file size
    fseek(in, 0,SEEK_END);
    long fsize = ftell(in);
    rewind(in);

    // Check that the file size is consistent with what we want
    if (fsize != matrix.size() * sizeof(T))
    {
        LOG_ERROR("Il file %s è %ld byte ma dovrebbe invece essere %ld byte\n",
                fsize, matrix.size() * sizeof(T));
        throw std::runtime_error("La dimensione della mappa statica non è quello che mi aspetto");
    }
    LOG_INFO ("DIMENSIONE MAPPA STATICA %u %u", matrix.SY, matrix.SX);

    for (unsigned i = 0; i < matrix.SY; ++i)
        if (fread(matrix.data + i * matrix.SX, matrix.SX, 1, in) != 1)
        {
            std::string errmsg("Error reading ");
            errmsg += fname;
            errmsg += ": ";
            errmsg += strerror(errno);
            fclose(in);
            throw std::runtime_error(errmsg);
        }

    fclose(in);
}

}
