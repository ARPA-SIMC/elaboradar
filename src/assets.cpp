#include "assets.h"
#include "utils.h"
#include "geo_par.h"
#include "vpr_par.h"
#include "matrix.h"
#include "site.h"
#include "image.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <stdint.h>
#include <stdexcept>

using namespace std;

namespace {

struct InputFile
{
    log4c_category_t* logging_category;
    const char* varname;
    const char* desc;
    const char* fname;
    FILE* file;

    InputFile(log4c_category_t* logging_category, const char* varname, const char* desc=0)
        : logging_category(logging_category), varname(varname), desc(desc ? desc : varname), fname(0), file(0)
    {
    }
    ~InputFile()
    {
        if (file) fclose(file);
    }

    bool open(const char* mode)
    {
        fname = getenv(varname);
        if (!fname)
        {
            LOG_ERROR("$%s is not set", varname);
            return false;
        }

        file = fopen(fname, mode);
        if (!file)
        {
            LOG_ERROR("Cannot open $%s=%s: %s", varname, fname, strerror(errno));
            return false;
        }
        return true;
    }
};

}


namespace cumbac {

Assets::Assets()
    : logging_category(log4c_category_get("radar.assets")), outfile_devel_data(0)
{
    gdal_init_once();
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

    if (fp) fclose(fp);

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

void Assets::load_dem(Matrix2D<float>& matrix)
{
    load_ascii(conf_site->get_dem_file_name(), "file dem", matrix);
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
    InputFile in(logging_category, "FILE_T");
    if (!in.open("rt"))
        return NODATAVPR;

    float media_t = 0;
    int icount = 0;
    float lon, lat, t;

    while (1) {
        if(fscanf(in.file, "%f %f %f \n",&lon,&lat,&t) == EOF) break;
        if (fabs(conf_site->radar_lat-lat)<=maxdlat && fabs(conf_site->radar_lon-lon)<=maxdlon) {
            ++icount;
            media_t += t - 273.15;
        }
    }

    if (icount == 0)
    {
        LOG_ERROR("Temperature data not found in $FILE_T=%s", in.fname);
        return NODATAVPR;
    }

    media_t /= (float)icount;
    LOG_INFO("ho %i stazioni dati affidabili e la t media è %f\n", icount, media_t);
    return media_t;
}

long int Assets::read_profile_gap() const
{
    LOG_CATEGORY("radar.vpr");
    InputFile in(logging_category, "LAST_VPR");
    if (!in.open("rb"))
        return 100;

    // FIXME: time_t può essere 64 bit, qui viene sempre troncato.
    // FIXME: l'ideale sarebbe, in questo caso, usare fprintf/fscanf invece di
    // FIXME: fread/fwrite
    uint32_t last_time;
    fread(&last_time, 4, 1, in.file);

    long int gap1 = abs(conf_acq_time - last_time)/900;
    LOG_INFO("old_data_header.norm.maq.acq_date last_time gap %ld %u %ld", conf_acq_time, last_time, gap1);

    return gap1;
}

int Assets::read_vpr_heating() const
{
    LOG_CATEGORY("radar.vpr");
    InputFile in(logging_category, "VPR_HEATING");
    if (!in.open("rt"))
        return 0;

    int heating;
    if (fscanf(in.file, "%i ", &heating) != 1)
    {
        LOG_ERROR("Cannot read $VPR_HEATING=%s: %s", in.fname, strerror(errno));
        return 0;
    }

    return heating;
}

void Assets::write_vpr_heating(int value) const
{
    LOG_CATEGORY("radar.vpr");
    InputFile out(logging_category, "VPR_HEATING");
    if (!out.open("wt"))
        return;

    if (fprintf(out.file, " %i \n", value) < 0);
        LOG_ERROR("Cannot write $VPR_HEATING=%s: %s", out.fname, strerror(errno));
}

bool Assets::read_0term(float& zeroterm)
{
    LOG_CATEGORY("radar.class");
    InputFile in(logging_category, "FILE_ZERO_TERMICO");
    if (!in.open("rt"))
        return false;

    if (fscanf(in.file, "%f", &zeroterm) != 1)
    {
        LOG_ERROR("$FILE_ZERO_TERMICO=%s cannot be read: %s", in.fname, strerror(errno));
        return false;
    }

    return true;
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

int Assets::read_vpr_hmax()
{
    InputFile in(logging_category, "VPR_HMAX");
    if (!in.open("rt"))
        return -9999;

    int value;
    if (fscanf(in.file, "%i", &value) != 1)
    {
        LOG_ERROR("$VPR_HMAX=%s cannot be read: %s", in.fname, strerror(errno));
        return -9999;
    }

    LOG_INFO("fatta lettura hmax vpr = %i", value);
    return value;
}

void Assets::write_vpr_hmax(int hvprmax)
{
    const char* fname = getenv("VPR_HMAX");
    if (!fname) throw runtime_error("$VPR_HMAX is not set");
    FILE* out = fopen_checked(fname, "wt", "hmax VPR");
    fprintf(out, "%d", hvprmax);
    fclose(out);
}

bool Assets::read_vpr0(std::vector<float>& vpr0, std::vector<long int>& area)
{
    InputFile in(logging_category, "VPR0_FILE");
    if (!in.open("rt")) return false;

    for (unsigned i = 0; i < vpr0.size(); ++i)
        //-----leggo vpr e area per ogni strato----
        if (fscanf(in.file, "%f %li\n", &vpr0[i], &area[i]) != 2)
        {
            LOG_ERROR("$VPR0_FILE=%s cannot be read: %s", in.fname, strerror(errno));
            throw std::runtime_error("cannot read $VPR0_FILE");
        }

    return true;
}

void Assets::write_vpr0(std::vector<float>& vpr, std::vector<long int>& area)
{
    const char* fname = getenv("VPR0_FILE");
    if (!fname) throw runtime_error("$VPR0_FILE (ultimo vpr) is not set");
    FILE* out = fopen_checked(fname, "wt", "ultimo vpr");
    for (unsigned i = 0; i < vpr.size(); ++i)
        if (fprintf(out, " %10.3f %li\n", vpr[i], area[i]) < 0)
        {
            LOG_ERROR("$VPR0_FILE=%s cannot be written: %s", fname, strerror(errno));
            fclose(out);
            throw std::runtime_error("cannot write to $VPR0_FILE");
        }
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
                fname.c_str(), fsize, matrix.size() * sizeof(T));
        throw std::runtime_error("La dimensione della mappa statica non è quello che mi aspetto");
    }
    LOG_INFO ("DIMENSIONE MAPPA STATICA %ld %ld", matrix.rows(), matrix.cols());

    for (unsigned i = 0; i < matrix.rows(); ++i)
        if (fread(matrix.data() + i * matrix.cols(), matrix.cols(), 1, in) != 1)
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

void Assets::load_ascii(const std::string& fname, const char* desc, Matrix2D<float>& matrix)
{
    LOG_INFO("Opening %s %s", desc, fname.c_str());
    FILE* in = fopen_checked(fname.c_str(), "rt", desc);

    for (unsigned x = 0; x < matrix.cols(); ++x)
        for (unsigned y = 0; y < matrix.rows(); ++y)
        {
            float val;
            fscanf(in, "%f ", &val);
            matrix(y, x) = val;
        }

    fclose(in);
}

std::string Assets::fname_from_acq_time() const
{
    const unsigned buf_size = 16;
    char buf[buf_size];

    /*
    if( do_medium){
    	tempo = gmtime(&load_info.acq_date);
    	time = NormalizzoData(load_info.acq_date);
    	tempo = gmtime(&time);
    } else {
    	time = NormalizzoData(load_info.acq_date);
    	tempo = gmtime(&time);
    }
    */

    struct tm *tempo = gmtime(&conf_acq_time);

    snprintf(buf, buf_size, "%04d%02d%02d%02d%02d",
            tempo->tm_year+1900, tempo->tm_mon+1, tempo->tm_mday,
            tempo->tm_hour, tempo->tm_min);

    return buf;
}

void Assets::write_image(const cumbac::Matrix2D<unsigned char>& image, const char* dir_env_var, const char* ext, const char* desc)
{
    const char* dir = getenv(dir_env_var);
    if (!dir)
    {
        LOG_INFO("$%s not set", dir_env_var);
        throw runtime_error("required env var is not set");
    }

    string fname = string(dir) + "/" + fname_from_acq_time() + ext;
    FILE* out = fopen_checked(fname.c_str(), "wb", desc);

    LOG_INFO("aperto file %s dimensione matrice %zd\n", fname.c_str(), image.size());

    if (fwrite(image.data(), image.size(), 1, out) != 1)
    {
        LOG_WARN("cannot write to %s: %s", fname.c_str(), strerror(errno));
        fclose(out);
        throw std::runtime_error("cannot write to image file");
    }

    fclose(out);
}

template<typename T>
void Assets::write_gdal_image(const cumbac::Matrix2D<T>& image, const char* dir_env_var, const char* name, const char* format)
{
    const char* dir = getenv(dir_env_var);
    if (!dir)
    {
        LOG_INFO("$%s not set", dir_env_var);
        throw runtime_error("required env var is not set");
    }

    string fname = string(dir) + "/" + fname_from_acq_time() + "-" + name + "." + gdal_extension_for_format(format);

    cumbac::write_image(image, fname, format);
}

template void Assets::write_gdal_image(const cumbac::Matrix2D<unsigned char>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const cumbac::Matrix2D<unsigned short>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const cumbac::Matrix2D<double>&, const char*, const char*, const char*);


}
