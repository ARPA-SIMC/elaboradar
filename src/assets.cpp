#include "assets.h"
#include "config.h"
#include <radarelab/utils.h>
#include <radarelab/matrix.h>
#include <radarelab/image.h>
#include <radarelab/algo/dbz.h>
#include <radarelab/algo/vpr.h>
#include <radarelab/vpr_par.h>
#include "site.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <stdint.h>
#include <stdexcept>

using namespace std;
using namespace radarelab;

namespace elaboradar {

Assets::Assets(const Config& cfg)
    : logging_category(log4c_category_get("radar.assets")), cfg(cfg), outfile_devel_data(0)
{
#ifdef HAVE_GDAL
    gdal_init_once();
#endif
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
    if (acq_time == 0) acq_time = conf_acq_time;

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
    if (!fname){
        LOG_WARN("leggo datipath da conf_site..");
        fname = conf_site->get_first_level_file_name(conf_month);
    }
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

namespace {

double parse_hray(File& fd, std::function<void (unsigned el, unsigned bin, double value)> on_sample)
{
    size_t line_no = 0;
    double dtrs;
    fd.read_lines([&](char* line, size_t len) {
        if (line_no == 0)
        {
            // Read dtrs in the first line
            dtrs = strtod(line, NULL);
        } else {
            char* s = line;
            int el = 0;
            while (true)
            {
                char* next;
                double val = strtod(s, &next);
                if (next == s) break;
                on_sample(el, line_no - 1, val);
                s = next;
                ++el;
            }
        }
        ++line_no;
    });
    if (line_no == 0)
        throw std::runtime_error("hray/hray_inf file is empty");
    return dtrs;
}

}

double Assets::read_file_hray(std::function<void (unsigned el, unsigned bin, double value)> on_sample)
{
    File fd(logging_category);
    if (!fd.open(fname_out_pp_bloc("h_ray.txt"), "rb", "hray"))
        throw std::runtime_error("cannot open hray file");
    return parse_hray(fd, on_sample);
}

double Assets::read_file_hray_inf(std::function<void (unsigned el, unsigned bin, double value)> on_sample)
{
    File fd(logging_category);
    if (!fd.open(fname_out_pp_bloc("h_rayinf.txt"), "rb", "hray inf"))
        throw std::runtime_error("cannot open hray inf file");
    return parse_hray(fd, on_sample);
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
    File in(logging_category);
    if (!in.open_from_env("FILE_T", "rt"))
        return NODATAVPR;

    float media_t = 0;
    int icount = 0;
    float lon, lat, t;

    while (1) {
        if(fscanf(in, "%f %f %f \n",&lon,&lat,&t) == EOF) break;
        if (fabs(conf_site->radarSite.lat_r-lat)<=maxdlat && fabs(conf_site->radarSite.lon_r-lon)<=maxdlon) {
            ++icount;
            media_t += t - 273.15;
        }
    }

    if (icount == 0)
    {
        LOG_ERROR("Temperature data not found in $FILE_T=%s", in.name());
        return NODATAVPR;
    }

    media_t /= (float)icount;
    LOG_INFO("ho %i stazioni dati affidabili e la t media è %f\n", icount, media_t);
    return media_t;
}

long int Assets::read_profile_gap() const
{
    LOG_CATEGORY("radar.vpr");
    File in(logging_category);
    if (!in.open_from_env("LAST_VPR", "rb"))
        return 100;

    // FIXME: time_t può essere 64 bit, qui viene sempre troncato.
    // FIXME: l'ideale sarebbe, in questo caso, usare fprintf/fscanf invece di
    // FIXME: fread/fwrite
    uint32_t last_time;
    fread(&last_time, 4, 1, in);

    long int gap1 = abs(conf_acq_time - last_time)/900;
    LOG_INFO("old_data_header.norm.maq.acq_date last_time gap %ld %u %ld", conf_acq_time, last_time, gap1);

    return gap1;
}

int Assets::read_vpr_heating() const
{
    LOG_CATEGORY("radar.vpr");
    File in(logging_category);
    if (!in.open_from_env("VPR_HEATING", "rt"))
        return 0;

    int heating;
    if (fscanf(in, "%i ", &heating) != 1)
    {
        LOG_ERROR("Cannot read $VPR_HEATING=%s: %s", in.name(), strerror(errno));
        return 0;
    }

    return heating;
}

void Assets::write_vpr_heating(int value) const
{
    LOG_CATEGORY("radar.vpr");
    File out(logging_category);
    if (!out.open_from_env("VPR_HEATING", "wt"))
        return;

    if (fprintf(out, " %i \n", value) < 0)
        LOG_ERROR("Cannot write $VPR_HEATING=%s: %s", out.name(), strerror(errno));
}

bool Assets::read_0term(float& zeroterm)
{
    LOG_CATEGORY("radar.class");
    File in(logging_category);
    if (!in.open_from_env("FILE_ZERO_TERMICO", "rt"))
        return false;

    if (fscanf(in, "%f", &zeroterm) != 1)
    {
        LOG_ERROR("$FILE_ZERO_TERMICO=%s cannot be read: %s", in.name(), strerror(errno));
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
    File in(logging_category);
    if (!in.open_from_env("VPR_HMAX", "rt"))
        return -9999;

    int value;
    if (fscanf(in, "%i", &value) != 1)
    {
        LOG_ERROR("$VPR_HMAX=%s cannot be read: %s", in.name(), strerror(errno));
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

bool Assets::read_vpr0(algo::VPR& vpr0)
{
    File in(logging_category);
    if (!in.open_from_env("VPR0_FILE", "rt")) return false;

    for (unsigned i = 0; i < vpr0.size(); ++i)
        //-----leggo vpr e area per ogni strato----
        if (fscanf(in, "%f %li\n", &vpr0.val[i], &vpr0.area[i]) != 2)
        {
            LOG_ERROR("$VPR0_FILE=%s cannot be read: %s", in.name(), strerror(errno));
            throw std::runtime_error("cannot read $VPR0_FILE");
        }

    return true;
}

bool Assets::read_archived_vpr(const algo::DBZ& dbz, time_t time, radarelab::algo::VPR& vpr)
{
    const char* dir = getenv("DIR_STORE_VPR"); //--questa non sarebbe una dir_arch più che store?... contesto il nome...
    if (!dir) return false;

    struct tm t;
    gmtime_r(&time, &t);

    char fname[64];
    snprintf(fname, 64, "%04d%02d%02d%02d%02d_vpr_%s",
            t.tm_year + 1900, t.tm_mon + 1, t.tm_mday,
            t.tm_hour, t.tm_min, conf_site->name.c_str());

    string pathname = dir;
    pathname += "/";
    pathname += fname;

    File in(logging_category);
    if (!in.open(pathname, "r", "archived VPR file"))
        return false;

    // TODO: check the actual format of the file and make the parsing safe:
    // currently if one of these strings is longer than 99, we crash or worse.
    char stringa[100];
    fscanf(in, " %s %s %s %s" ,stringa ,stringa,stringa,stringa);
    for (unsigned ilay=0; ilay < vpr.size(); ++ilay){
        float vpr_dbz;
        long int ar;
        int il;
        fscanf(in, " %i %f %li", &il, &vpr_dbz, &ar);  //---NB il file in archivio è in dBZ e contiene anche la quota----

        //---- converto in R il profilo vecchio--
        if (vpr_dbz > 0)
        {
            vpr.val[ilay] = dbz.DBZtoR(vpr_dbz);
            vpr.area[ilay] = ar;
        }
        else
            vpr.val[ilay] = NODATAVPR;
    }

    return true;
}

bool Assets::find_vpr0(const radarelab::algo::DBZ& dbz, radarelab::algo::VPR& vpr0, long int& gap)
{
    /*------calcolo la distanza temporale che separa l'ultimo profilo calcolato dall'istante attuale--*/
    /* (dentro il file LAST_VPR c'è una data che contiene la data cui si riferisce il vpr in n0 di secondi dall'istante di riferimento)*/
    gap = read_profile_gap();

    /*------leggo il profilo vecchio più recente di MEMORY ----*/
    /*------nota bene: è in R ovvero  pioggia!! ----*/

    if (!read_vpr0(vpr0))
    {
        LOG_WARN("non esiste file vpr vecchio: %s",getenv("VPR0_FILE"));

        //----se file non esiste assegno gap=100----
        gap = 100;
    }

    //------------se gap < MEMORY leggo vpr e area per ogni strato-----------

    if (gap <= MEMORY)
        return true;

    //-----Se gap > MEMORY

    //a)----- tento .. sono in POST-ELABORAZIONE:----

    //-----devo andare a ricercare tra i profili 'buoni' in archivio quello con cui combinare il dato----
    //---- trattandosi di profili con data nel nome del file, costruisco il nome a partire dall'istante corrente ciclando su un numero di quarti d'ora
    //---- pari a memory finchè non trovo un profilo. se non lo trovo gap resta=100

    /* questo per fare ciclo sul vpr vecchio*/
    time_t Time = conf_acq_time;

    // TODO: cerca in archivio se esiste un VPR piú recente del
    // last_vpr: togliere dal calcolo VPR generico e spostarlo nel
    // punto dove viene caricato il VPR precedente
    for (unsigned i=0;i<MEMORY;i++)
        if (read_archived_vpr(dbz, Time + i * 900, vpr0))
        {
            gap = 0;
            return true;
        }

    return false;
}

void Assets::write_vpr0(const algo::VPR& vpr)
{
    const char* fname = getenv("VPR0_FILE");
    if (!fname) throw runtime_error("$VPR0_FILE (ultimo vpr) is not set");
    FILE* out = fopen_checked(fname, "wt", "ultimo vpr");
    for (unsigned i = 0; i < vpr.size(); ++i)
        if (fprintf(out, " %10.3f %li\n", vpr.val[i], vpr.area[i]) < 0)
        {
            LOG_ERROR("$VPR0_FILE=%s cannot be written: %s", fname, strerror(errno));
            fclose(out);
            throw std::runtime_error("cannot write to $VPR0_FILE");
        }
    fclose(out);
}

void Assets::write_dbz_coefficients(const algo::DBZ& dbz)
{
    const char* dirname = getenv("OUTPUT_Z_LOWRIS_DIR");
    if (!dirname) throw runtime_error("OUTPUT_Z_LOWRIS_DIR is not set");
    string fname(dirname);
    fname += "/MP_coeff";
    File out(logging_category);
    out.open(fname, "wb", "MP coefficients");

    unsigned char MP_coeff[2]; /* a/10 e b*10 per scrivere come 2 byte */
    MP_coeff[0]=(unsigned char)(dbz.aMP/10);
    MP_coeff[1]=(unsigned char)(dbz.bMP*10);

    fwrite(MP_coeff, sizeof(MP_coeff), 1, out);
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
    if ((unsigned)fsize != matrix.size() * sizeof(T))
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
    struct tm *tempo = gmtime(&conf_acq_time);
    strftime(buf, buf_size, "%Y%m%d%H%M", tempo);
    return buf;
}

void Assets::write_image(const Matrix2D<unsigned char>& image, const char* dir_env_var, const char* ext, const char* desc)
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

    // Convert to south-north columns scanned west to east
    Matrix2D<unsigned char> transformed(image.cols(), image.rows());
    for (unsigned y = 0; y < image.cols(); ++y)
        for (unsigned x = 0; x < image.rows(); ++x)
            transformed(x, image.cols()-1-y) = image(y, x);
    if (fwrite(transformed.data(), transformed.size(), 1, out) != 1)
    {
        LOG_WARN("cannot write to %s: %s", fname.c_str(), strerror(errno));
        fclose(out);
        throw std::runtime_error("cannot write to image file");
    }

    fclose(out);
}

void Assets::write_subimage(const Matrix2D<unsigned char>& image,unsigned image_side, const char* dir_env_var, const char* ext, const char* desc)
{
    const char* dir = getenv(dir_env_var);
    if (!dir)
    {
        LOG_INFO("$%s not set", dir_env_var);
        throw runtime_error("required env var is not set");
    }

    string fname = string(dir) + "/" + fname_from_acq_time() + "_" + std::to_string(image_side) + ext;
    FILE* out = fopen_checked(fname.c_str(), "wb", desc);

    LOG_INFO("aperto file %s dimensione matrice %zd\n", fname.c_str(), image.size());

    // Convert to south-north columns scanned west to east
    unsigned xofs = (image.cols() - image_side) / 2;
    unsigned yofs = (image.rows() - image_side) / 2;
    //LOG_INFO(" Image_size %4d , Image.cols %4d Image.Rows %4d -- xofs %d yofs %d", image_side, image.cols(), image.rows(), xofs, yofs);
    Matrix2D<unsigned char> transformed(image_side, image_side);
    for (unsigned y = 0; y < image_side; ++y)
        for (unsigned x = 0; x < image_side; ++x)
            transformed(x, image_side-1-y) = image(y + yofs, x + xofs);

    if (fwrite(transformed.data(), transformed.size(), 1, out) != 1)
    {
        LOG_WARN("cannot write to %s: %s", fname.c_str(), strerror(errno));
        fclose(out);
        throw std::runtime_error("cannot write to image file");
    }

    fclose(out);
}

void Assets::write_subimage(const Matrix2D<unsigned char>& image, unsigned image_side, string algos,const char* dir_env_var, const char* ext, const char* desc)
{
    const char* dir = getenv(dir_env_var);
    if (!dir)
    {
        LOG_INFO("$%s not set", dir_env_var);
        throw runtime_error("required env var is not set");
    }
    string fname = string(dir) + "/" + fname_from_acq_time() + "_" + std::to_string(image_side) + "_"+algos+ext;
    FILE* out = fopen_checked(fname.c_str(), "wb", desc);

    LOG_INFO("aperto file %s dimensione matrice %zd\n", fname.c_str(), image.size());

    // Convert to south-north columns scanned west to east
    unsigned xofs = (image.cols() - image_side) / 2;
    unsigned yofs = (image.rows() - image_side) / 2;
    LOG_INFO(" Image_size %4d , Image.cols %4d Image.Rows %4d -- xofs %d yofs %d", image_side, (int)image.cols(), (int)image.rows(), xofs, yofs);
    Matrix2D<unsigned char> transformed(image_side, image_side);
    for (unsigned y = 0; y < image_side; ++y)
        for (unsigned x = 0; x < image_side; ++x)
            transformed(x, image_side-1-y) = image(y + yofs, x + xofs);

    if (fwrite(transformed.data(), transformed.size(), 1, out) != 1)
    {
        LOG_WARN("cannot write to %s: %s", fname.c_str(), strerror(errno));
        fclose(out);
        throw std::runtime_error("cannot write to image file");
    }

    fclose(out);
}

template<typename T>
void Assets::write_gdal_image(const Matrix2D<T>& image, const char* dir_env_var, const char* name, const char* format)
{
#ifdef HAVE_GDAL
    const char* dir = getenv(dir_env_var);
    if (!dir)
    {
        LOG_INFO("$%s not set", dir_env_var);
        throw runtime_error("required env var is not set");
    }

    string fname = string(dir) + "/" + fname_from_acq_time() + "-" + name + "." + gdal_extension_for_format(format);

    radarelab::write_image(image, fname, format);
#else
    throw std::runtime_error("GDAL support was not enabled at compile time");
#endif
}

template void Assets::write_gdal_image(const Matrix2D<unsigned char>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const Matrix2D<unsigned short>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const Matrix2D<int>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const Matrix2D<unsigned>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const Matrix2D<short>&, const char*, const char*, const char*);
template void Assets::write_gdal_image(const Matrix2D<double>&, const char*, const char*, const char*);

   time_t Assets::getAcqTime() { return conf_acq_time;} ;
   RadarSite Assets::getRadarSite() { return conf_site->radarSite; };

}
