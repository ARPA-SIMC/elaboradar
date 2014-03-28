#include "volume.h"
#include "utils.h"
#include "site.h"
#include "volume_cleaner.h"
#include "logging.h"
#include <radarlib/radar.hpp>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <cstring>
#include <cerrno>
#include <H5Cpp.h>

//#define IMPRECISE_AZIMUT

#ifdef __cplusplus
extern "C" {
#endif
// libreria radar
#include <func_SP20read.h>
#ifdef __cplusplus
}
#endif

#ifdef NEL
#undef NEL
#endif

using namespace std;

/// This is not used anymore, but it is here to satisfy libSP20 linking needs
int elev_array[15];

namespace cumbac {

void LoadLog::print(FILE* out) const
{
    if (empty())
    {
        fprintf(out, "no beams loaded\n");
        return;
    }

    for (const_iterator i = begin(); i != end(); ++i)
    {
        if (i != begin())
            fprintf(out, ", ");
        fprintf(out, "ϑ%.2f α%.2f", i->theta, i->alpha);
    }
    fprintf(out, "\n");
}

PolarScan::PolarScan(unsigned beam_size)
    : beam_count(NUM_AZ_X_PPI), beam_size(beam_size), elevation(0)
{
    if (beam_size > 0)
    {
        scan = gsl_matrix_alloc(NUM_AZ_X_PPI, beam_size);
        gsl_matrix_set_all(scan, BYTEtoDB(1));
    } else {
        scan = 0;
    }

    beam_info.resize(beam_count);
}

PolarScan::~PolarScan()
{
    if (scan)
        gsl_matrix_free(scan);
}

unsigned char PolarScan::get_raw(unsigned az, unsigned beam) const
{
    return DBtoBYTE(get_db(az, beam));
}

void PolarScan::set_raw(unsigned az, unsigned beam, unsigned char val)
{
    gsl_matrix_set(scan, az, beam, BYTEtoDB(val));
}

unsigned PolarScan::count_rays_filled() const
{
    unsigned count = 0;
    for (vector<BeamInfo>::const_iterator i = beam_info.begin(); i != beam_info.end(); ++i)
        if (!i->load_log.empty())
            ++count;
    return count;
}

void PolarScan::read_beam_db(unsigned az, float* out, unsigned out_size, float missing) const
{
    // Prima riempio il minimo tra ray.size() e out_size
    size_t set_count = min(beam_size, out_size);

    for (unsigned i = 0; i < set_count; ++i)
        out[i] = get_db(az, i);

    for (unsigned i = set_count; i < out_size; ++i)
        out[i] = missing;
}

const LoadLog& PolarScan::get_beam_load_log(unsigned az) const
{
    return beam_info[az].load_log;
}

void PolarScan::fill_beam(int el_num, double theta, double alpha, unsigned size, const double* data)
{
    int alfa = alpha / FATT_MOLT_AZ;
    if (alfa >= 4096) return;

    int az_num = azimut_index_MDB(alfa);

    /*
    if (az_num == 0)
    {
        printf("fbeam ϑ%f→%d α%f→%d %u", theta, el_num, alpha, az_num, size);
        for (unsigned i = 0; i < 20; ++i)
            printf(" %d", (int)data[i]);
        printf("\n");
    }
    */

    merge_beam(el_num, az_num, theta, alpha, size, data);
    if(az_num*0.9 - alpha < 0.)
    {
        int new_az_num = (az_num + 1) % 400;
        if (new_az_num != az_num)
            merge_beam(el_num, new_az_num, theta, alpha, size, data);
    }
    else if(az_num*0.9 - alpha > 0.)
    {
        int new_az_num = (az_num -1+400) %400;
        if (new_az_num != az_num)
            merge_beam(el_num, new_az_num, theta, alpha, size, data);
    }
}

void PolarScan::merge_beam(int el_num, int az_num, double theta, double alpha, unsigned size, const double* dati)
{
    //LOG_CATEGORY("radar.io");
    // if (az_num >= vol_pol[el_num].size())
    //     vol_pol[el_num].resize(az_num + 1);
    beam_info[az_num].load_log.log(theta, alpha);

    unsigned overlap = min(beam_size, size);
    for (unsigned i = 0; i < overlap; ++i)
        if (get_db(az_num, i) < dati[i])
            set_db(az_num, i, dati[i]);

    beam_info[az_num].elevation = theta;
    //raggio.b_header.tipo_gran = INDEX_Z;  // FIXME: to be changed when we load different quantities
}

void VolumeStats::print(FILE* out)
{
    fprintf(out, "Nel    Zeros     Ones   Others      Sum\n");
    for (size_t iel =0; iel<count_zeros.size(); ++iel){
        fprintf(out, "%4zu %8u %8u %8u %8u\n",iel,count_zeros[iel],count_ones[iel],count_others[iel],sum_others[iel]);
    }
}

Volume::Volume()
    : acq_date(0), size_cell(0), declutter_rsp(false), NEL(0)
{
}

Volume::~Volume()
{
    for (vector<PolarScan*>::iterator i = scans.begin(); i != scans.end(); ++i)
        if (*i)
            delete *i;
}

Volume::LoadOptions::LoadOptions(const Site& site, bool medium, bool clean)
    : site(site), medium(medium), clean(clean), elev_array(site.get_elev_array(medium))
{
}

int Volume::LoadOptions::elevation_index(double elevation) const
{
    for (unsigned i=0; i < elev_array.size(); ++i)
        if (elevation >= (elev_array[i]-0.5) && elevation < (elev_array[i]+0.5))
            return i;
    return -1;
}


PolarScan& Volume::make_scan(const LoadOptions& opts, unsigned idx, unsigned beam_size)
{
    // Enlarge the scans vector if needed
    if (idx >= scans.size())
    {
        scans.resize(idx + 1, 0);
        NEL = idx + 1;
    }

    // Create the PolarScan if needed
    if (!scans[idx])
    {
        scans[idx] = new PolarScan(beam_size);
        scans[idx]->elevation = opts.elev_array[idx];
    }
    else if (beam_size != scans[idx]->beam_size)
    {
        LOG_CATEGORY("radar.io");
        LOG_ERROR("make_scan(idx=%u, beam_size=%u) called, but the scan already existed with beam_size=%u", idx, beam_size, scans[idx]->beam_size);
        throw runtime_error("beam_size mismatch");
    }

    // Return it
    return *scans[idx];
}

const unsigned Volume::max_beam_count() const
{
    unsigned res = 0;
    for (size_t i = 0; i < scans.size(); ++i)
        res = max(res, scans[i]->beam_count);
    return res;
}

const unsigned Volume::max_beam_size() const
{
    unsigned res = 0;
    for (size_t i = 0; i < scans.size(); ++i)
        res = max(res, scans[i]->beam_size);
    return res;
}

double Volume::elevation_min() const
{
    return scan(0).elevation;
}

double Volume::elevation_max() const
{
    return scan(NEL - 1).elevation;
}

void Volume::read_sp20(const char* nome_file, const LoadOptions& opts)
{
 LOG_CATEGORY("Volume");
    // dimensioni cella a seconda del tipo di acquisizione
    static const float size_cell_by_resolution[]={62.5,125.,250.,500.,1000.,2000.};
    //LOG_CATEGORY("radar.io");
    HD_DBP_SP20_RAW hd_char;

    filename = nome_file;

    // Replicato qui la read_dbp_SP20, per poi metterci mano e condividere codice con la lettura di ODIM

    FILE* sp20_in = fopen_checked(nome_file, "rb", "input sp20 file");

    // Read and decode file header
    if (fread(&hd_char, sizeof(hd_char), 1, sp20_in) != 1)
    {
      fclose(sp20_in);
      throw std::runtime_error("errore lettura header SP20");
    }
    HD_DBP_SP20_DECOD hd_file;
    decode_header_DBP_SP20(&hd_char, &hd_file);

    /*--------
      ATTENZIONE PRENDO LA DATA DAL NOME DEL FILE
      -------*/
    struct tm data_nome;
    acq_date = get_date_from_name(0, &data_nome, nome_file);
    size_cell = size_cell_by_resolution[(int)hd_file.cell_size];
    declutter_rsp = (bool)hd_file.filtro_clutter;

    BeamCleaner cleaner;
    cleaner.bin_wind_magic_number = opts.site.get_bin_wind_magic_number(acq_date);

    auto_ptr<Beams> b(new Beams);

    /* Ciclo su tutti i raggi del volume */
    while (true)
    {
      char beam_header[40];
      size_t res = fread(beam_header, sizeof(beam_header), 1, sp20_in);
      if (res != 1)
      {
        if (feof(sp20_in))
          break;
        else
        {
          string errmsg("Error reading ");
          errmsg += filename;
          errmsg += ": ";
          errmsg += strerror(errno);
          fclose(sp20_in);
          throw runtime_error(errmsg);
        }
      }

      BEAM_HD_SP20_INFO beam_info;
      decode_header_sp20_date_from_name(beam_header, &beam_info, (char*)nome_file);

      // Salta beam non validi
      if (!beam_info.valid_data) continue;

      // Calcola la nuova dimensione dei raggi
      float my_max_range = 123500;
      unsigned max_range;
      if (opts.clean)
          max_range = get_new_cell_num(beam_info.cell_num, my_max_range / size_cell);
      else
          max_range = beam_info.cell_num;

      // TODO: controllare il valore di ritorno delle fread
      if (beam_info.flag_quantities[0]) fread(b->data_z, 1, beam_info.cell_num, sp20_in);
      if (beam_info.flag_quantities[1]) fread(b->data_d, 1, beam_info.cell_num, sp20_in);
      if (beam_info.flag_quantities[2]) fread(b->data_v, 1, beam_info.cell_num, sp20_in);
      if (beam_info.flag_quantities[3]) fread(b->data_w, 1, beam_info.cell_num, sp20_in);

      vector<bool> cleaned(max_range, false);
      if (opts.clean)
          cleaner.clean_beams(*b, max_range, cleaned);

      int el_num = opts.elevation_index(beam_info.elevation);
      if (el_num < 0) continue;
      PolarScan& scan = make_scan(opts, el_num, max_range);

      // Convert to DB
      double* dbs = new double[max_range];
      for (unsigned i = 0; i < max_range; ++i)
          dbs[i] = BYTEtoDB(b->data_z[i]);
      //scan.elevation = beam_info.elevation;
#ifdef IMPRECISE_AZIMUT
      scan.fill_beam(el_num, beam_info.elevation, (int)(beam_info.azimuth / FATT_MOLT_AZ)*FATT_MOLT_AZ, max_range, dbs);
#else
      scan.fill_beam(el_num, beam_info.elevation, beam_info.azimuth, max_range, dbs);
#endif
      delete[] dbs;
    }

    fclose(sp20_in);

    LOG_DEBUG ("Nel volume ci sono %d scan",scans.size());
    for (int i=0; i<scans.size(); i++)  LOG_DEBUG (" Scan %2d - dimensione beam %5d",i, scans[i]->beam_size);
    // printf("NEL %d\n", (int)old_data_header.norm.maq.num_el);  // TODO: usare questo invece di NEL
    // for (int i = 0; i < old_data_header.norm.maq.num_el; ++i)
    //     printf("VALUE %d %d\n", i, old_data_header.norm.maq.value[i]); // Questi non so se ci servono

    // Initialize the rest of the volume after all beams are loaded
    resize_elev_fin();
}

double eldes_converter_azimut(double start, double stop)
{
    unsigned short azStart = (int)(start / 360.0 * 8192);
    unsigned short azStop = (int)(stop / 360.0 * 8192);

    //calcolo la media di azimuth e elevazione
    unsigned short azAvg = (azStart + azStop) / 2;

    //4096 e' la meta (cioe' 180 gradi) di un intero a 13 bit
    //13 bit a 1 (360 gradi) corrispondono al valore 8191
    //se si supera 4096 si sottrae 4096 sia per az che per ele
    //se sono a cavallo di 180 gradi·
    if( (azStart > (4096 + 2048) && azStop < 2048) ||
        (azStart < 2048 && azStop > (4096 + 2048)) )
        azAvg += 4095;

    azAvg &= 0x1FFF;

    double res = (double)azAvg * 360.0 / 8192.0;
#ifdef IMPRECISE_AZIMUT
    // Further truncate to 12 bits
    res = (int)(res / FATT_MOLT_AZ) * FATT_MOLT_AZ;
#endif
    return res;
}

unsigned char eldes_counter_to_db(unsigned short val)
{
    const int minScale = -20;
    float fVal = -31.5f + val * ((96.0f + 31.5f) / 65535);
    float ret = (fVal - minScale) / 0.3125f;
    if( ret < 0 )
        return 0;
    else if( ret > 255 )
        return 255;
    else
        return (unsigned char)(ret + 0.5f);
}

namespace {

unsigned int_to_unsigned(int val, const char* desc)
{
    if (val < 0)
    {
        string errmsg(desc);
        errmsg += " has a negative value";
        throw std::runtime_error(errmsg);
    }
    return (unsigned)val;
}

}

void Volume::read_odim(const char* nome_file, const LoadOptions& opts)
{
    LOG_CATEGORY("radar.io");

    namespace odim = OdimH5v21;
    using namespace Radar;
    using namespace std;

    filename = nome_file;

    auto_ptr<odim::OdimFactory> factory(new odim::OdimFactory());
    auto_ptr<odim::PolarVolume> volume(factory->openPolarVolume(nome_file));

    acq_date = volume->getDateTime();

    std::vector<double> elevationAngles = volume->getElevationAngles();

    // Check that the levels match what we want
    for (unsigned i = 0; i < elevationAngles.size(); ++i)
    {
        double v = elevationAngles[i];
        if (v <= opts.elev_array[i] - 0.5 || opts.elev_array[i] + 0.5 <= v)
        {
            LOG_ERROR("elevation %f does not match our expectation: we want %d but we got %f", i, opts.elev_array[i], v);
            throw runtime_error("elevation mismatch");
        }
    }

    double range_scale;

    // Iterate all scans
    unsigned scan_count = int_to_unsigned(volume->getScanCount(), "scan count");
    for (unsigned src_elev = 0; src_elev < scan_count; ++src_elev)
    {
        auto_ptr<odim::PolarScan> scan(volume->getScan(src_elev));
        double elevation = scan->getEAngle();

        // Read and and validate resolution information
        if (src_elev == 0)
            range_scale = scan->getRangeScale();
        else {
            double rs = scan->getRangeScale();
            if (rs != range_scale)
            {
                LOG_ERROR("scan %d (elevation %f) has rangeScale %f that is different from %f in the previous scans",
                        src_elev, elevation, rs, range_scale);
                throw runtime_error("rangeScale mismatch");
            }
        }

        // Pick the best quantity among the ones available
        auto_ptr<odim::PolarScanData> data;
        if (scan->hasQuantityData(odim::PRODUCT_QUANTITY_DBZH))
            data.reset(scan->getQuantityData(odim::PRODUCT_QUANTITY_DBZH));
        else if (scan->hasQuantityData(odim::PRODUCT_QUANTITY_TH))
        {
            LOG_WARN("no DBZH found for elevation angle %f: using TH", elevation);
            data.reset(scan->getQuantityData(odim::PRODUCT_QUANTITY_TH));
        }
        else
        {
            LOG_WARN("no DBZH or TH found for elevation angle %f", elevation);
            continue;
        }

        // Get and validate the azimuth angles for this scan
        std::vector<odim::AZAngles> azangles = scan->getAzimuthAngles();
        //int rpm_sign = scan->getDirection();

        unsigned nrays = int_to_unsigned(data->getNumRays(), "number of rays");
        if (azangles.size() != nrays)
        {
            LOG_ERROR("elevation %f has %zd azimuth angles and %d rays", elevation, azangles.size(), nrays);
            throw runtime_error("mismatch between number of azumuth angles and number of rays");
        }

        std::vector<double> elevation_angles = scan->getElevationAngles();

        unsigned beam_size = int_to_unsigned(data->getNumBins(), "beam size");
        if (beam_size > MAX_DIM)
        {
            LOG_ERROR("elevation %f has beams of %d elements, but the maximum we can store is %d", elevation, beam_size, MAX_DIM);
            throw runtime_error("beam size too big");
        }

        // Read all scan beam data

        // RayMatrix<float> matrix;
        // data->readTranslatedData(matrix);
        odim::RayMatrix<unsigned short> matrix;
        matrix.resize(nrays, beam_size);
        //data->readData(matrix);
        data->readData(const_cast<unsigned short*>(matrix.get()));

        double* beam = new double[beam_size];

        int el_num = opts.elevation_index(elevation);
        if (el_num < 0) continue;
        PolarScan& vol_pol_scan = make_scan(opts, el_num, beam_size);
        //vol_pol_scan.elevation = elevation;
        std::vector<bool> angles_seen(400, false);
        for (unsigned src_az = 0; src_az < nrays; ++src_az)
        {
            // FIXME: reproduce a bad truncation from the eldes sp20 converte
            //double azimut = azangles[src_az].averagedAngle(rpm_sign);
            double azimut = eldes_converter_azimut(azangles[src_az].start, azangles[src_az].stop);

            // Convert back to bytes, to fit into vol_pol as it is now
            for (unsigned i = 0; i < beam_size; ++i){
                // FIXME: QUESTO PEZZO DI CODICE E' STATO INSERITO PER EMULARE LA CONVERSIONE ELDES IN FORMATO SP20
                // DEVE ESSERE RIMOSSO A FINE LAVORO E RIATTIVATA QUESTA LINEA DI CODICE ORA COMMENTATA
                // beam[i] = DBtoBYTE(matrix.elem(src_az, i));
                beam[i] = BYTEtoDB(eldes_counter_to_db(matrix.elem(src_az, i)));
            }
            //vol_pol_scan.fill_beam(el_num, elevation, azimut, beam_size, beam);
            vol_pol_scan.fill_beam(el_num, elevation_angles[src_az], azimut, beam_size, beam);
        }

        delete[] beam;
    }

    size_cell = range_scale;

    resize_elev_fin();
}

void Volume::compute_stats(VolumeStats& stats) const
{
    stats.count_zeros.resize(scans.size());
    stats.count_ones.resize(scans.size());
    stats.count_others.resize(scans.size());
    stats.sum_others.resize(scans.size());

    for (int iel = 0; iel < scans.size(); ++iel)
    {
        stats.count_zeros[iel] = 0;
        stats.count_ones[iel] = 0;
        stats.count_others[iel] = 0;
        stats.sum_others[iel] = 0;

        for (unsigned iaz = 0; iaz < scan(iel).beam_count; ++iaz)
        {
            for (size_t i = 0; i < scan(iel).beam_size; ++i)
            {
                int val = scan(iel).get_raw(iaz, i);
                switch (val)
                {
                    case 0: stats.count_zeros[iel]++; break;
                    case 1: stats.count_ones[iel]++; break;
                    default:
                            stats.count_others[iel]++;
                            stats.sum_others[iel] += val;
                            break;
                }
            }
        }
    }
}

void Volume::resize_elev_fin()
{
    // FIXME: set to 0 to have the right size. We start from 512 (MAX_BIN)
    // to allocate enough memory for legacy code that iterates on MAX_BIN
    // to successfully read zeroes
    unsigned max_size = 512;
    for (unsigned iel = 0; iel < scans.size(); ++iel)
    {
        if (scan(iel).beam_size && scan(iel).beam_size > max_size)
            max_size = scan(iel).beam_size;
    }

    for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
    {
        elev_fin[i].resize(max_size, 0);
    }
}

void Volume::write_info_to_debug_file(H5::H5File out)
{
    using namespace H5;

    // Compute dimensions
    hsize_t dims[2] = { NUM_AZ_X_PPI, 0 };
    for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
        if (dims[1] < elev_fin[i].size())
            dims[1] = elev_fin[i].size();

    DataSpace file_data_space(2, dims);

    // Dataset data type
    IntType datatype( PredType::NATIVE_UCHAR );

    // Dataset fill value
    DSetCreatPropList props;
    unsigned char fill_value(0);
    props.setFillValue(datatype, &fill_value);

    // Create the dataset
    DataSet ds = out.createDataSet("/elev_fin", datatype, file_data_space, props);

    // Write elev_fin to it
    for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
    {
        hsize_t mdims[1] = { elev_fin[i].size() };
        DataSpace memory_data_space(1, mdims);

        hsize_t count[] = { 1, elev_fin[i].size() };
        hsize_t start[] = { i, 0 };
        file_data_space.selectHyperslab(H5S_SELECT_SET, count, start);

        ds.write(elev_fin[i].data(), datatype, memory_data_space, file_data_space);
    }
}

}
