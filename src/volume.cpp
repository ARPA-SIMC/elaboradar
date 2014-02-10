#include "volume.h"
#include "utils.h"
#include "logging.h"
#include <radarlib/radar.hpp>
#include <stdexcept>
#include <memory>
#include <cstring>

#ifdef __cplusplus
extern "C" {
#endif
// libreria radar
#include <func_SP20read.h>
#ifdef __cplusplus
}
#endif

using namespace std;

/// This needs to be a global variable, as it is expected by libsp20
int elev_array[NEL];

namespace cumbac {

Ray::Ray()
    : alfa_true(0), teta_true(0), teta(0), alfa(0)
{
}

void Ray::print_load_log(FILE* out) const
{
    if (load_log.empty())
    {
        fprintf(out, "no beams loaded\n");
        return;
    }

    for (vector<LoadLogEntry>::const_iterator i = load_log.begin(); i != load_log.end(); ++i)
    {
        if (i != load_log.begin())
            fprintf(out, ", ");
        fprintf(out, "ϑ%.2f α%.2f", i->theta, i->alpha);
    }
    fprintf(out, "\n");
}

PolarScan::PolarScan()
{
    //resize(NUM_AZ_X_PPI);
}

void VolumeStats::print(FILE* out)
{
    fprintf(out, "Nel   Zeros    Ones  Others     Sum\n"); 
    for (int iel =0; iel<NEL; ++iel){
        fprintf(out, "%4d%8d%8d%8d%8d\n",iel,count_zeros[iel],count_ones[iel],count_others[iel],sum_others[iel]);
    }
}

Volume::Volume()
    : acq_date(0), size_cell(0), declutter_rsp(false)
{
    memset(nbeam_elev, 0, sizeof(nbeam_elev));
}

void Volume::fill_beam(double theta, double alpha, unsigned size, const unsigned char* data)
{
    int teta = theta / FATT_MOLT_EL;

    int el_num = elevation_index_MDB(teta);
    if (el_num >= NEL) return;

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
        merge_beam(el_num, new_az_num, theta, alpha, size, data);
    }
    else if(az_num*0.9 - alpha > 0.)
    {
        int new_az_num = (az_num -1+400) %400;
        merge_beam(el_num, new_az_num, theta, alpha, size, data);
    }
}

void Volume::merge_beam(int el_num, int az_num, double theta, double alpha, unsigned size, const unsigned char* dati)
{
    //LOG_CATEGORY("radar.io");
    // if (az_num >= vol_pol[el_num].size())
    //     vol_pol[el_num].resize(az_num + 1);

    Ray& raggio = vol_pol[el_num][az_num];
    raggio.log(theta, alpha);

    if (raggio.ray.empty())
    {
        raggio.ray.reserve(size);

        for (unsigned i = 0; i < size; i++)
        {
            if (dati[i])
                raggio.ray.push_back(dati[i]);
            else
                raggio.ray.push_back(1);
        }
        nbeam_elev[el_num]++;
    }
    else
    {
        if (raggio.ray.size() < size)
        {
            raggio.ray.resize(size, 1);
            //LOG_WARN("volume[%d][%d]: old ray size: %zd, new ray size: %u", el_num, az_num, raggio.ray.size(), size);
            //throw runtime_error("attempted to merge two beams of different size");
        }

        for (unsigned i = 0; i < size; i++)
            if(raggio.ray[i] < dati[i])
                raggio.ray[i] = dati[i];
    }

    raggio.alfa =(short)(az_num*.9/FATT_MOLT_AZ);
    raggio.teta = elev_array[el_num];
    raggio.alfa_true = alpha / FATT_MOLT_AZ;
    raggio.teta_true = theta / FATT_MOLT_EL;
    //raggio.b_header.tipo_gran = INDEX_Z;  // FIXME: to be changed when we load different quantities
}

void Volume::read_sp20(const char* nome_file)
{
    // dimensioni cella a seconda del tipo di acquisizione
    static const float size_cell_by_resolution[]={62.5,125.,250.,500.,1000.,2000.};
    //LOG_CATEGORY("radar.io");
    T_MDB_data_header   old_data_header;

    filename = nome_file;

    // Replicato qui la read_dbp_SP20, per poi metterci mano e condividere codice con la lettura di ODIM

    FILE* sp20_in = fopen_checked(nome_file, "rb", "input sp20 file");

    if (ReadHeaderSP20toMDB(&old_data_header, sp20_in) == NO_OK)
    {
        fclose(sp20_in);
        throw std::runtime_error("errore lettura header SP20");
    }

    /*--------
      ATTENZIONE PRENDO LA DATA DAL NOME DEL FILE
      -------*/
    struct tm data_nome;
    old_data_header.norm.maq.acq_date=get_date_from_name(&old_data_header,&data_nome,nome_file);

    T_MDB_ap_beam_header  old_beam_header;
    unsigned char dati[MAX_DIM];
    const int tipo_dati = INDEX_Z;
    while(1)
    {
        if(read_ray_SP20(&old_beam_header,dati,sp20_in,tipo_dati)==NO_OK) break;
        fill_beam(old_beam_header.teta * FATT_MOLT_EL, old_beam_header.alfa * FATT_MOLT_AZ, old_beam_header.max_bin, dati);
    }

    fclose(sp20_in);

    // printf("NEL %d\n", (int)old_data_header.norm.maq.num_el);  // TODO: usare questo invece di NEL
    // for (int i = 0; i < old_data_header.norm.maq.num_el; ++i)
    //     printf("VALUE %d %d\n", i, old_data_header.norm.maq.value[i]); // Questi non so se ci servono

    acq_date = old_data_header.norm.maq.acq_date;
    size_cell = size_cell_by_resolution[old_data_header.norm.maq.resolution];
    declutter_rsp = (bool)old_data_header.norm.maq.declutter_rsp;

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
    // Further truncate to 12 bits
    res = (int)(res / FATT_MOLT_AZ) * FATT_MOLT_AZ;
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

void Volume::read_odim(const char* nome_file)
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

    // Make sure that we can store all the levels in the scan
    if (elevationAngles.size() > NEL)
    {
        LOG_INFO("%zd elevation angles found, but we can only store %d", elevationAngles.size(), NEL);
        throw runtime_error("number of elevation angles too big");
    }

    // Check that the levels match what we want
    for (unsigned i = 0; i < elevationAngles.size(); ++i)
    {
        double v = elevationAngles[i] * 4096.0 / 360.0;
        if (v <= elev_array[i] - 5 || elev_array[i] + 5 <= v)
        {
            LOG_ERROR("elevation %d does not match our expectation: we want %d but we got %f", i, elev_array[i], v);
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

        unsigned char* beam = new unsigned char[beam_size];

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
                beam[i] = eldes_counter_to_db(matrix.elem(src_az, i));
            }
            fill_beam(elevation, azimut, beam_size, beam);
        }

        delete[] beam;
    }

    size_cell = range_scale;

    resize_elev_fin();
}

void Volume::compute_stats(VolumeStats& stats) const
{
    for (int iel = 0; iel < NEL; ++iel)
    {
        stats.count_zeros[iel] = 0;
        stats.count_ones[iel] = 0;
        stats.count_others[iel] = 0;
        stats.sum_others[iel] = 0;

        for (unsigned ibeam = 0; ibeam < nbeam_elev[iel]; ++ibeam)
        {
            for (size_t i = 0; i < vol_pol[iel][ibeam].ray.size(); ++i)
            {
                int val = vol_pol[iel][ibeam].ray[i];
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
    for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
    {
        // FIXME: set to 0 to have the right size. We start from 512 (MAX_BIN)
        // to allocate enough memory for legacy code that iterates on MAX_BIN
        // to successfully read zeroes
        unsigned max_size = 512;
        for (unsigned iel = 0; iel < NEL; ++iel)
        {
            if (vol_pol[iel][i].ray.size() && vol_pol[iel][i].ray.size() > max_size)
                max_size = vol_pol[iel][i].ray.size();
        }
        elev_fin[i].resize(max_size, 0);
    }
}

}
