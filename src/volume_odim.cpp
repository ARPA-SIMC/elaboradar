#include "volume/odim.h"
#include "logging.h"
#include <radarlib/radar.hpp>
#include <memory>

using namespace std;

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

}


namespace cumbac {
namespace volume {

void ODIMLoader::make_scan(unsigned idx, unsigned beam_size)
{
    Loader::make_scan(idx, beam_size);
    if (vol_z) vol_z->make_scan(idx, beam_size, elev_array[idx]);
}

void ODIMLoader::load(const std::string& pathname)
{
    LOG_CATEGORY("radar.io");

    namespace odim = OdimH5v21;
    using namespace Radar;
    using namespace std;

    if (load_info) load_info->filename = pathname;

    auto_ptr<odim::OdimFactory> factory(new odim::OdimFactory());
    auto_ptr<odim::PolarVolume> volume(factory->openPolarVolume(pathname));

    if (load_info) load_info->acq_date = volume->getDateTime();

    std::vector<double> elevationAngles = volume->getElevationAngles();

    // Check that the levels match what we want
    for (unsigned i = 0; i < elevationAngles.size(); ++i)
    {
        double v = elevationAngles[i];
        if (v <= elev_array[i] - 0.5 || elev_array[i] + 0.5 <= v)
        {
            LOG_ERROR("elevation %d does not match our expectation: we want %f but we got %f", i, elev_array[i], v);
            throw runtime_error("elevation mismatch");
        }
    }

    double range_scale = 0;

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


/* 
 *  per permettere al cleaner di funzionare per dati ODIm bisogna fare i seguenti passi
 *  
 *  1) leggere Z come valore fisico  data->readTranslatedData(matrix);
 *
 *  2) verificare che siano presenti le grandezze PRODUCT_QUANTITY_VRAD e PRODUCT_QUANTITY_WRAD
 *
 *  3) leggere VRAD (valore fisico)
 *
 *  4) leggere WRAD valore fisico)
 *
 *  5) calcolare bin_wind_magic_number  -> = VRAD.undetect*VRAD.gain+VRAD.offset
 *
 *  6) calcolare z_missing  -> noData
 *
 *  7) calcolare v_missing -> noData
 *
 *  8) soglia W = 0.
 *
 *  9) riempire la struttura Beams
 *
 * 10) eseguire il cleaner sulla struttura
 *
 * 11) copiare i dati ripuliti in volume
 *
 * 12) nuovo raggio si ricomincia da 1 fino a fine PolarScan
 *
 * */


        // Read all scan beam data

        // RayMatrix<float> matrix;
        odim::RayMatrix<unsigned short> matrix;
        matrix.resize(nrays, beam_size);
        // data->readTranslatedData(matrix);
        data->readData(const_cast<unsigned short*>(matrix.get()));

        double* beam = new double[beam_size];

        int el_num = elevation_index(elevation);
        if (el_num < 0) continue;

        make_scan(el_num, beam_size);
        PolarScan<double>& vol_pol_scan = vol_z->scan(el_num);

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
            fill_beam(vol_pol_scan, el_num, elevation_angles[src_az], azimut, beam_size, beam);
        }

        delete[] beam;
    }

    if (load_info) load_info->size_cell = range_scale;
}


}
}
