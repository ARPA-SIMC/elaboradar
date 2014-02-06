#include "volume.h"
#include "utils.h"
#include "logging.h"
#include <radarlib/radar.hpp>
#include <stdexcept>
#include <memory>
#include <cstring>

/// This needs to be a global variable, as it is expected by libsp20
int elev_array[NEL];

namespace cumbac {

Volume::Volume()
    : acq_date(0), size_cell(0), declutter_rsp(false)
{
    memset(vol_pol, 0, sizeof(vol_pol));
    memset(nbeam_elev, 0, sizeof(nbeam_elev));
}

void Volume::fill_beam(double theta, double alpha, unsigned size, const unsigned char* data)
{
    int teta = theta / FATT_MOLT_EL;

    int el_num = elevation_index_MDB(teta);
    if (el_num > NEL) return;

    int alfa = alpha / FATT_MOLT_AZ;
    if (alfa >= 4096) return;

    int az_num = azimut_index_MDB(alfa);

    if (az_num == 0)
    {
        printf("fbeam ϑ%f→%d α%f→%d %u", theta, el_num, alpha, az_num, size);
        for (unsigned i = 0; i < 20; ++i)
            printf(" %d", (int)data[i]);
        printf("\n");
    }

    merge_beam(&vol_pol[el_num][az_num], theta, alpha, az_num, el_num, size, data);
    if(az_num*0.9 - alpha < 0.)
    {
        int new_az_num = (az_num + 1) % 400;
        merge_beam(&vol_pol[el_num][new_az_num], theta, alpha, new_az_num, el_num, size, data);
    }
    else if(az_num*0.9 - alpha > 0.)
    {
        int new_az_num = (az_num -1+400) %400;
        merge_beam(&vol_pol[el_num][new_az_num], theta, alpha, new_az_num, el_num, size, data);
    }
}

void Volume::merge_beam(VOL_POL* raggio, double theta, double alpha, int az_num, int el_num, unsigned size, const unsigned char* dati)
{
    if (raggio->flag == 0)
    {
        for (unsigned i = 0; i < size; i++)
        {
            if(dati[i])
                raggio->ray[i] = dati[i];
            else
                raggio->ray[i] = 1;
        }
        nbeam_elev[el_num]++;
    }
    else
        for (unsigned i = 0; i < size; i++)
            if(raggio->ray[i]<dati[i])
                raggio->ray[i]=dati[i];

    raggio->flag=1;
    raggio->b_header.alfa =(short)(az_num*.9/FATT_MOLT_AZ);
    raggio->b_header.teta = elev_array[el_num];
    raggio->alfa_true = alpha / FATT_MOLT_AZ;
    raggio->teta_true = theta / FATT_MOLT_EL;
    raggio->b_header.tipo_gran = INDEX_Z;  // FIXME: to be changed when we load different quantities
    raggio->b_header.max_bin = size;
}

void Volume::read_sp20(const char* nome_file)
{
    // dimensioni cella a seconda del tipo di acquisizione
    static const float size_cell_by_resolution[]={62.5,125.,250.,500.,1000.,2000.};

    LOG_CATEGORY("radar.io");

    T_MDB_data_header   old_data_header;

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
}

void Volume::read_odim(const char* nome_file)
{
    LOG_CATEGORY("radar.io");

    using namespace OdimH5v21;
    using namespace Radar;
    using namespace std;

    auto_ptr<OdimFactory> factory(new OdimFactory());
    auto_ptr<PolarVolume> volume(factory->openPolarVolume(nome_file));

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
    int scan_count = volume->getScanCount();
    for (unsigned src_elev = 0; src_elev < scan_count; ++src_elev)
    {
        auto_ptr<PolarScan> scan(volume->getScan(src_elev));
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
        auto_ptr<PolarScanData> data;
        if (scan->hasQuantityData(PRODUCT_QUANTITY_DBZH))
            data.reset(scan->getQuantityData(PRODUCT_QUANTITY_DBZH));
        else if (scan->hasQuantityData(PRODUCT_QUANTITY_TH))
        {
            LOG_WARN("no DBZH found for elevation angle %f: using TH", elevation);
            data.reset(scan->getQuantityData(PRODUCT_QUANTITY_TH));
        }
        else
        {
            LOG_WARN("no DBZH or TH found for elevation angle %f", elevation);
            continue;
        }

        // Get and validate the azimuth angles for this scan
        std::vector<AZAngles> azangles = scan->getAzimuthAngles();
        int rpm_sign = scan->getDirection();

        int nrays = data->getNumRays();
        if (azangles.size() != nrays)
        {
            LOG_ERROR("elevation %f has %zd azimuth angles and %d rays", elevation, azangles.size(), nrays);
            throw runtime_error("mismatch between number of azumuth angles and number of rays");
        }

        int beam_size = data->getNumBins();
        if (beam_size > MAX_DIM)
        {
            LOG_ERROR("elevation %f has beams of %d elements, but the maximum we can store is %d", elevation, beam_size, MAX_DIM);
            throw runtime_error("beam size too big");
        }

        // Read all scan beam data
        RayMatrix<float> matrix;
        data->readTranslatedData(matrix);

        //printf("Offset %f, gain %f\n", data->getOffset(), data->getGain());
        //double offset = data->getOffset();
        //double gain = data->getGain();

        unsigned char* beam = new unsigned char[beam_size];

        std::vector<bool> angles_seen(400, false);
        for (int src_az = 0; src_az < nrays; ++src_az)
        {
            double azimut = azangles[src_az].averagedAngle(rpm_sign);
            // printf("fbeam ϑ%5.1f α1%6.1f α2%6.1f α%6.1f sign %2d\n", elevation, azangles[src_az].start,  azangles[src_az].stop, azimut, rpm_sign);
            // Convert back to bytes, to fit into vol_pol as it is now
            for (unsigned i = 0; i < beam_size; ++i){
                // QUESTO PEZZO DI CODICE E' STATO INSERITO PER EMULARE LA CONVERSIONE ELDES IN FORMATO SP20
                // DEVE ESSERE RIMOSSO A FINE LAVORO E RIATTIVATA QUESTA LINEA DI CODICE ORA COMMENTATA
                //   beam[i] = DBtoBYTE(matrix.elem(src_az, i));
                typedef unsigned char byte;
                float ret = (matrix.elem(src_az, i)- (-20.0) ) / 0.3125f;
                if( ret < 0 )
                    beam[i] = 0;
                else if( ret > 255 )
                    beam[i] = 255;
                else
                    beam[i] = byte(ret + 0.5f);
            }
            fill_beam(elevation, azimut, beam_size, beam);
        }

        delete[] beam;
    }

    size_cell = range_scale;
}

}
