#include "volume/sp20.h"
#include "logging.h"
#include "utils.h"
#include "volume_cleaner.h"
#include "site.h"
#include <memory>
#include <cerrno>
#include <cstring>

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
namespace volume {

namespace sp20 {

struct Beam
{
    BEAM_HD_SP20_INFO beam_info;
    Beams<unsigned char> beams;
    unsigned beam_size = 0;

    // Read the next beam in the file. Returns false on end of file.
    bool read(File& in)
    {
        char beam_header[40];
        if (!in.fread(beam_header, sizeof(beam_header)))
            return false;

        decode_header_sp20_date_from_name(beam_header, &beam_info, (char*)in.name());

        // Salta beam non validi
        if (!beam_info.valid_data)
            return true;

        beam_size = beam_info.cell_num;

        if (has_z()) in.fread(beams.data_z, beam_info.cell_num);
        if (has_d()) in.fread(beams.data_d, beam_info.cell_num);
        if (has_v()) in.fread(beams.data_v, beam_info.cell_num);
        if (has_w()) in.fread(beams.data_w, beam_info.cell_num);

        return true;
    }

    bool has_z() const { return beam_info.flag_quantities[0]; }
    bool has_d() const { return beam_info.flag_quantities[1]; }
    bool has_v() const { return beam_info.flag_quantities[2]; }
    bool has_w() const { return beam_info.flag_quantities[3]; }
};

}

void SP20Loader::make_scan(unsigned idx, unsigned beam_size, double cell_size)
{
    Loader::make_scan(idx, beam_size);
    if (vol_z) vol_z->make_scan(idx, beam_size, elev_array[idx], cell_size);
    if (vol_d) vol_d->make_scan(idx, beam_size, elev_array[idx], cell_size);
    if (vol_v) vol_v->make_scan(idx, beam_size, elev_array[idx], cell_size);
    if (vol_w) vol_w->make_scan(idx, beam_size, elev_array[idx], cell_size);
}

namespace {

struct Beams : public std::vector<sp20::Beam*>
{
    // Elevation number
    unsigned el_num;

    Beams(unsigned el_num) : el_num(el_num) {}
    Beams(const Beams&) = delete;
    Beams(const Beams&&) = delete;
    ~Beams()
    {
        for (auto i: *this)
            delete i;
    }
    Beams& operator=(const Beams&) = delete;

    // Compute the minimum beam size
    unsigned min_beam_size() const
    {
        if (empty()) return 0;
        unsigned res = (*this)[0]->beam_size;
        for (auto i: *this)
            res = min(res, i->beam_size);
        return res;
    }

    // Set the beam size on each beam
    void set_beam_size(unsigned val)
    {
        for (auto i: *this)
            i->beam_size = val;
    }
};

struct Elevations : public std::vector<unique_ptr<Beams>>
{
};

}

void SP20Loader::load(const std::string& pathname)
{
    LOG_CATEGORY("Volume");
    // dimensioni cella a seconda del tipo di acquisizione
    static const float size_cell_by_resolution[]={62.5,125.,250.,500.,1000.,2000.};
    //LOG_CATEGORY("radar.io");
    HD_DBP_SP20_RAW hd_char;

    if (load_info) load_info->filename = pathname;

    // Replicato qui la read_dbp_SP20, per poi metterci mano e condividere codice con la lettura di ODIM

    File sp20_in(logging_category);
    sp20_in.open(pathname, "rb", "input sp20 file");
    if (!sp20_in)
        throw std::runtime_error("failed to open sp20 input file");

    // Read and decode file header
    sp20_in.fread(&hd_char, sizeof(hd_char));

    HD_DBP_SP20_DECOD hd_file;
    decode_header_DBP_SP20(&hd_char, &hd_file);

    /*--------
      ATTENZIONE PRENDO LA DATA DAL NOME DEL FILE
      -------*/
    struct tm data_nome;
    time_t acq_date = get_date_from_name(0, &data_nome, pathname.c_str());
    if (load_info) load_info->acq_date = acq_date;
    if (load_info) load_info->declutter_rsp = (bool)hd_file.filtro_clutter;
    double size_cell = size_cell_by_resolution[(int)hd_file.cell_size];

    BeamCleaner<unsigned char> cleaner(site.get_bin_wind_magic_number(acq_date), 0,0);
//    cleaner.bin_wind_magic_number = site.get_bin_wind_magic_number(acq_date);

    // Read all beams from the file
    Elevations elevations;

    while (true)
    {
        unique_ptr<sp20::Beam> beam(new sp20::Beam);

        if (!beam->read(sp20_in))
            break;

        int el_num = elevation_index(beam->beam_info.elevation);
        if (el_num < 0) continue;

        // Calcola la nuova dimensione dei raggi
        if (max_bin)
            beam->beam_size = min(beam->beam_size, max_bin);

        vector<bool> cleaned(beam->beam_size, false);
        if (clean)
            cleaner.clean_beams(beam->beams, beam->beam_size, cleaned);

        while ((unsigned)el_num >= elevations.size())
        {
            elevations.push_back(unique_ptr<Beams>(new Beams(elevations.size())));
        }

        elevations[el_num]->push_back(beam.release());
    }

    for (auto& beams: elevations)
    {
        // Set the beam size on each elevation to the minimum beam size
        unsigned beam_size = beams->min_beam_size();
        beams->set_beam_size(beam_size);

        // TODO: check if an elevation has been skipped
        make_scan(beams->el_num, beam_size, size_cell);

        for (auto& beam: *beams)
            beam_to_volumes(*beam, beams->el_num);
    }

    LOG_DEBUG ("Nel volume ci sono %zd scan", vol_z->size());
    for (size_t i = 0; i < vol_z->size(); ++i)
        LOG_DEBUG (" Scan %2zd - dimensione beam %5d", i, vol_z->scan(i).beam_size);
    // printf("NEL %d\n", (int)old_data_header.norm.maq.num_el);  // TODO: usare questo invece di NEL
    // for (int i = 0; i < old_data_header.norm.maq.num_el; ++i)
    //     printf("VALUE %d %d\n", i, old_data_header.norm.maq.value[i]); // Questi non so se ci servono
}

void SP20Loader::beam_to_volumes(const sp20::Beam& beam, unsigned el_num)
{
    const unsigned max_range = beam.beam_size;

    if (vol_z && beam.has_z()) // Riflettività Z
    {
        // Convert to DB
        double* dbs = new double[max_range];
        for (unsigned i = 0; i < max_range; ++i)
            dbs[i] = BYTEtoDB(beam.beams.data_z[i]);
        //f_ray[p]=data->beam[p]*RANGE_Z/255. + min_zeta;

        PolarScan<double>& scan = vol_z->scan(el_num);
#ifdef IMPRECISE_AZIMUT
        fill_beam(scan, el_num, beam.beam_info.elevation, (int)(beam.beam_info.azimuth / FATT_MOLT_AZ)*FATT_MOLT_AZ, max_range, dbs);
#else
        fill_beam(scan, el_num, beam.beam_info.elevation, beam.beam_info.azimuth, max_range, dbs);
#endif
        delete[] dbs;
    }

    if (vol_d && beam.has_d()) // Riflettività differenziale ZDR
    {
        // range variabilita ZDR
        const double range_zdr = 16.;
        const double min_zdr = -6.;

        // Convert to DB 
        double* dbs = new double[max_range];
        for (unsigned i = 0; i < max_range; ++i)
            dbs[i] = beam.beams.data_d[i] * range_zdr / 255. + min_zdr;

        PolarScan<double>& scan = vol_d->scan(el_num);
#ifdef IMPRECISE_AZIMUT
        fill_beam(scan, el_num, beam.beam_info.elevation, (int)(beam.beam_info.azimuth / FATT_MOLT_AZ)*FATT_MOLT_AZ, max_range, dbs);
#else
        fill_beam(scan, el_num, beam.beam_info.elevation, beam.beam_info.azimuth, max_range, dbs);
#endif
    }

    if (vol_v && beam.has_v()) // Velocità V
    {
        // Convert to m/s 
        double* ms = new double[max_range];
        if (beam.beam_info.PRF == 'S')
        {
            // range variabilita V - velocità radiale
            const double range_v = 33.;
            for (unsigned i = 0; i < max_range; ++i)
                if (beam.beams.data_v[i] == -128)
                    ms[i] = -range_v / 2;
                else
                    ms[i] = beam.beams.data_v[i] * range_v / 254.;
        } else {
            // range variabilita V - velocità radiale
            const double range_v = 99.;
            for (unsigned i = 0; i < max_range; ++i)
                if (beam.beams.data_v[i] == -128)
                    ms[i] = -range_v / 2;
                else
                    ms[i] = beam.beams.data_v[i] * range_v / 254.;
        }
        // if (data->beam_w[p] == -128) data->beam_w[p] = -127;
        // if ( beam_info->PRF == 'S')
        //   f_ray[p] = data->beam_w[p] * RANGE_V / 127.*.5;
        // else
        //   f_ray[p] = data->beam_w[p] * RANGE_V2 / 127.*.5;
        PolarScan<double>& scan = vol_v->scan(el_num);
#ifdef IMPRECISE_AZIMUT
        fill_beam(scan, el_num, beam.beam_info.elevation, (int)(beam.beam_info.azimuth / FATT_MOLT_AZ)*FATT_MOLT_AZ, max_range, ms);
#else
        fill_beam(scan, el_num, beam.beam_info.elevation, beam.beam_info.azimuth, max_range, ms);
#endif
    }

    if (vol_w && beam.has_w()) // Spread - Sigma V
    {
        // range variabilita Sigma V - Spread velocità
        const double range_sig_v = 10.;
        // Convert to m/s 
        double* ms = new double[max_range];
        for (unsigned i = 0; i < max_range; ++i)
            ms[i] = beam.beams.data_w[i] * range_sig_v / 255.0;

        PolarScan<double>& scan = vol_w->scan(el_num);
#ifdef IMPRECISE_AZIMUT
        fill_beam(scan, el_num, beam.beam_info.elevation, (int)(beam.beam_info.azimuth / FATT_MOLT_AZ)*FATT_MOLT_AZ, max_range, ms);
#else
        fill_beam(scan, el_num, beam.beam_info.elevation, beam.beam_info.azimuth, max_range, ms);
#endif
    }
}

}
}
