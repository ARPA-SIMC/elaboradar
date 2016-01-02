#include "sp20.h"
#include "logging.h"
#include "utils.h"
#include <radarlib/radar.hpp>
#include <memory>
#include <cerrno>
#include <cstring>
#include "sp20/func_SP20read.h"

#ifdef NEL
#undef NEL
#endif

using namespace std;

/// This is not used anymore, but it is here to satisfy libSP20 linking needs
int elev_array[15];

namespace elaboradar {
namespace volume {

namespace sp20 {

struct Beam
{
    BEAM_HD_SP20_INFO beam_info;
    unsigned char data_z[1024];
    unsigned char data_d[1024];
    unsigned char data_v[1024];
    unsigned char data_w[1024];
    unsigned beam_size = 0;

    // Read the next beam in the file. Returns false on end of file.
    bool read(File& in)
    {
        unsigned char beam_header[40];
        if (!in.fread(beam_header, sizeof(beam_header)))
            return false;

        decode_header_sp20_date_from_name(beam_header, &beam_info, (char*)in.name());

        // Salta beam non validi
        if (!beam_info.valid_data){
            unsigned to_be_skipped=0;
            if (has_z()) to_be_skipped=to_be_skipped+beam_info.cell_num;
            if (has_d()) to_be_skipped=to_be_skipped+beam_info.cell_num;
            if (has_v()) to_be_skipped=to_be_skipped+beam_info.cell_num;
            if (has_w()) to_be_skipped=to_be_skipped+beam_info.cell_num;
            in.fseek(to_be_skipped, SEEK_CUR);
            return true;
        }

        beam_size = beam_info.cell_num;

        if (has_z()) in.fread(data_z, beam_info.cell_num);
        if (has_d()) in.fread(data_d, beam_info.cell_num);
        if (has_v()) in.fread(data_v, beam_info.cell_num);
        if (has_w()) in.fread(data_w, beam_info.cell_num);

        return true;
    }

    bool has_z() const { return beam_info.flag_quantities[0]; }
    bool has_d() const { return beam_info.flag_quantities[1]; }
    bool has_v() const { return beam_info.flag_quantities[2]; }
    bool has_w() const { return beam_info.flag_quantities[3]; }
};

}

namespace {

struct Beams : public std::vector<sp20::Beam*>
{
    // Elevation
    double elevation;
    // Beam size for the polar scan
    unsigned beam_size = 0;

    Beams(double elevation) : elevation(elevation) {}
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
};

struct Elevations : public std::vector<unique_ptr<Beams>>
{
    Elevations(float* elevations, unsigned el_count)
    {
        // Create a sorted vector of elevations
        vector<double> els;
        for (unsigned i = 0; i < el_count; ++i)
            els.push_back(elevations[i]);
        std::sort(els.begin(), els.end());

        // Use it to initialize our Beams structures
        for (double el: els)
            emplace_back(new Beams(el));
    }

    Beams* get_closest(double elevation)
    {
        for (unsigned i=0; i < size(); ++i)
            if (elevation >= (at(i)->elevation-0.5) && elevation < (at(i)->elevation+0.5))
                return at(i).get();
        return nullptr;
    }
};

}

void SP20Loader::load(const std::string& pathname)
{
    namespace odim = OdimH5v21;
    LOG_CATEGORY("Volume");
    // dimensioni cella a seconda del tipo di acquisizione
    static const float size_cell_by_resolution[]={62.5,125.,250.,500.,1000.,2000.};
    //LOG_CATEGORY("radar.io");
    HD_DBP_SP20_RAW hd_char;

    shared_ptr<LoadInfo> load_info = make_shared<LoadInfo>();
    load_info->filename = pathname;

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
    load_info->acq_date = acq_date;
    load_info->declutter_rsp = (bool)hd_file.filtro_clutter;
    double size_cell = size_cell_by_resolution[(int)hd_file.cell_size];
    bool has_dual_prf = hd_file.Dual_PRF;

    // Read all beams from the file
    Elevations elevations(hd_file.ele, hd_file.num_ele);
    // TODO for (int i = 0; i < hd_file.num_ele; ++i)
    // TODO     fprintf(stderr, "%d: %f\n", i, (double)hd_file.ele[i]);



    while (true)
    {
        unique_ptr<sp20::Beam> beam(new sp20::Beam);

        if (!beam->read(sp20_in))
            break;

        Beams* beams = elevations.get_closest(beam->beam_info.elevation);
        if (!beams) continue;

        beams->push_back(beam.release());
    }

    if (elevations.empty())
        throw std::runtime_error("sp20 file has no beams");

    // Set the beam size on each elevation to the minimum beam size
    for (auto& beams: elevations)
    {
        unsigned beam_size = beams->min_beam_size();
        beams->beam_size = beam_size;
    }
//for (unsigned i=0; i<elevations.size(); ++i)  std::cout<<elevations[i]->size()<<std::endl;

    if (elevations.back()->beam_size == 0)
        throw std::runtime_error("last elevation beam size is 0");

    // For elevations that have beam_size = 0, use the beam size from the
    // elevation above
    for (int i = elevations.size() - 2; i >= 0; --i)
    {
        if (elevations[i]->beam_size == 0)
            elevations[i]->beam_size = elevations[i + 1]->beam_size;
    }

    // Create the polarscans and fill them with data
    for (unsigned idx = 0; idx < elevations.size(); ++idx)
    {
        const Beams& beams = *elevations[idx];

        // Create structures for a new PolarScan
        if (vol_z) vol_z->append_scan(beams.size(), beams.beam_size, beams.elevation, size_cell);
        if (vol_d) vol_d->append_scan(beams.size(), beams.beam_size, beams.elevation, size_cell);
        if (vol_v) vol_v->append_scan(beams.size(), beams.beam_size, beams.elevation, size_cell);
        if (vol_w) vol_w->append_scan(beams.size(), beams.beam_size, beams.elevation, size_cell);

        // Fill the new scan with data
        for (unsigned i = 0; i < beams.size(); ++i)
            beam_to_volumes(*beams.at(i), i, beams.beam_size, idx);
    }

    if (vol_z)
    {
        vol_z->load_info = load_info;
        if (load_info->declutter_rsp)
            vol_z->quantity = odim::PRODUCT_QUANTITY_DBZH;
        else
            vol_z->quantity = odim::PRODUCT_QUANTITY_TH;
        for (auto& scan : *vol_z)
        {
            scan.nodata = DBtoBYTE(0);
            scan.undetect = DBtoBYTE(1);
            scan.gain = 80.0 / 255.0;
            scan.offset = -20;
        }
    }
    if (vol_d)
    {
        vol_d->load_info = load_info;
        vol_v->quantity = odim::PRODUCT_QUANTITY_ZDR;
        for (auto& scan : *vol_z)
        {
            scan.nodata = -7;
            scan.undetect = -6;
            scan.gain = 16.0 / 255.0;
            scan.offset = -6;
        }
    }
    if (vol_v)
    {
        vol_v->load_info = load_info;
        vol_v->quantity = odim::PRODUCT_QUANTITY_VRAD;
        for (auto& scan : *vol_v)
        {
            if (has_dual_prf)
            {
                scan.nodata = -49.5;
                scan.undetect = -49.5;
                scan.gain = 99.0 / 254;       // TODO: check why it's 254 and not 255
                scan.offset = -49.5;
            } else {
                scan.nodata = -16.5;
                scan.undetect = -16.5;
                scan.gain = 33.0 / 254;       // TODO: check why it's 254 and not 255
                scan.offset = -16.5;
            }
        }
    }
    if (vol_w)
    {
        vol_w->load_info = load_info;
        vol_w->quantity = odim::PRODUCT_QUANTITY_WRAD;
        for (auto& scan : *vol_w)
        {
            scan.nodata = -1;
            scan.undetect = 0;
            scan.gain = 10.0 / 255.0;
            scan.offset = 0;
        }
    }

    LOG_DEBUG ("Nel volume ci sono %zd scan", vol_z->size());
    for (size_t i = 0; i < vol_z->size(); ++i)
        LOG_DEBUG (" Scan %2zd - dimensione beam %5d -  numero raggi %3d", i, vol_z->at(i).beam_size, vol_z->at(i).beam_count);
    // printf("NEL %d\n", (int)old_data_header.norm.maq.num_el);  // TODO: usare questo invece di NEL
    // for (int i = 0; i < old_data_header.norm.maq.num_el; ++i)
    //     printf("VALUE %d %d\n", i, old_data_header.norm.maq.value[i]); // Questi non so se ci servono
}

void SP20Loader::beam_to_volumes(const sp20::Beam& beam, unsigned az_idx, unsigned beam_size, unsigned el_num)
{
    const unsigned max_range = min(beam_size, beam.beam_size);

    if (vol_z && beam.has_z()) // Riflettività Z
    {
        // Convert to DB
        Eigen::VectorXd dbs(max_range);
        for (unsigned i = 0; i < max_range; ++i)
            dbs(i) = BYTEtoDB(beam.data_z[i]);
        PolarScan<double>& scan = vol_z->at(el_num);
        scan.row(az_idx) = dbs;
        scan.elevations_real(az_idx) = beam.beam_info.elevation;
        scan.azimuths_real(az_idx) = beam.beam_info.azimuth;
    }

    if (vol_d && beam.has_d()) // Riflettività differenziale ZDR
    {
        // range variabilita ZDR
        const double range_zdr = 16.;
        const double min_zdr = -6.;

        // Convert to DB 
        Eigen::VectorXd dbs(max_range);
        for (unsigned i = 0; i < max_range; ++i)
            dbs(i) = beam.data_d[i] * range_zdr / 255. + min_zdr;

        PolarScan<double>& scan = vol_d->at(el_num);
        scan.row(az_idx) = dbs;
        scan.elevations_real(az_idx) = beam.beam_info.elevation;
        scan.azimuths_real(az_idx) = beam.beam_info.azimuth;
    }

    if (vol_v && beam.has_v()) // Velocità V
    {
        // Convert to m/s 
        Eigen::VectorXd ms(max_range);
        if (beam.beam_info.PRF == 'S')
        {
            // range variabilita V - velocità radiale
            const double range_v = 33.;
            for (unsigned i = 0; i < max_range; ++i)
                if (beam.data_v[i] == -128)
                    ms(i) = -range_v / 2;
                else
                    ms(i) = beam.data_v[i] * range_v / 254.;
        } else {
            // range variabilita V - velocità radiale
            const double range_v = 99.;
            for (unsigned i = 0; i < max_range; ++i)
                if (beam.data_v[i] == -128)
                    ms(i) = -range_v / 2;
                else
                    ms(i) = beam.data_v[i] * range_v / 254.;
        }
        // if (data->beam_w[p] == -128) data->beam_w[p] = -127;
        // if ( beam_info->PRF == 'S')
        //   f_ray[p] = data->beam_w[p] * RANGE_V / 127.*.5;
        // else
        //   f_ray[p] = data->beam_w[p] * RANGE_V2 / 127.*.5;
        PolarScan<double>& scan = vol_v->at(el_num);
        scan.row(az_idx) = ms;
        scan.elevations_real(az_idx) = beam.beam_info.elevation;
        scan.azimuths_real(az_idx) = beam.beam_info.azimuth;
    }

    if (vol_w && beam.has_w()) // Spread - Sigma V
    {
        // range variabilita Sigma V - Spread velocità
        const double range_sig_v = 10.;
        // Convert to m/s 
        Eigen::VectorXd ms(max_range);
        for (unsigned i = 0; i < max_range; ++i)
            ms[i] = beam.data_w[i] * range_sig_v / 255.0;

        PolarScan<double>& scan = vol_w->at(el_num);
        scan.row(az_idx) = ms;
        scan.elevations_real(az_idx) = beam.beam_info.elevation;
        scan.azimuths_real(az_idx) = beam.beam_info.azimuth;
    }
}

}
}
