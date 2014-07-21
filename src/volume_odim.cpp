#include "volume/odim.h"
#include "logging.h"
#include <radarlib/radar.hpp>
#include "utils.h"
#include "volume_cleaner.h"
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

}


namespace cumbac {
namespace volume {

void ODIMLoader::request_quantity(const std::string& name, Scans<double>* volume)
{
    to_load.insert(make_pair(name, volume));
}

void ODIMLoader::make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double size_cell)
{
    if (azimuth_maps.size() <= idx) azimuth_maps.resize(idx + 1);
    for (auto& i: to_load)
        i.second->make_scan(idx, beam_count, beam_size, elev_array[idx], size_cell);
}

void ODIMLoader::load(const std::string& pathname)
{
    LOG_CATEGORY("radar.io");

    namespace odim = OdimH5v21;
    using namespace Radar;
    using namespace std;

    shared_ptr<LoadInfo> load_info = make_shared<LoadInfo>();
    load_info->filename = pathname;

    unique_ptr<odim::OdimFactory> factory(new odim::OdimFactory());
    unique_ptr<odim::PolarVolume> volume(factory->openPolarVolume(pathname));

    load_info->acq_date = volume->getDateTime();

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
        unique_ptr<odim::PolarScan> scan(volume->getScan(src_elev));
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

        // Get and validate the azimuth angles for this scan
        std::vector<odim::AZAngles> azangles = scan->getAzimuthAngles();
        int rpm_sign = scan->getDirection();

        std::vector<double> elevation_angles = scan->getElevationAngles();

        unsigned beam_count = int_to_unsigned(scan->getNumRays(), "number of rays");
        if (azangles.size() != beam_count)
        {
            LOG_ERROR("elevation %f has %zd azimuth angles and %d rays", elevation, azangles.size(), beam_count);
            throw runtime_error("mismatch between number of azumuth angles and number of rays");
        }

        unsigned beam_size = int_to_unsigned(scan->getNumBins(), "beam size");

        int el_num = elevation_index(elevation);
        if (el_num < 0) continue;

        // Create PolarScan objects for this elevation
        make_scan(el_num, beam_count, beam_size, range_scale);

        // Fill in the azimuth map for this elevation
        for (unsigned src_az = 0; src_az < beam_count; ++src_az)
            azimuth_maps[el_num].add(azangles[src_az].averagedAngle(rpm_sign), src_az);

        // Read all quantities that have been requested
        for (auto& todo : to_load)
        {
            const string& name = todo.first;
            Scans<double>& target = *todo.second;
            PolarScan<double>& vol_pol_scan = target.at(el_num);

            // Pick the best quantity among the ones available
            if (!scan->hasQuantityData(name))
            {
                LOG_WARN("no %s found for elevation angle %f: skipping", name.c_str(), elevation);
                continue;
            }

            unique_ptr<odim::PolarScanData> data(scan->getQuantityData(odim::PRODUCT_QUANTITY_DBZH));

            // Fill variable metadata
            target.quantity.name = name;
            target.quantity.nodata = data->getNodata() * data->getGain() + data->getOffset();
            target.quantity.undetect = data->getUndetect() * data->getGain() + data->getOffset();
            target.quantity.gain = data->getGain();
            target.quantity.offset = data->getOffset();

            // Read actual data from ODIM
            odim::RayMatrix<double> matrix;
            matrix.resize(beam_count, beam_size);
            data->readTranslatedData(matrix);

            for (unsigned src_az = 0; src_az < beam_count; ++src_az)
            {
                Eigen::VectorXd beam(beam_size);

                for (unsigned i = 0; i < beam_size; ++i)
                    beam(i) = matrix.elem(src_az, i);

                vol_pol_scan.row(src_az) = beam;
                vol_pol_scan.elevations_real(src_az) = elevation_angles[src_az];
            }

        }
    }

    for (auto& i: to_load)
        i.second->load_info = load_info;
}

#if 0
bool ODIMLoader::load(const std::string& pathname, const std::string& quantity)
{
	LOG_CATEGORY("radar.io");
	namespace odim = OdimH5v21;
	using namespace Radar;
	using namespace std;

    shared_ptr<LoadInfo> load_info = make_shared<LoadInfo>();
	load_info->filename = pathname;

	unique_ptr<odim::OdimFactory> factory(new odim::OdimFactory());
	unique_ptr<odim::PolarVolume> volume(factory->openPolarVolume(pathname));

	load_info->acq_date = volume->getDateTime();

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
		unique_ptr<odim::PolarScan> scan(volume->getScan(src_elev));
		double elevation = scan->getEAngle();

		// Read and and validate resolution information
		if (src_elev == 0) range_scale = scan->getRangeScale();
		else 
		{
			double rs = scan->getRangeScale();
			if (rs != range_scale)
			{
				LOG_ERROR("scan %d (elevation %f) has rangeScale %f that is different from %f in the previous scans",
				src_elev, elevation, rs, range_scale);
				throw runtime_error("rangeScale mismatch");
			}
		}
		// Pick the best quantity among the ones available
		unique_ptr<odim::PolarScanData> data;
		if (scan->hasQuantityData(quantity))
		{
			data.reset(scan->getQuantityData(quantity));
			vol_z->quantity.name=quantity;
			vol_z->quantity.nodata=data->getNodata() * data->getGain()+data->getOffset();
			vol_z->quantity.undetect=data->getUndetect() * data->getGain()+data->getOffset();
		}
		else
		{
			LOG_WARN("no %s found for elevation angle %f", quantity.c_str(),elevation);
			return false;
		}
		// Get and validate the azimuth angles for this scan
		std::vector<odim::AZAngles> azangles = scan->getAzimuthAngles();
		int rpm_sign = scan->getDirection();	///! USELESS??		
		unsigned beam_count = int_to_unsigned(data->getNumRays(), "number of rays");
		
		if (azangles.size() != beam_count)
		{
			LOG_ERROR("elevation %f has %zd azimuth angles and %d rays", elevation, azangles.size(), beam_count);
			throw runtime_error("mismatch between number of azumuth angles and number of rays");
		}

		std::vector<double> elevation_angles = scan->getElevationAngles();
		unsigned beam_size = int_to_unsigned(data->getNumBins(), "beam size");

		// Read all scan beam data
		odim::RayMatrix<double> matrix;
		matrix.resize(beam_count, beam_size);
		data->readTranslatedData(matrix); // 1)

		int el_num = elevation_index(elevation);
		if (el_num < 0) continue;
		make_scan(el_num, beam_count, beam_size, range_scale);
		PolarScan<double>& vol_pol_scan = vol_z->at(el_num); // PolarScan<double>& vol_pol_scan = vol_z->scan(el_num);
									// vol_z ora è un Scans* e non possiede il metodo scan(unsigned)

        //vol_pol_scan.elevation = elevation;
	std::vector<bool> angles_seen(400, false);
	for (unsigned src_az = 0; src_az < beam_count; ++src_az)
	{
		azimuth_maps[el_num].add(azangles[src_az].averagedAngle(rpm_sign), src_az);
		Eigen::VectorXd beam(beam_size);
/*
            if (clean) {
		//std::cout<<"PASSO PER IL CLEANER - Raggio "<<src_az<<std::endl;
              BeamCleaner<double> cleaner(bin_wind_magic_number, Z_missing, W_threshold, V_missing);
              unique_ptr<Beams<double>> b(new Beams<double>);
              for (unsigned i = 0; i < beam_size; ++i){
                         b->data_z[i]=matrix.elem(src_az,i);         
                         b->data_v[i]=VRAD_matrix.elem(src_az,i);         
                         b->data_w[i]=WRAD_matrix.elem(src_az,i);         
               }
               vector <bool> cleaned(beam_size,false);
               cleaner.clean_beams(*b,beam_size,cleaned);
               for (unsigned i = 0; i < beam_size; ++i)
                 beam(i) = b->data_z[i];
            } else*/ {
             // Convert back to bytes, to fit into vol_pol as it is now
               for (unsigned i = 0; i < beam_size; ++i){
                 // FIXME: QUESTO PEZZO DI CODICE E' STATO INSERITO PER EMULARE LA CONVERSIONE ELDES IN FORMATO SP20
                 // DEVE ESSERE RIMOSSO A FINE LAVORO E RIATTIVATA QUESTA LINEA DI CODICE ORA COMMENTATA
                 // beam[i] = DBtoBYTE(matrix.elem(src_az, i));
                 beam(i) = matrix.elem(src_az, i);
                 //beam[i] = BYTEtoDB(eldes_counter_to_db(matrix.elem(src_az, i)));
              }
            }

            vol_pol_scan.row(src_az) = beam;
            vol_pol_scan.elevations_real(src_az) = elevation_angles[src_az];
        }

    }


////////////////////////////////////////////////////////////////////////////////7

        /*
        for (unsigned src_az = 0; src_az < beam_count; ++src_az)
        {
            // FIXME: reproduce a bad truncation from the eldes sp20 converte
            //double azimut = azangles[src_az].averagedAngle(rpm_sign);
            double azimut = eldes_converter_azimut(azangles[src_az].start, azangles[src_az].stop);
            double* beam = new double[beam_size];

            for (unsigned i = 0; i < beam_size; ++i){
              // FIXME: QUESTO PEZZO DI CODICE E' STATO INSERITO PER EMULARE LA CONVERSIONE ELDES IN FORMATO SP20
              // DEVE ESSERE RIMOSSO A FINE LAVORO E RIATTIVATA QUESTA LINEA DI CODICE ORA COMMENTATA
              // beam[i] = DBtoBYTE(matrix.elem(src_az, i));
              beam[i] = matrix.elem(src_az, i);
              //beam[i] = BYTEtoDB(eldes_counter_to_db(matrix.elem(src_az, i)));
            }
            fill_beam(vol_pol_scan, el_num, elevation_angles[src_az], azimut, beam_size, beam);
            delete[] beam;
        }
	*/

    if (vol_z) vol_z->load_info = load_info;

    return true;
}
#endif


}
}
