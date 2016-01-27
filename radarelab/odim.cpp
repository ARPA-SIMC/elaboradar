#include "odim.h"
#include "logging.h"
#include "utils.h"
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

}


namespace radarelab {
namespace volume {

void ODIMLoader::request_quantity(const std::string& name, Scans<double>* volume)
{
    to_load.insert(make_pair(name, volume));
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

    double range_scale = 0;

    // Iterate all scans
    unsigned scan_count = int_to_unsigned(volume->getScanCount(), "scan count");
    double old_elevation = -1000.;
    for (unsigned src_elev = 0; src_elev < scan_count; ++src_elev)
    {
        unique_ptr<odim::PolarScan> scan(volume->getScan(src_elev));
        double elevation = scan->getEAngle();
	Available_Elevations.push_back(elevation);
        // FIXME: please add a comment on why this is needed: are there faulty
        // ODIM files out there which repeat elevations? Is it correct that
        // only the first elevation is used and not, say, the last one?
        if( elevation == old_elevation ) continue;
        old_elevation=elevation;
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

        // Read all quantities that have been requested
        for (auto& todo : to_load)
        {
            // Create azimuth maps and PolarScan objects for this elevation
            const string& name = todo.first;
            Scans<double>& target = *todo.second;

            // Pick the best quantity among the ones available
            if (!scan->hasQuantityData(name))
            {
                LOG_WARN("no %s found for elevation angle %f: skipping", name.c_str(), elevation);
                continue;
            }
            PolarScan<double>& vol_pol_scan = target.append_scan(beam_count, beam_size, elevation, range_scale);

            unique_ptr<odim::PolarScanData> data(scan->getQuantityData(name));

            // Fill/overwrite variable metadata at volume level
            target.quantity = name;
	    target.h_radar = volume->getAltitude()/1000.;
            // Fill variable metadata at scan level
            vol_pol_scan.nodata = data->getNodata() * data->getGain() + data->getOffset();
            vol_pol_scan.undetect = data->getUndetect() * data->getGain() + data->getOffset();
            vol_pol_scan.gain = data->getGain();
            vol_pol_scan.offset = data->getOffset();

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
                vol_pol_scan.azimuths_real(src_az) = azangles[src_az].averagedAngle(rpm_sign);
           }
        }
    }
    for (auto& i: to_load)
	{ 
        i.second->load_info = load_info;
	}
}

void ODIMStorer::store(const std::string& pathname)
{
//    LOG_CATEGORY("radar.io");

    namespace odim = OdimH5v21;
    using namespace Radar;
    using namespace std;

    shared_ptr<LoadInfo> load_info = make_shared<LoadInfo>();
    load_info->filename = pathname;

    unique_ptr<odim::OdimFactory> factory(new odim::OdimFactory());
    unique_ptr<odim::PolarVolume> volume(factory->openPolarVolume(pathname));

// cout<<"aperto file"<<endl;
//    unsigned max_elev=0;
    for(unsigned i=0;i<to_store_int.size();i++)
	for(unsigned j=0;j<to_store_int[i]->size();j++)
	{
		vector<odim::PolarScan*> scans;
		scans = volume->getScans(to_store_int[i]->scan(j).elevation,0.1);
		shared_ptr<odim::PolarScan> scan;
		if(scans.size()) scan.reset(scans[0]);
		else
		{
			scan.reset(volume->createScan());
			scan->setEAngle(to_store_int[i]->scan(j).elevation);
		}
		unique_ptr<odim::PolarScanData> data(scan->createQuantityData(to_store_int[i]->quantity));
		data->setGain((double)to_store_int[i]->scan(j).gain);
		data->setOffset((double)to_store_int[i]->scan(j).offset);
		data->setNodata((double)to_store_int[i]->scan(j).nodata);
		data->setUndetect((double)to_store_int[i]->scan(j).undetect);
		odim::RayMatrix<unsigned short> matrix;
        	matrix.resize(to_store_int[i]->scan(j).beam_count,to_store_int[i]->scan(j).beam_size);
		for(unsigned ii=0;ii<to_store_int[i]->scan(j).beam_count;ii++)
			for(unsigned jj=0;jj<to_store_int[i]->scan(j).beam_size;jj++)
				matrix.elem(ii,jj) = to_store_int[i]->scan(j)(ii,jj);
		data->writeData(matrix);
	}
//cout<<"vado coi double"<<endl;
    for(unsigned i=0;i<to_store_fp.size();i++)
	for(unsigned j=0;j<to_store_fp[i]->size();j++)
	{
//cout<<"vol "<<i<<"scan "<<j<<endl;
		vector<odim::PolarScan*> scans;
		scans = volume->getScans(to_store_fp[i]->scan(j).elevation,0.1);
		shared_ptr<odim::PolarScan> scan;
		if(scans.size()) scan.reset(scans[0]);
		else
		{
			scan.reset(volume->createScan());
			scan->setEAngle(to_store_fp[i]->scan(j).elevation);
		}
//cout<<"settato puntatore a scan"<<endl;
		unique_ptr<odim::PolarScanData> data(scan->createQuantityData(to_store_fp[i]->quantity));
//cout<<"settato puntatore a data "<<to_store_fp[i]->quantity<<endl;
		data->setGain((double)to_store_fp[i]->scan(j).gain);
		data->setOffset((double)to_store_fp[i]->scan(j).offset);
		data->setNodata((double)to_store_fp[i]->scan(j).nodata);
		data->setUndetect((double)to_store_fp[i]->scan(j).undetect);
//cout<<"settati metadati"<<endl;
		odim::RayMatrix<float> matrix;
        	matrix.resize(to_store_fp[i]->scan(j).beam_count,to_store_fp[i]->scan(j).beam_size);
//cout<<"settate dimensioni matrice"<<endl;
		for(unsigned ii=0;ii<to_store_fp[i]->scan(j).beam_count;ii++)
			for(unsigned jj=0;jj<to_store_fp[i]->scan(j).beam_size;jj++)
				matrix.elem(ii,jj) = to_store_fp[i]->scan(j)(ii,jj);
//cout<<"scrivo matrice"<<endl;
		data->writeAndTranslate(matrix,(float)data->getOffset(),(float)data->getGain(),H5::PredType::NATIVE_UINT16);
//cout<<"scritto"<<endl;
	}
}


}	// namespace volume
}	// namespace radarelab
