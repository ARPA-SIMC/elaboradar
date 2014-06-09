#ifndef ARCHIVIATORE_VOLUME_LOADER_CLASS_H
#define ARCHIVIATORE_VOLUME_LOADER_CLASS_H

#include <volume.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#undef FATT_MOLT_AZ

namespace cumbac {
struct Site;

namespace volume {

static const double FATT_MOLT_AZ = (double) 360./(double)4096.;

// Indice per un beam dato il suo azimut
static inline int azimut_index_MDB(short az)
{
  float azimut = az * FATT_MOLT_AZ / .9;
  if (azimut - (int)azimut <= .5)
    return (int)azimut % 400;
  else
    return ((int)azimut + 1) % 400;
}


struct LoadLogEntry
{
    double theta;
    double alpha;

    LoadLogEntry(double theta, double alpha)
        : theta(theta), alpha(alpha)
    {
    }

    bool operator==(const LoadLogEntry& e) const
    {
        return theta == e.theta && alpha == e.alpha;
    }
};

struct LoadLog : public std::vector<LoadLogEntry>
{
    void log(double theta, double alpha)
    {
        push_back(LoadLogEntry(theta, alpha));
    }
    void print(FILE* out) const;
};

struct BeamInfo
{
    /// Real beam elevation in degrees
    LoadLog load_log;
    double elevation;
};

struct PolarScanLoadInfo
{
    std::vector<volume::BeamInfo> beam_info;

    PolarScanLoadInfo(unsigned beam_count)
    {
        beam_info.resize(beam_count);
    }

    /// Return the number of beams that have been filled with data while loading
    unsigned count_rays_filled() const
    {
        unsigned count = 0;
        for (std::vector<volume::BeamInfo>::const_iterator i = beam_info.begin(); i != beam_info.end(); ++i)
            if (!i->load_log.empty())
                ++count;
        return count;
    }

    /// Return the load log for the given beam
    const volume::LoadLog& get_beam_load_log(unsigned az) const
    {
        return beam_info[az].load_log;
    }

    inline double get_elevation(unsigned az) const
    {
        return beam_info[az].elevation;
    }

    inline double get_elevation_rad(unsigned az) const
    {
        return beam_info[az].elevation * M_PI / 180.;
    }
};

struct LoadInfo
{
    std::vector<PolarScanLoadInfo> scans;
    std::string filename;
    // Acquisition date
    time_t acq_date;
    bool declutter_rsp; // ?

    LoadInfo()
        : declutter_rsp(false)
    {
    }

    const PolarScanLoadInfo& scan(unsigned idx) const
    {
        return scans[idx];
    }

    void make_scan(unsigned idx, unsigned beam_count)
    {
        while (idx >= scans.size())
            scans.push_back(PolarScanLoadInfo(beam_count));
    }
};

// Base class for volume loaders
struct Loader
{
    const Site& site;
    bool medium;
    bool clean;
    std::vector<double> elev_array;

    bool coherent_loader;

    /**
     * If this is greather than zero, truncate each beam to this number of
     * samples
     */
    unsigned max_bin;

    /// If set to non-zero, it will be filled with volume load information
    LoadInfo* load_info;

    Loader(const Site& site, bool medium=false, bool clean=false, unsigned max_bin=0);

    /**
     * Compute the vol_pol index of an elevation angle
     * @returns -1 if no suitable index was found, else the index
     */
    int elevation_index(double elevation) const;

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_size)
    {
        if (load_info) load_info->make_scan(idx, 400);
    }

    template<typename T>
    void fill_beam(PolarScan<T>& scan, int el_num, double theta, double alpha, unsigned size, const T* data)
    {
        int alfa = alpha / FATT_MOLT_AZ;
        if (alfa >= 4096) return;

        int az_num = azimut_index_MDB(alfa);

/*        
        if (az_num == 0)
        {
            printf("fbeam ϑ%f→%d α%f→%d %u", theta, el_num, alpha, az_num, size);
            for (unsigned i = 0; i < 8; ++i)
                printf(" %f\t", data[i]);
            printf("\n");
        }
*/        
// if informations on beam azimuths has been recorded 
// check if the current beam is closer to the nominal beam azimtuh
// if no informations on beam azimuths has been recorded 
// always update to ensure homogeneous beam data against different radar quantities (volumes)

	if(coherent_loader)
	{
	 	unsigned overlap = std::min(scan.beam_size, size);
		if(load_info) 
		{
			PolarScanLoadInfo* li = 0;
	        	li = &(load_info->scans[el_num]);
			if(li->beam_info[az_num].load_log.size())
			{
				double delta_alpha;
				double delta_alpha_old=361.;
				for(unsigned i=0;i<li->beam_info[az_num].load_log.size();i++)
				{
					delta_alpha=std::fabs(az_num*0.9-li->beam_info[az_num].load_log[0].alpha);
					if(delta_alpha<delta_alpha_old) delta_alpha_old=delta_alpha;
				}
				delta_alpha=std::fabs(az_num*0.9-alpha);
				if(delta_alpha<delta_alpha_old)
				{
					for (unsigned i = 0; i < overlap; ++i) scan.set(az_num,i,data[i]);
				}			
			}
			else for (unsigned i = 0; i < overlap; ++i) scan.set(az_num,i,data[i]);
			li->beam_info[az_num].load_log.log(theta,alpha);
		}
		else for (unsigned i = 0; i < overlap; ++i) scan.set(az_num,i,data[i]);
	}
	else  // Old code  										
	{
        	merge_beam(scan, el_num, az_num, theta, alpha, size, data);
	        if(az_num*0.9 - alpha < 0.)
	        {
        	    int new_az_num = (az_num + 1) % 400;
	            if (new_az_num != az_num)
        	        merge_beam(scan, el_num, new_az_num, theta, alpha, size, data);
	        }
	        else if(az_num*0.9 - alpha > 0.)
        	{
	            int new_az_num = (az_num -1+400) %400;
        	    if (new_az_num != az_num)
                	merge_beam(scan, el_num, new_az_num, theta, alpha, size, data);
	        }
	}
	
    }

    template<typename T>
    void merge_beam(PolarScan<T>& scan, int el_num, int az_num, double theta, double alpha, unsigned size, const T* dati)
    {
        PolarScanLoadInfo* li = 0;
        if (load_info) li = &(load_info->scans[el_num]);

        if (li) li->beam_info[az_num].load_log.log(theta, alpha);

        unsigned overlap = std::min(scan.beam_size, size);
        for (unsigned i = 0; i < overlap; ++i)
            if (scan.get(az_num, i) < dati[i])
                scan.set(az_num, i, dati[i]);

        if (li) li->beam_info[az_num].elevation = theta;
    }
};

}
}

#endif
