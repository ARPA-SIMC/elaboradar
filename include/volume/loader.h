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
static inline int azimut_index_MDB(short az, unsigned beam_count)
{
  float azimut = az * FATT_MOLT_AZ / (360. / beam_count);
  if (azimut - (int)azimut <= .5)
    return (int)azimut % beam_count;
  else
    return ((int)azimut + 1) % beam_count;
}

template<typename T>
struct Scans : public std::vector<PolarScan<T>>
{
    Variable<T> quantity;

    // Create or reuse a scan at position idx, with the given beam size
    PolarScan<T>& make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double elevation, double cell_size)
    {
        if (idx < this->size())
        {
            if (beam_count != (*this)[idx].beam_count)
            {
                LOG_CATEGORY("radar.io");
                LOG_ERROR("make_scan(idx=%u, beam_count=%u, beam_size=%u) called, but the scan already existed with beam_count=%u", idx, beam_count, beam_size, (*this)[idx].beam_count);
                throw std::runtime_error("beam_size mismatch");
            }
            if (beam_size != (*this)[idx].beam_size)
            {
                LOG_CATEGORY("radar.io");
                LOG_ERROR("make_scan(idx=%u, beam_count=%u, beam_size=%u) called, but the scan already existed with beam_size=%u", idx, beam_count, beam_size, (*this)[idx].beam_size);
                throw std::runtime_error("beam_size mismatch");
            }
        } else {
            // If some elevation has been skipped, fill in the gap
            if (idx > this->size())
            {
                if (this->empty())
                    this->push_back(PolarScan<T>(beam_count, beam_size));
                while (this->size() < idx)
                    this->push_back(PolarScan<T>(beam_count, this->back().beam_size));
            }

            // Add the new polar scan
            this->push_back(PolarScan<T>(beam_count, beam_size));
            this->back().elevation = elevation;
            this->back().cell_size = cell_size;
        }

        // Return it
        return (*this)[idx];
    }
};

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
        // FIXME: hardcoding 450 to have enough space to test new SP20 loading.
        // Long term plan: get rid of this, and let volume_resample mergers do
        // accounting if they feel like it.
        if (load_info) load_info->make_scan(idx, 460);
    }
};

}
}

#endif
