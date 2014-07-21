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

struct LoadInfo
{
    std::string filename;
    // Acquisition date
    time_t acq_date;
    bool declutter_rsp; // ?

    LoadInfo()
        : declutter_rsp(false)
    {
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
};

}
}

#endif
