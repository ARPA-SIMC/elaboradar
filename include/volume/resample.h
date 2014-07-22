#ifndef ARCHIVIATORE_VOLUME_RESAMPLE_H
#define ARCHIVIATORE_VOLUME_RESAMPLE_H

#include "volume.h"
#include "volume/loader.h"
#include "volume/azimuthmap.h"
#include <functional>
#include <cmath>

namespace cumbac {
namespace volume {

template<typename T>
void polarscan_resample(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst,
        std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger);

/**
 * Fill dst with data from src, coping with the two volumes having a different
 * number of beams per scan.
 *
 * Merger is the function used to merge beams from src into dst. It takes the
 * source PolarScan, the destination PolarScan and a vector with the indices of
 * the beams of src that need to be used.
 */
template<typename T>
void volume_resample(const Volume<T>& src, const AzimuthMap& azmap, Volume<T>& dst,
        std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger)
{
    // Copy volume metadata
    dst.quantity = src.quantity;
    dst.load_info = src.load_info;

    for (unsigned iel = 0; iel < src.size(); ++iel)
    {
        const PolarScan<T>& src_scan = src.scan(iel);
        PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
        polarscan_resample(src_scan, azmap, dst_scan, merger);
    }
}

/**
 * Fill dst with data from src, coping with the two volumes having a different
 * number of beams per scan.
 *
 * Merger is the function used to merge beams from src into dst. It takes the
 * source PolarScan, the destination PolarScan and a vector with the indices of
 * the beams of src that need to be used.
 */
template<typename T, typename AZMAP>
void volume_resample(const volume::Scans<T>& src, const std::vector<AZMAP>& azmaps, Volume<T>& dst,
        std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger)
{
    // Copy volume metadata
    dst.load_info = src.load_info;
    ds.quantity = src.quantity;

    for (unsigned iel = 0; iel < src.size(); ++iel)
    {
        const PolarScan<T>& src_scan = src.at(iel);
        PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
        polarscan_resample(src_scan, azmaps[iel], dst_scan, merger);
        dst_scan.quantity = src_scan.quantity;
    }
}

template<typename T>
void merger_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx);

template<typename T>
void merger_max_of_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx);

}
}

#endif
