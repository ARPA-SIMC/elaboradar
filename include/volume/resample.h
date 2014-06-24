#ifndef ARCHIVIATORE_VOLUME_RESAMPLE_H
#define ARCHIVIATORE_VOLUME_RESAMPLE_H

#include "volume.h"
#include "volume/loader.h"
#include <functional>
#include <cmath>

namespace cumbac {

template<typename T>
void polarscan_resample(const PolarScan<T>& src, PolarScan<T>& dst,
        std::function<void(const PolarScan<T>&, double, PolarScan<T>&, unsigned)> merger);

/**
 * Fill dst with data from src, coping with the two volumes having a different
 * number of beams per scan.
 *
 * Merger is the function used to merge beams from src into dst. It takes the
 * source PolarScan, the destination PolarScan and a vector with the indices of
 * the beams of src that need to be used.
 */
template<typename T>
void volume_resample(const Volume<T>& src, Volume<T>& dst,
        std::function<void(const PolarScan<T>&, double, PolarScan<T>&, unsigned)> merger)
{
    // Copy quantity information
    dst.quantity = src.quantity;

    for (unsigned iel = 0; iel < src.size(); ++iel)
    {
        const PolarScan<T>& src_scan = src.scan(iel);
        PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
        polarscan_resample(src_scan, dst_scan, merger);
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
template<typename T>
void volume_resample(const volume::Scans<T>& src, Volume<T>& dst,
        std::function<void(const PolarScan<T>&, double, PolarScan<T>&, unsigned)> merger)
{
    // Copy quantity information
    dst.quantity = src.quantity;

    for (unsigned iel = 0; iel < src.size(); ++iel)
    {
        const PolarScan<T>& src_scan = src.at(iel);
        PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
        polarscan_resample(src_scan, dst_scan, merger);
    }
}

template<typename T>
void merger_max_of_closest(const PolarScan<T>& src, double src_idx, PolarScan<T>& dst, unsigned dst_idx);

}

#endif
