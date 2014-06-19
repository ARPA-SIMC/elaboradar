#ifndef ARCHIVIATORE_VOLUME_RESAMPLE_H
#define ARCHIVIATORE_VOLUME_RESAMPLE_H

#include "volume.h"
#include <functional>
#include <cmath>

namespace cumbac {

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
        std::function<void(const PolarScan<T>&, double, PolarScan<T>&, unsigned)> merger);

}

#endif
