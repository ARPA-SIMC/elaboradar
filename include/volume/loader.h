#ifndef ARCHIVIATORE_VOLUME_LOADER_CLASS_H
#define ARCHIVIATORE_VOLUME_LOADER_CLASS_H

#include <volume.h>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <cmath>

#undef FATT_MOLT_AZ

namespace elaboradar {
struct Site;

namespace volume {

static const double FATT_MOLT_AZ = (double) 360./(double)4096.;

// Base class for volume loaders
struct Loader
{
    const Site& site;
    std::vector<double> elev_array;

    Loader(const Site& site, bool medium=false);

    /**
     * Compute the vol_pol index of an elevation angle
     * @returns -1 if no suitable index was found, else the index
     */
    int elevation_index(double elevation) const;
};

}
}

#endif
