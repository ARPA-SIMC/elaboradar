#ifndef ARCHIVIATORE_VOLUME_RESAMPLE_H
#define ARCHIVIATORE_VOLUME_RESAMPLE_H

#include "volume.h"

namespace cumbac {

/// Fill dst with data from src, coping with the two volumes having a different
/// number of beams per scan
template<typename T>
void volume_resample(const Volume<T>& src, Volume<T>& dst);

}

#endif
