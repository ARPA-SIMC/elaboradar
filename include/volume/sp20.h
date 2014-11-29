#ifndef ARCHIVIATORE_VOLUME_SP20_CLASS_H
#define ARCHIVIATORE_VOLUME_SP20_CLASS_H

#include "elaboradar/volume.h"
#include "elaboradar/azimuthmap.h"

namespace elaboradar {
namespace volume {

namespace sp20 {
struct Beam;
}

struct SP20Loader
{
    std::vector<NonuniformAzimuthMap> azimuth_maps;
    bool clean;
    Scans<double>* vol_z = 0;
    Scans<double>* vol_d = 0;
    Scans<double>* vol_v = 0;
    Scans<double>* vol_w = 0;

    void load(const std::string& pathname);

    void beam_to_volumes(const sp20::Beam& beam, unsigned az_idx, unsigned beam_size, unsigned el_num);
};

}
}

#endif
