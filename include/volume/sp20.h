#ifndef ARCHIVIATORE_VOLUME_SP20_CLASS_H
#define ARCHIVIATORE_VOLUME_SP20_CLASS_H

#include "volume.h"
#include "volume/loader.h"
#include "volume/azimuthmap.h"

namespace elaboradar {
namespace volume {

namespace sp20 {
struct Beam;
}

struct SP20Loader : public Loader
{
    std::vector<NonuniformAzimuthMap> azimuth_maps;
    Scans<double>* vol_z;
    Scans<double>* vol_d;
    Scans<double>* vol_v;
    Scans<double>* vol_w;

    SP20Loader(const Site& site, bool medium=false, bool clean=false)
        : Loader(site, medium, clean), vol_z(0), vol_d(0), vol_v(0), vol_w(0)
    {
    }

    void load(const std::string& pathname);

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double cell_size);

    void beam_to_volumes(const sp20::Beam& beam, unsigned az_idx, unsigned beam_size, unsigned el_num);
};

}
}

#endif
