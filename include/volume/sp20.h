#ifndef ARCHIVIATORE_VOLUME_SP20_CLASS_H
#define ARCHIVIATORE_VOLUME_SP20_CLASS_H

#include "volume.h"
#include "volume/loader.h"

namespace cumbac {
namespace volume {

namespace sp20 {
struct Beam;
}

struct SP20Loader : public Loader
{
    Volume<double>* vol_z;
    Volume<double>* vol_d;
    Volume<double>* vol_v;
    Volume<double>* vol_w;

    SP20Loader(const Site& site, bool medium=false, bool clean=false, unsigned max_bin=0)
        : Loader(site, medium, clean, max_bin), vol_z(0), vol_d(0), vol_v(0), vol_w(0)
    {
    }

    void load(const std::string& pathname);

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_size, double cell_size);

    void beam_to_volumes(const sp20::Beam& beam, unsigned beam_size, unsigned el_num);
};

}
}

#endif
