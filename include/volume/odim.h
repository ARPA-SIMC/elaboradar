#ifndef ARCHIVIATORE_VOLUME_ODIM_CLASS_H
#define ARCHIVIATORE_VOLUME_ODIM_CLASS_H

#include "volume.h"
#include "volume/loader.h"
#include "volume/azimuthmap.h"

namespace cumbac {
namespace volume {

struct ODIMLoader : public volume::Loader
{
    std::vector<NonuniformAzimuthMap> azimuth_maps;
    Scans<double>* vol_z;

    ODIMLoader(const Site& site, bool medium=false, bool clean=false, unsigned max_bin=0)
        : Loader(site, medium, clean, max_bin), vol_z(0)
    {
    }

    void load(const std::string& pathname);
    bool load(const std::string& pathname,const std::string& quantity);

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double size_cell);
};

}
}

#endif

