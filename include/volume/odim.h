#ifndef ARCHIVIATORE_VOLUME_ODIM_CLASS_H
#define ARCHIVIATORE_VOLUME_ODIM_CLASS_H

#include "volume.h"
#include "volume/loader.h"

namespace cumbac {
namespace volume {

struct ODIMLoader : public volume::Loader
{
    Volume<double>* vol_db;

    ODIMLoader(const Site& site, bool medium=false, bool clean=false, unsigned max_bin=0)
        : Loader(site, medium, clean, max_bin), vol_db(0)
    {
    }

    void load(const std::string& pathname);

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_size);
};

}
}

#endif

