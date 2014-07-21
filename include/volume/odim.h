#ifndef ARCHIVIATORE_VOLUME_ODIM_CLASS_H
#define ARCHIVIATORE_VOLUME_ODIM_CLASS_H

#include "volume.h"
#include "volume/loader.h"
#include "volume/azimuthmap.h"
#include <string>
#include <map>

namespace cumbac {
namespace volume {

struct ODIMLoader : public volume::Loader
{
    std::vector<NonuniformAzimuthMap> azimuth_maps;
    std::map<std::string, Scans<double>*> to_load;

    ODIMLoader(const Site& site, bool medium=false, bool clean=false, unsigned max_bin=0)
        : Loader(site, medium, clean, max_bin)
    {
    }

    void request_quantity(const std::string& name, Scans<double>* volume);

    void load(const std::string& pathname);
    //bool load(const std::string& pathname,const std::string& quantity);

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double size_cell);
};

}
}

#endif

