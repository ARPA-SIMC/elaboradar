#ifndef ARCHIVIATORE_VOLUME_ODIM_CLASS_H
#define ARCHIVIATORE_VOLUME_ODIM_CLASS_H

#include "volume.h"
#include "volume/loader.h"
#include "volume/azimuthmap.h"
#include <string>
#include <map>
#include <vector>

namespace elaboradar {
namespace volume {

struct ODIMLoader : public volume::Loader
{
    std::vector<NonuniformAzimuthMap> azimuth_maps;
    std::map<std::string, Scans<double>*> to_load;

    ODIMLoader(const Site& site, bool medium=false, unsigned max_bin=0)
        : Loader(site, medium, false, max_bin)
    {
    }

    void request_quantity(const std::string& name, Scans<double>* volume);

    void load(const std::string& pathname);

    // Create or reuse a scan at position idx, with the given beam size
    void make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double size_cell);
};

/*
union bit64
{
	double fp64;
	int64_t int64;
};
*/

class ODIMStorer : public volume::Loader
{
public:
    ODIMStorer(const Site& site, bool medium=false, unsigned max_bin=0)
        : Loader(site, medium, false, max_bin)
    {
    }
    std::vector<Volume<double>*> to_store_fp;
    std::vector<Volume<int>*> to_store_int;
    void store_quantity_fp(Volume<double>* vol_fp) {to_store_fp.push_back(vol_fp);}
    void store_quantity_int(Volume<int>* vol_int) {to_store_int.push_back(vol_int);}
    void store(const std::string& pathname);
};

} // volume
} // elaboradars

#endif

