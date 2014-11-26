#include "volume/loader.h"
#include "logging.h"
#include "site.h"

namespace elaboradar {
namespace volume {

Loader::Loader(const Site& site, bool medium)
    : site(site), elev_array(site.get_elev_array(medium))
{
}

int Loader::elevation_index(double elevation) const
{
    for (unsigned i=0; i < elev_array.size(); ++i)
        if (elevation >= (elev_array[i]-0.5) && elevation < (elev_array[i]+0.5))
            return i;
    return -1;
}


}
}
