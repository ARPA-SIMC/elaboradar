#include "volume/loader.h"
#include "logging.h"
#include "site.h"

namespace cumbac {
namespace volume {

void LoadLog::print(FILE* out) const
{
    if (empty())
    {
        fprintf(out, "no beams loaded\n");
        return;
    }

    for (const_iterator i = begin(); i != end(); ++i)
    {
        if (i != begin())
            fprintf(out, ", ");
        fprintf(out, "ϑ%.2f α%.2f", i->theta, i->alpha);
    }
    fprintf(out, "\n");
}

Loader::Loader(const Site& site, bool medium, bool clean, unsigned max_bin)
    : site(site), medium(medium), clean(clean), elev_array(site.get_elev_array(medium)), max_bin(max_bin), load_info(0), coherent_loader(false)
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
