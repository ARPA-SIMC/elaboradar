#ifndef RADARELAB_SP20_H
#define RADARELAB_SP20_H

#include <radarelab/volume.h>

namespace radarelab {
namespace volume {

namespace sp20 {
struct Beam;
}

struct SP20Loader
{
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
