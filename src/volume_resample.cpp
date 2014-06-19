#include "volume/resample.h"
#include <cmath>

using namespace std;

namespace cumbac {

template<typename T>
void volume_resample(const Volume<T>& src, Volume<T>& dst,
       std::function<void(const PolarScan<T>&, double, PolarScan<T>&, unsigned)> merger)
{
    if (src.beam_count < dst.beam_count)
        throw std::runtime_error("volume_resample currently only work for resampling to smaller volumes");

    // Copy quantity information
    dst.quantity = src.quantity;

    // Do the merge
    for (unsigned i = 0; i < dst.beam_count; ++i)
    {
        // Fractional index specifying the precise position in src that should go to dst
        double src_idx = (double)i * (double)src.beam_count / (double)dst.beam_count;

        // Merge this beam on all scans
        for (unsigned iel = 0; iel < src.size(); ++iel)
        {
            const PolarScan<T>& src_scan = src.scan(iel);
            PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);

            merger(src_scan, src_idx, dst_scan, i);
        }
    }
}

template void volume_resample<double>(const Volume<double>& src, Volume<double>& dst,
        std::function<void(const PolarScan<double>&, double, PolarScan<double>&, unsigned)>);
template void volume_resample<unsigned char>(const Volume<unsigned char>& src, Volume<unsigned char>& dst,
        std::function<void(const PolarScan<unsigned char>&, double, PolarScan<unsigned char>&, unsigned)>);

}
