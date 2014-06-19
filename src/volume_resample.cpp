#include "volume/resample.h"
#include <cmath>

using namespace std;

namespace cumbac {

template<typename T>
void volume_resample(const Volume<T>& src, Volume<T>& dst,
       std::function<void(const PolarScan<T>&, PolarScan<T>&, std::vector<unsigned>)> merger)
{
    if (src.beam_count < dst.beam_count)
        throw std::runtime_error("volume_resample currently only work for resampling to smaller volumes");

    // Copy quantity information
    dst.quantity = src.quantity;

    // Compute what src beams are used for each dst beam
    vector<vector<unsigned>> merge_map(dst.beam_count);
    for (unsigned i = 0; i < src.beam_count; ++i)
    {
        unsigned dst_idx = round((double)i * (double)dst.beam_count / (double)src.beam_count);
        merge_map[dst_idx].push_back(i);
    }

    // Do the merge
    for (auto src_beams : merge_map)
    {
        if (src_beams.empty())
            throw std::runtime_error("no source beams found");

        // Merge this beam on all scans
        for (unsigned iel = 0; iel < src.size(); ++iel)
        {
            const PolarScan<T>& src_scan = src.scan(iel);
            PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);

            merger(src_scan, dst_scan, src_beams);
        }
    }
}

template void volume_resample<double>(const Volume<double>& src, Volume<double>& dst,
        std::function<void(const PolarScan<double>&, PolarScan<double>&, std::vector<unsigned>)>);
template void volume_resample<unsigned char>(const Volume<unsigned char>& src, Volume<unsigned char>& dst,
        std::function<void(const PolarScan<unsigned char>&, PolarScan<unsigned char>&, std::vector<unsigned>)>);

}
