#include "volume/resample.h"
#include <cmath>

using namespace std;

namespace cumbac {

template<typename T>
void polarscan_resample(const PolarScan<T>& src, PolarScan<T>& dst,
        std::function<void(const PolarScan<T>&, double, PolarScan<T>&, unsigned)> merger)
{
    if (src.beam_count < dst.beam_count)
        throw std::runtime_error("polarscan_resample currently only work for resampling to smaller volumes");

    // Do the merge
    for (unsigned i = 0; i < dst.beam_count; ++i)
    {
        // Fractional index specifying the precise position in src that should go to dst
        double src_idx = (double)i * (double)src.beam_count / (double)dst.beam_count;

        merger(src, src_idx, dst, i);
    }
}

template<typename T>
void merger_max_of_closest(const PolarScan<T>& src, double src_idx, PolarScan<T>& dst, unsigned dst_idx)
{
    double rounded_idx = round(src_idx);

    // Copy the closest beam
    unsigned idx = (unsigned)rounded_idx % src.beam_count;
    dst.row(dst_idx) = src.row(idx);

    // Copy the previous or next beam in case they overlap
    if (rounded_idx < src_idx)
    {
        unsigned idx = ((unsigned)rounded_idx + 1) % src.beam_count;
        for (unsigned i = 0; i < dst.beam_size; ++i)
            if (dst(dst_idx, i) < src(idx, i))
                dst(dst_idx, i) = src(idx, i);
    } else {
        unsigned idx = ((unsigned)rounded_idx - 1) % src.beam_count;
        for (unsigned i = 0; i < dst.beam_size; ++i)
            if (dst(dst_idx, i) < src(idx, i))
                dst(dst_idx, i) = src(idx, i);
    }
}


template void polarscan_resample<double>(const PolarScan<double>& src, PolarScan<double>& dst,
        std::function<void(const PolarScan<double>&, double, PolarScan<double>&, unsigned)>);
template void polarscan_resample<unsigned char>(const PolarScan<unsigned char>& src, PolarScan<unsigned char>& dst,
        std::function<void(const PolarScan<unsigned char>&, double, PolarScan<unsigned char>&, unsigned)>);

template void merger_max_of_closest<double>(const PolarScan<double>&, double, PolarScan<double>&, unsigned);
template void merger_max_of_closest<unsigned char>(const PolarScan<unsigned char>&, double, PolarScan<unsigned char>&, unsigned);

}
