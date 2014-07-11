#include "volume/resample.h"
#include <cmath>
#include<iostream>

using namespace std;
using namespace cumbac::volume::azimuthmap;

namespace cumbac {
namespace volume {

template<typename T>
void polarscan_resample(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst,
        std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger)
{
    if (src.beam_count < dst.beam_count)
        throw std::runtime_error("polarscan_resample currently only work for resampling to smaller volumes");

    // Do the merge
    for (unsigned i = 0; i < dst.beam_count; ++i)
        merger(src, azmap, dst, i);
}

template<typename T>
void merger_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx)
{
    double dst_azimuth = (double)dst_idx / (double)dst.beam_count * 360.0;
    Position pos = azmap.closest(dst_azimuth);

    //// Copy the closest beam
    dst.row(dst_idx) = src.row(pos.index);
}

template<typename T>
void merger_max_of_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx)
{
    double dst_azimuth = (double)dst_idx / (double)dst.beam_count * 360.0;
    vector<Position> positions = azmap.intersecting(dst_azimuth, 360.0 / dst.beam_count);
    if (positions.empty()) throw std::runtime_error("no source beams found");

    // Copy the first beam
    auto i = positions.cbegin();
    dst.row(dst_idx) = src.row(i->index);

    // Take the maximum of all beam values
    for (++i; i != positions.end(); ++i)
    {
        for (unsigned bi = 0; bi < dst.beam_size; ++bi)
            if (dst(dst_idx, bi) < src(i->index, bi))
                dst(dst_idx, bi) = src(i->index, bi);
    }
}


template void polarscan_resample<double>(const PolarScan<double>& src, const AzimuthMap&, PolarScan<double>& dst,
        std::function<void(const PolarScan<double>&, const AzimuthMap&, PolarScan<double>&, unsigned)>);
template void polarscan_resample<unsigned char>(const PolarScan<unsigned char>& src, const AzimuthMap&, PolarScan<unsigned char>& dst,
        std::function<void(const PolarScan<unsigned char>&, const AzimuthMap&, PolarScan<unsigned char>&, unsigned)>);

template void merger_closest<double>(const PolarScan<double>&, const AzimuthMap&, PolarScan<double>&, unsigned);
template void merger_closest<unsigned char>(const PolarScan<unsigned char>&, const AzimuthMap&, PolarScan<unsigned char>&, unsigned);

template void merger_max_of_closest<double>(const PolarScan<double>&, const AzimuthMap&, PolarScan<double>&, unsigned);
template void merger_max_of_closest<unsigned char>(const PolarScan<unsigned char>&, const AzimuthMap&, PolarScan<unsigned char>&, unsigned);

}
}
