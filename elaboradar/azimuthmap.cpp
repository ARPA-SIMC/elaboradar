#include "azimuthmap.h"

using namespace std;
using namespace elaboradar::azimuthmap;

namespace elaboradar {

namespace azimuthmap {
std::ostream& operator<<(std::ostream& o, const Position& pos)
{
    return o << pos.index << ":" << pos.azimuth;
}
}

inline double angle_distance(double first, double second)
{
    return 180.0 - std::fabs(std::fmod(std::fabs(first - second), 360.0) - 180.0);
}

azimuthmap::Position UniformAzimuthMap::closest(double azimuth) const
{
    unsigned index = (unsigned)round(azimuth * beam_count / 360.0) % beam_count;
    return Position((double)index * 360.0 / (double)beam_count, index);
}

std::vector<azimuthmap::Position> UniformAzimuthMap::intersecting(double azimuth, double amplitude) const
{
    double my_amplitude = 360.0 / (double)beam_count;
    // Angles closer than this amount are considered the same for overlap detection
    static const double precision = 0.000000001;

    double lowest_angle = azimuth - amplitude / 2 - my_amplitude / 2 + precision;
    int lowest_index = round(lowest_angle * beam_count / 360.0);

    double highest_angle = azimuth + amplitude / 2 + my_amplitude / 2 - precision;
    int highest_index = round(highest_angle * beam_count / 360.0);

    std::vector<Position> res;
    for (int i = lowest_index; i <= highest_index; ++i)
    {
        unsigned index = (i + 360) % 360;
        res.push_back(Position(index * 360.0 / beam_count, index));
    }
    return res;
}

template<typename ITER1, typename ITER2>
azimuthmap::Position NonuniformAzimuthMap::closest_of_two(double azimuth, const ITER1& i1, const ITER2& i2) const
{
    if (angle_distance(i1->first, azimuth) < angle_distance(i2->first, azimuth))
        return i1->second;
    else
        return i2->second;
}

azimuthmap::Position NonuniformAzimuthMap::closest(double azimuth) const
{
    if (by_angle.empty()) throw std::runtime_error("closest() called on empty NonuniformAzimuthMap");

    auto i = by_angle.lower_bound(azimuth);

    // Result between the end and the beginning: assume it falls between
    // first and last
    if (i == by_angle.end() || i == by_angle.begin())
        return closest_of_two(azimuth, by_angle.rbegin(), by_angle.begin());

    // Exact match: return the Position
    if (i->first == azimuth)
        return i->second;

    // Return the closest between the previous element and this one
    std::map<double, Position>::const_iterator prev = i;
    --prev;
    return closest_of_two(azimuth, prev, i);
}

std::vector<azimuthmap::Position> NonuniformAzimuthMap::intersecting(double azimuth, double amplitude) const
{
    // Approximate our amplitude assuming the angles we have are close to
    // evenly spaced
    double my_amplitude = 360.0 / (double)by_angle.size();
    // Angles closer than this amount are considered the same for overlap detection
    static const double precision = 0.000000001;

    double lowest_azimuth = azimuth - amplitude / 2 - my_amplitude + precision;
    while (lowest_azimuth < 0) lowest_azimuth += 360;
    lowest_azimuth = fmod(lowest_azimuth, 360);
    double highest_azimuth = azimuth + amplitude / 2 + my_amplitude - precision;
    while (highest_azimuth < 0) highest_azimuth += 360;
    highest_azimuth = fmod(highest_azimuth, 360);

    std::vector<Position> res;

    if (lowest_azimuth <= highest_azimuth)
    {
        auto begin = by_angle.lower_bound(lowest_azimuth);
        auto end = by_angle.lower_bound(highest_azimuth);
        for (auto i = begin; i != end; ++i)
            res.push_back(i->second);
    } else {
        auto begin = by_angle.upper_bound(lowest_azimuth);
        auto end = by_angle.lower_bound(highest_azimuth);
        for (auto i = begin; i != by_angle.end(); ++i)
            res.push_back(i->second);
        for (auto i = by_angle.begin(); i != end; ++i)
            res.push_back(i->second);
    }

    return res;
}

template<typename T>
void AzimuthMap::polarscan_resample(const PolarScan<T>& src, PolarScan<T>& dst,
        std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger)
{
    if (src.beam_count < dst.beam_count)
        throw std::runtime_error("polarscan_resample currently only work for resampling to smaller volumes");

    // Do the merge
    for (unsigned i = 0; i < dst.beam_count; ++i)
        merger(src, *this, dst, i);
}

template<typename T>
void AzimuthMap::merger_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx)
{
    double dst_azimuth = (double)dst_idx / (double)dst.beam_count * 360.0;
    Position pos = azmap.closest(dst_azimuth);

    /// Copy the real elevation value
    dst.elevations_real(dst_idx) = src.elevations_real(pos.index);

    /// Copy the closest beam
    dst.row(dst_idx) = src.row(pos.index);
}

template<typename T>
void AzimuthMap::merger_max_of_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx)
{
    double dst_azimuth = (double)dst_idx / (double)dst.beam_count * 360.0;
    vector<Position> positions = azmap.intersecting(dst_azimuth, 360.0 / dst.beam_count);
    if (positions.empty()) throw std::runtime_error("no source beams found");

    double el_sum = 0;

    // Copy the first beam
    auto i = positions.cbegin();
    dst.row(dst_idx) = src.row(i->index);
    el_sum += src.elevations_real(i->index);

    // Take the maximum of all beam values
    for (++i; i != positions.end(); ++i)
    {
        for (unsigned bi = 0; bi < dst.beam_size; ++bi)
            if (dst(dst_idx, bi) < src(i->index, bi))
                dst(dst_idx, bi) = src(i->index, bi);
        el_sum += src.elevations_real(i->index);
    }

    // The real elevation of this beam is the average of the beams we used
    dst.elevations_real(dst_idx) = el_sum / positions.size();
}

template void AzimuthMap::polarscan_resample<double>(const PolarScan<double>& src, PolarScan<double>& dst,
        std::function<void(const PolarScan<double>&, const AzimuthMap&, PolarScan<double>&, unsigned)>);
template void AzimuthMap::polarscan_resample<unsigned char>(const PolarScan<unsigned char>& src, PolarScan<unsigned char>& dst,
        std::function<void(const PolarScan<unsigned char>&, const AzimuthMap&, PolarScan<unsigned char>&, unsigned)>);

template void AzimuthMap::merger_closest<double>(const PolarScan<double>&, const AzimuthMap&, PolarScan<double>&, unsigned);
template void AzimuthMap::merger_closest<unsigned char>(const PolarScan<unsigned char>&, const AzimuthMap&, PolarScan<unsigned char>&, unsigned);

template void AzimuthMap::merger_max_of_closest<double>(const PolarScan<double>&, const AzimuthMap&, PolarScan<double>&, unsigned);
template void AzimuthMap::merger_max_of_closest<unsigned char>(const PolarScan<unsigned char>&, const AzimuthMap&, PolarScan<unsigned char>&, unsigned);

}
