#include "azimuth_resample.h"
#include <cmath>

using namespace std;

namespace radarelab {
namespace algo {
namespace azimuthresample {

namespace {

inline double angle_distance(double first, double second)
{
    return 180.0 - std::fabs(std::fmod(std::fabs(first - second), 360.0) - 180.0);
}

inline std::pair<double, unsigned> closest_of_two(double azimuth,
        const std::pair<double, unsigned>& i1,
        const std::pair<double, unsigned>& i2)
{
    if (angle_distance(i1.first, azimuth) < angle_distance(i2.first, azimuth))
        return i1;
    else
        return i2;
}

}


AzimuthIndex::AzimuthIndex(const Eigen::VectorXd& azimuths)
{
    for (unsigned i = 0; i < azimuths.size(); ++i)
    {
        by_angle.insert(make_pair(azimuths(i), i));
        /*
        map<double, unsigned>::iterator old;
        bool inserted;
        tie(old, inserted) = by_angle.insert(make_pair(azimuths(i), i));
        if (!inserted)
        {
            fprintf(stderr, "IDX new %u old %u, val %f\n", i, old->second, old->first);
            throw std::runtime_error("source PolarScan has two beams with the same azimuth");
        }
        */
    }
}

pair<double, unsigned> AzimuthIndex::closest(double azimuth) const
{
    auto i = by_angle.lower_bound(azimuth);

    // Result between the end and the beginning: assume it falls between
    // first and last
    if (i == by_angle.end() || i == by_angle.begin())
        return closest_of_two(azimuth, *by_angle.rbegin(), *by_angle.begin());

    // Exact match: return the Position
    if (i->first == azimuth)
        return *i;

    // Return the closest between the previous element and this one
    std::map<double, unsigned>::const_iterator prev = i;
    --prev;
    return closest_of_two(azimuth, *prev, *i);
}

std::vector<pair<double, unsigned>> AzimuthIndex::intersecting(double dst_azimuth, double dst_amplitude, double src_amplitude) const
{
    // Compute the amplitude between our beams assuming the angles we have are
    // close to evenly spaced
    double my_semi_amplitude = src_amplitude / 2.0;
    // Angles closer than this amount are considered the same for overlap detection
    static const double precision = 0.000000001;

    double lowest_azimuth = dst_azimuth - dst_amplitude / 2 - my_semi_amplitude + precision;
    while (lowest_azimuth < 0) lowest_azimuth += 360;
    lowest_azimuth = fmod(lowest_azimuth, 360);
    double highest_azimuth = dst_azimuth + dst_amplitude / 2 + my_semi_amplitude - precision;
    while (highest_azimuth < 0) highest_azimuth += 360;
    highest_azimuth = fmod(highest_azimuth, 360);

    std::vector<pair<double, unsigned>> res;

    if (lowest_azimuth <= highest_azimuth)
    {
        auto begin = by_angle.lower_bound(lowest_azimuth);
        auto end = by_angle.lower_bound(highest_azimuth);
        for (auto i = begin; i != end; ++i)
            res.push_back(*i);
    } else {
        auto begin = by_angle.upper_bound(lowest_azimuth);
        auto end = by_angle.lower_bound(highest_azimuth);
        for (auto i = begin; i != by_angle.end(); ++i)
            res.push_back(*i);
        for (auto i = by_angle.begin(); i != end; ++i)
            res.push_back(*i);
    }

    return res;
}


template<typename T>
void Closest<T>::resample_polarscan(const PolarScan<T>& src, PolarScan<T>& dst, double src_beam_width) const
{
    AzimuthIndex index(src.azimuths_real);

    double max_distance = src_beam_width/2 + 360.0 / dst.beam_count / 2;

    // Do the merge
    for (unsigned dst_idx = 0; dst_idx < dst.beam_count; ++dst_idx)
    {
        double dst_azimuth = (double)dst_idx / (double)dst.beam_count * 360.0;

        pair<double, unsigned> pos = index.closest(dst_azimuth);

        if (angle_distance(dst_azimuth, pos.first) > max_distance)
        {
            /// Use the nominal elevation as the real one
            dst.elevations_real(dst_idx) = src.elevation;

            /// Use the nominal azimuth as the real one
            dst.azimuths_real(dst_idx) = dst_azimuth;
        } else {
            /// Copy the real elevation value
            dst.elevations_real(dst_idx) = src.elevations_real(pos.second);

            /// Copy the real azimuth value
            dst.azimuths_real(dst_idx) = pos.first;

            /// Copy the closest beam
            dst.row(dst_idx) = src.row(pos.second);
        }
    }
}


template<typename T>
void MaxOfClosest<T>::resample_polarscan(const PolarScan<T>& src, PolarScan<T>& dst, double src_beam_width) const
{
    AzimuthIndex index(src.azimuths_real);

    double max_distance = src_beam_width/2 + 360.0 / dst.beam_count / 2;

    // Do the merge
    for (unsigned dst_idx = 0; dst_idx < dst.beam_count; ++dst_idx)
    {
        double dst_azimuth = (double)dst_idx / (double)dst.beam_count * 360.0;
        vector<pair<double, unsigned>> positions = index.intersecting(dst_azimuth, 360.0 / dst.beam_count, src_beam_width);
        if (positions.empty())
        {
            /// Use the nominal elevation as the real one
            dst.elevations_real(dst_idx) = src.elevation;

            /// Use the nominal azimuth as the real one
            dst.azimuths_real(dst_idx) = dst_azimuth;
        } else {
            double el_sum = 0;
            double az_sum = 0;

            // Copy the first beam
            auto i = positions.cbegin();
            dst.row(dst_idx) = src.row(i->second);
            el_sum += src.elevations_real(i->second);
            az_sum += src.azimuths_real(i->second);

            // Take the maximum of all beam values
            for (++i; i != positions.end(); ++i)
            {
                for (unsigned bi = 0; bi < dst.beam_size; ++bi)
                    if (dst(dst_idx, bi) < src(i->second, bi))
                        dst(dst_idx, bi) = src(i->second, bi);
                el_sum += src.elevations_real(i->second);
                az_sum += src.azimuths_real(i->second);
            }

            // The real elevation of this beam is the average of the beams we used
            dst.elevations_real(dst_idx) = el_sum / positions.size();

            // The real azimuth of this beam is the average of the azimuths we used
            dst.azimuths_real(dst_idx) = az_sum / positions.size();
        }
    }
}

template class Closest<double>;
template class Closest<unsigned char>;
template class MaxOfClosest<double>;
template class MaxOfClosest<unsigned char>;

}
}
}

