#ifndef ELABORADAR_AZIMUTHMAP_H
#define ELABORADAR_AZIMUTHMAP_H

#include <map>
#include <vector>
#include <ostream>
#include <elaboradar/volume.h>

namespace elaboradar {
namespace azimuthmap {

struct Position
{
    double azimuth;
    unsigned index;
    Position(double azimuth, unsigned index)
        : azimuth(azimuth), index(index) {}
    bool operator==(const Position& pos) const
    {
        return azimuth == pos.azimuth && index == pos.index;
    }
    bool operator!=(const Position& pos) const
    {
        return azimuth != pos.azimuth || index != pos.index;
    }
};

std::ostream& operator<<(std::ostream& o, const Position& pos);

}

/**
 * Map azimut angles to PolarScan indices, and vice-versa
 */
struct AzimuthMap
{
    virtual ~AzimuthMap() {}

    /// Get the closest position to an azimuth angle
    virtual azimuthmap::Position closest(double azimuth) const = 0;

    /**
     * Get all the positions intersecting an angle centered on azimuth and with the given amplitude
     */
    virtual std::vector<azimuthmap::Position> intersecting(double azimuth, double amplitude) const = 0;

    template<typename T>
    void polarscan_resample(const PolarScan<T>& src, PolarScan<T>& dst,
            std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger);

    /**
     * Fill dst with data from src, coping with the two volumes having a different
     * number of beams per scan.
     *
     * Merger is the function used to merge beams from src into dst. It takes the
     * source PolarScan, the destination PolarScan and a vector with the indices of
     * the beams of src that need to be used.
     */
    template<typename T>
    void volume_resample(const Volume<T>& src, Volume<T>& dst,
            std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger)
    {
        // Copy volume metadata
        dst.quantity = src.quantity;
        dst.units = src.units;
        dst.load_info = src.load_info;

        for (unsigned iel = 0; iel < src.size(); ++iel)
        {
            const PolarScan<T>& src_scan = src.scan(iel);
            PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
            polarscan_resample(src_scan, dst_scan, merger);
        }
    }

    /**
     * Fill dst with data from src, coping with the two volumes having a different
     * number of beams per scan.
     *
     * Merger is the function used to merge beams from src into dst. It takes the
     * source PolarScan, the destination PolarScan and a vector with the indices of
     * the beams of src that need to be used.
     */
    template<typename T, typename AZMAP>
    static void volume_resample(const volume::Scans<T>& src, const std::vector<AZMAP>& azmaps, Volume<T>& dst,
            std::function<void(const PolarScan<T>&, const AzimuthMap&, PolarScan<T>&, unsigned)> merger)
    {
        // Copy volume metadata
        dst.load_info = src.load_info;
        dst.quantity = src.quantity;
        dst.units = src.units;

        for (unsigned iel = 0; iel < src.size(); ++iel)
        {
            const PolarScan<T>& src_scan = src.at(iel);
            PolarScan<T>& dst_scan = dst.make_scan(iel, src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
            azmaps[iel].polarscan_resample(src_scan, dst_scan, merger);
            dst_scan.nodata = src_scan.nodata;
            dst_scan.undetect = src_scan.undetect;
            dst_scan.gain = src_scan.gain;
            dst_scan.offset = src_scan.offset;
        }
    }

    template<typename T>
    static void merger_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx);

    template<typename T>
    static void merger_max_of_closest(const PolarScan<T>& src, const AzimuthMap& azmap, PolarScan<T>& dst, unsigned dst_idx);
};

/// Azimuth map that assumes equally spaced beams
struct UniformAzimuthMap : public AzimuthMap
{
    unsigned beam_count;

    template<typename T>
    UniformAzimuthMap(const PolarScan<T>& scan) : beam_count(scan.beam_count) {}
    UniformAzimuthMap(unsigned beam_count) : beam_count(beam_count) {}

    azimuthmap::Position closest(double azimuth) const override;
    std::vector<azimuthmap::Position> intersecting(double azimuth, double amplitude) const override;
};

/// Azimuth map that assumes that each beam has its own arbitrary angle
class NonuniformAzimuthMap : public AzimuthMap
{
protected:
    std::map<double, azimuthmap::Position> by_angle;

    template<typename ITER1, typename ITER2>
    azimuthmap::Position closest_of_two(double azimuth, const ITER1& i1, const ITER2& i2) const;

public:
    void add(double azimuth, unsigned idx)
    {
        by_angle.insert(std::make_pair(azimuth, azimuthmap::Position(azimuth, idx)));
    }

    azimuthmap::Position closest(double azimuth) const override;
    std::vector<azimuthmap::Position> intersecting(double azimuth, double amplitude) const override;
};

}

#endif
