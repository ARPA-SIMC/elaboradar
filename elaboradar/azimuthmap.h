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
