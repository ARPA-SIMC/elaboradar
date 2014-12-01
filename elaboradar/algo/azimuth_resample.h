#ifndef ELABORADAR_ALGO_AZIMUTHRESAMPLE_H
#define ELABORADAR_ALGO_AZIMUTHRESAMPLE_H

#include <map>
#include <vector>
#include <ostream>
#include <elaboradar/volume.h>

namespace elaboradar {
namespace algo {
namespace azimuthresample {

/// Index beam positions by azimuth angles
class AzimuthIndex
{
protected:
    std::map<double, unsigned> by_angle;

public:
    /// Build an index with the azimuths of a PolarScan
    AzimuthIndex(const Eigen::VectorXd& azimuths);

    /// Get the closest position to an azimuth angle
    std::pair<double, unsigned> closest(double azimuth) const;

    /**
     * Get all the positions intersecting an angle centered on azimuth and with the given amplitude
     */
    std::vector<std::pair<double, unsigned>> intersecting(double azimuth, double amplitude) const;
};

/// Resample a volume one level at a time
template<typename T>
struct LevelwiseResampler
{
    /**
     * Fill dst with data from src, using the given merger function
     */
    virtual void resample_polarscan(const PolarScan<T>& src, PolarScan<T>& dst) const = 0;


    /// Merge 
    //virtual void merger(const PolarScan<T>&, PolarScan<T>&, unsigned) = 0;

    /**
     * Fill dst with data from src, coping with the two volumes having a different
     * number of beams per scan.
     *
     * Merger is the function used to merge beams from src into dst. It takes the
     * source PolarScan, the destination PolarScan and a vector with the indices of
     * the beams of src that need to be used.
     */
    void resample_volume(const Volume<T>& src, Volume<T>& dst) const
    {
        // Copy volume metadata
        dst.quantity = src.quantity;
        dst.units = src.units;
        dst.load_info = src.load_info;

        for (unsigned iel = 0; iel < src.size(); ++iel)
        {
            const PolarScan<T>& src_scan = src.scan(iel);
            PolarScan<T>& dst_scan = dst.append_scan(src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
            resample_polarscan(src_scan, dst_scan);
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
    void resample_volume(const volume::Scans<T>& src, Volume<T>& dst) const
    {
        // Copy volume metadata
        dst.load_info = src.load_info;
        dst.quantity = src.quantity;
        dst.units = src.units;

        for (unsigned iel = 0; iel < src.size(); ++iel)
        {
            const PolarScan<T>& src_scan = src.at(iel);
            PolarScan<T>& dst_scan = dst.append_scan(src_scan.beam_size, src_scan.elevation, src_scan.cell_size);
            resample_polarscan(src_scan, dst_scan);
            dst_scan.nodata = src_scan.nodata;
            dst_scan.undetect = src_scan.undetect;
            dst_scan.gain = src_scan.gain;
            dst_scan.offset = src_scan.offset;
        }
    }
};

template<typename T>
struct Closest : public LevelwiseResampler<T>
{
    virtual void resample_polarscan(const PolarScan<T>& src, PolarScan<T>& dst) const override;
};

template<typename T>
struct MaxOfClosest : public LevelwiseResampler<T>
{
    virtual void resample_polarscan(const PolarScan<T>& src, PolarScan<T>& dst) const override;
};

}
}
}

#endif

