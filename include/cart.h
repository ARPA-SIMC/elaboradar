#ifndef ELABORADAR_CART_H
#define ELABORADAR_CART_H

#include <elaboradar/matrix.h>
#include <elaboradar/volume.h>
#include <limits>

namespace elaboradar {

/**
 * Mapping of cartesian coordinates to raw azimuth angles and range distances.
 */
struct CoordinateMapping
{
    /// Beam size of the volume that we are mapping to cartesian coordinates
    const unsigned beam_size;
    /// Azimuth indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<double> map_azimuth;
    /// Range indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<double> map_range;

    /**
     * Build a cartography mapping cartesian coordinates to volume polar
     * indices.
     *
     * The mapping is a 1 to 1 mapping, without scaling.
     */
    CoordinateMapping(unsigned beam_size);

    /// Generate all the (azimuth, range) indices corresponding to a map point
    void sample(unsigned beam_count, unsigned x, unsigned y, std::function<void(unsigned, unsigned)>& f);
};


/**
 * Mapping of cartesian coordinates to specific azimuth and range volume indices
 */
class IndexMapping
{
public:
    /// Missing value in the azimuth and range index mappings
    static const unsigned missing = 0xffffffff;

    const unsigned height;
    const unsigned width;

    /// Azimuth indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned> map_azimuth;
    /// Range indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned> map_range;

    IndexMapping(unsigned height, unsigned width);

    /// Copy data from the polar scan src to the cartesian map dst
    template<typename SRC, typename DST>
    void to_cart(const PolarScan<SRC>& src, Matrix2D<DST>& dst)
    {
        // In case dst is not a square with side beam_size*2, center it
        int dx = ((int)width - dst.cols()) / 2;
        int dy = ((int)height - dst.rows()) / 2;

        for (unsigned y = 0; y < dst.rows(); ++y)
        {
            if (y + dy < 0 || y + dy >= height) continue;

            for (unsigned x = 0; x < dst.cols(); ++x)
            {
                if (x + dx < 0 || x + dx >= width) continue;

                auto azimuth = map_azimuth(y + dy, x + dx);
                auto range = map_range(y + dy, x + dx);

                if (azimuth == missing || range == missing) continue;
                if (azimuth >= src.beam_count || range >= src.beam_size) continue;

                dst(y, x) = src(azimuth, range);
            }
        }
    }

    /// Fill the cartesian map dst with the output of the function src(azimuth, range)
    template<typename T>
    void to_cart(const std::function<T(unsigned, unsigned)>& src, Matrix2D<T>& dst)
    {
        // In case dst is not a square with side beam_size*2, center it
        int dx = ((int)width - dst.cols()) / 2;
        int dy = ((int)height - dst.rows()) / 2;

        for (unsigned y = 0; y < dst.rows(); ++y)
        {
            if (y + dy < 0 || y + dy >= height) continue;

            for (unsigned x = 0; x < dst.cols(); ++x)
            {
                if (x + dx < 0 || x + dx >= width) continue;

                auto azimuth = map_azimuth(y + dy, x + dx);
                auto range = map_range(y + dy, x + dx);

                if (azimuth == missing || range == missing) continue;
                dst(y, x) = src(azimuth, range);
            }
        }
    }
};


/**
 * Index mapping where the pixel size corresponds to the radar cell size
 */
struct FullsizeIndexMapping : public IndexMapping
{
    FullsizeIndexMapping(unsigned beam_size);

    /**
     * Map cartesian cardinates to polar volume indices. When a cartesian
     * coordinate maps to more than one polar value, take the one with the
     * maximum data value.
     */
    void map_max_sample(const PolarScan<double>& scan);

    /**
     * Same as map_max_sample(PolarScan), but reuse an existing
     * CoordinateMapping.
     */
    void map_max_sample(const PolarScan<double>& scan, const CoordinateMapping& mapping);
};


/**
 * Index mapping with arbitrary pixel size.
 *
 * The scaled image will still be a square image, and its side must be a
 * divisor of the full image size (beam_size * 2).
 */
struct ScaledIndexMapping : public IndexMapping
{
    unsigned fullsize_pixels_per_scaled_pixel;
    int image_offset;

    ScaledIndexMapping(unsigned beam_size, unsigned image_side, unsigned fullsize_pixels_per_scaled_pixel);

    /**
     * Map cartesian cardinates to polar volume indices. When a cartesian
     * coordinate maps to more than one polar value, take the one with the
     * maximum data value.
     *
     * It requires a FullsizeIndexMapping that has been initialized with
     * map_max_sample().
     */
    void map_max_sample(const PolarScan<double>& scan, const FullsizeIndexMapping& mapping);
};

}

#endif
