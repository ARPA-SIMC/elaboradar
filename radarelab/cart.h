/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_CART_H
#define RADARELAB_CART_H

#include <radarelab/matrix.h>
#include <radarelab/volume.h>
#include <limits>
#include <vector>

namespace radarelab {

/**
 * Mapping of cartesian coordinates to raw azimuth angles and range distances.
 *
 * Cartesian coordinates follow the pixels of a rendered image, with (0,0) in
 * the top left corner, X increasing west to east, Y increasing north to south.
 */
struct CoordinateMapping
{
    /// Beam size of the volume that we are mapping to cartesian coordinates
    const unsigned beam_size;

    /**
     * Azimuth indices to use to lookup a map point in a volume.
     *
     * -1 means no mapping.
     *
     * Each value in the matrix is the azimuth pointing to the center of the
     * corresponding cartesian pixel.
     */
    Matrix2D<double> map_azimuth;

    /**
     * Range indices to use to lookup a map point in a volume.
     *
     * -1 means no mapping.
     *
     * Each value in the matrix is the distance (in cells) from the radar to
     * the center of the corresponding cartesian pixel.
     */
    Matrix2D<double> map_range;

    /**
     * Build a cartography mapping cartesian coordinates to volume polar
     * indices.
     *
     * The mapping is a 1 to 1 mapping, without scaling.
     */
    CoordinateMapping(unsigned beam_size);

    /// Generate all the (azimuth, range) indices corresponding to a map point
    void sample(unsigned beam_count, unsigned x, unsigned y, std::function<void(unsigned, unsigned)>& f) const;
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
    void to_cart(const PolarScan<SRC>& src, Matrix2D<DST>& dst) const
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
    void to_cart(const std::function<T(unsigned, unsigned)>& src, Matrix2D<T>& dst) const
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
    const CoordinateMapping& mapping;
    unsigned fullsize_pixels_per_scaled_pixel;
    /// Image offset in full size pixels
    int image_offset;

    ScaledIndexMapping(const CoordinateMapping& mapping, unsigned image_side, unsigned fullsize_pixels_per_scaled_pixel);

    /// Generate all the (azimuth, range) indices corresponding to a map point
    void sample(unsigned beam_count, unsigned x, unsigned y, std::function<void(unsigned, unsigned)>& f);

    /**
     * Map cartesian cardinates to polar volume indices. When a cartesian
     * coordinate maps to more than one polar value, take the one with the
     * maximum data value.
     *
     * It requires a FullsizeIndexMapping that has been initialized with
     * map_max_sample().
     */
    void map_max_sample(const PolarScan<double>& scan, const FullsizeIndexMapping& mapping);

    /// Fill the cartesian map dst with the output of the function src(azimuth, range)
    template<typename T>
    void to_cart_average(const PolarScan<double>& src, std::function<T(const std::vector<double>&)>& convert, Matrix2D<T>& dst) const
    {
        // In case dst is not a square with side beam_size*2, center it
        int dx = ((int)width - dst.cols()) / 2;
        int dy = ((int)height - dst.rows()) / 2;
        std::vector<double> samples;

        for (unsigned y = 0; y < dst.rows(); ++y)
        {
            if (y + dy < 0 || y + dy >= height) continue;

            for (unsigned x = 0; x < dst.cols(); ++x)
            {
                if (x + dx < 0 || x + dx >= width) continue;

                samples.clear();
                std::function<void(unsigned, unsigned)> compute_average = [&samples, &src](unsigned azimuth, unsigned range) {
                    if (azimuth < 0 || azimuth > src.beam_count) return;
                    if (range < 0 || range > src.beam_size) return;
                    samples.push_back(src(azimuth, range));
                };

                for(unsigned sy = 0; sy < fullsize_pixels_per_scaled_pixel; ++sy)
                    for(unsigned sx = 0; sx < fullsize_pixels_per_scaled_pixel; ++sx)
                    {
                        int src_x = x * fullsize_pixels_per_scaled_pixel + sx + image_offset;
                        int src_y = y * fullsize_pixels_per_scaled_pixel + sy + image_offset;
                        if (src_x < 0 || src_x >= mapping.beam_size * 2 || src_y < 0 || src_y >= mapping.beam_size * 2)
                            continue;
                        mapping.sample(src.beam_count, src_x, src_y, compute_average);
                    }

                if (!samples.empty())
                    dst(y + dy, x + dx) = convert(samples);
            }
        }
    }
};

}

#endif
