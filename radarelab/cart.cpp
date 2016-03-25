#include <radarelab/cart.h>
#include <cmath>

namespace radarelab {

CoordinateMapping::CoordinateMapping(unsigned beam_size)
    : beam_size(beam_size),
      map_azimuth(beam_size * 2, beam_size * 2),
      map_range(beam_size * 2, beam_size * 2)
{
    for (unsigned y = 0; y < beam_size * 2; ++y)
        for (unsigned x = 0; x < beam_size * 2; ++x)
        {
            // x and y centered on the map center
            double absx = (double)x - beam_size;
            double absy = (double)y - beam_size;
            //absx += (absx < 0) ? -0.5 : 0.5;
            //absy += (absy < 0) ? -0.5 : 0.5;
            absx += 0.5;
            absy += 0.5;

            // Compute range
            map_range(y, x) = hypot(absx, absy);

            // Compute azimuth
            double az;
            if (absx > 0)
                if (absy < 0)
                    az = atan(absx / -absy) * M_1_PI * 180.;
                else
                    az = 90.0 + atan(absy / absx) * M_1_PI * 180.;
            else
                if (absy > 0)
                    az = 180 + atan(-absx / absy) * M_1_PI * 180.;
                else
                    az = 270 + atan(-absy / -absx) * M_1_PI * 180.;
            map_azimuth(y, x) = az;
        }
}

void CoordinateMapping::sample(unsigned beam_count, unsigned x, unsigned y, std::function<void(unsigned, unsigned)>& f) const
{
    // Map cartesian coordinates to angles and distances
    unsigned range_idx = floor(map_range(y, x));
    if (range_idx >= beam_size) return;

    // Exact angle
    double az = map_azimuth(y, x);

#if 0
    // Iterate indices 0.45° before and after
    int az_min = floor((az - .45) * beam_count / 360.0);
    int az_max = ceil((az + .45) * beam_count / 360.0);
    if (az_min < 0)
    {
        az_min += beam_count;
        az_max += beam_count;
    }
#endif

    // Iterate on angles that actually overlap with the map cell
    double d_az = M_1_PI * 180. / (range_idx + 0.5) / 2;
    int az_min = round((az - d_az) * beam_count / 360.);
    int az_max = round((az + d_az) * beam_count / 360.);

    // Iterate all points between az_min and az_max
    for (unsigned iaz = az_min; iaz <= (unsigned)az_max; ++iaz)
        f(iaz % beam_count, range_idx);
}

const unsigned IndexMapping::missing;

IndexMapping::IndexMapping(unsigned height, unsigned width)
    : height(height), width(width),
      map_azimuth(Matrix2D<unsigned>::Constant(height, width, missing)),
      map_range(Matrix2D<unsigned>::Constant(height, width, missing))
{
}

FullsizeIndexMapping::FullsizeIndexMapping(unsigned beam_size)
    : IndexMapping(beam_size * 2, beam_size * 2)
{
}

void FullsizeIndexMapping::map_max_sample(const PolarScan<double>& scan)
{
    CoordinateMapping raw_mapping(scan.beam_size);
    map_max_sample(scan, raw_mapping);
}

void FullsizeIndexMapping::map_max_sample(const PolarScan<double>& scan, const CoordinateMapping& mapping)
{
    for (unsigned y = 0; y < height; ++y)
        for (unsigned x = 0; x < width; ++x)
        {
            bool first = true;
            double maxval;
            unsigned maxval_azimuth;
            unsigned maxval_range;
            std::function<void(unsigned, unsigned)> f = [&](unsigned azimuth, unsigned range) {
                double sample = scan.get(azimuth, range);
                if (first || sample > maxval)
                {
                    maxval = sample;
                    maxval_azimuth = azimuth;
                    maxval_range = range;
                    first = false;
                }
            };
            mapping.sample(scan.beam_count, x, y, f);

            // If nothing has been found, skip this sample
            if (first) continue;

            // Store the indices for this mapping
            map_azimuth(y, x) = maxval_azimuth;
            map_range(y, x) = maxval_range;
        }
}


ScaledIndexMapping::ScaledIndexMapping(const CoordinateMapping& mapping, unsigned image_side, unsigned fullsize_pixels_per_scaled_pixel)
    : IndexMapping(image_side, image_side), mapping(mapping),
      fullsize_pixels_per_scaled_pixel(fullsize_pixels_per_scaled_pixel),
      image_offset(((int)mapping.beam_size * 2 - (int)image_side * (int)fullsize_pixels_per_scaled_pixel) / 2)
{
    if ((image_side * fullsize_pixels_per_scaled_pixel) % 2 != 0)
        throw std::runtime_error("the image cannot be properly centered on the full size image");
}

void ScaledIndexMapping::sample(unsigned beam_count, unsigned x, unsigned y, std::function<void(unsigned, unsigned)>& f)
{
    // Load each sample with a value from a 4x4 window on the original image
    for(unsigned sy = 0; sy < fullsize_pixels_per_scaled_pixel; ++sy)
        for(unsigned sx = 0; sx < fullsize_pixels_per_scaled_pixel; ++sx)
        {
            // Use the full size mapping to get the volume value at this point
            int src_x = x * fullsize_pixels_per_scaled_pixel + sx + image_offset;
            int src_y = y * fullsize_pixels_per_scaled_pixel + sy + image_offset;
            if (src_x < 0 || src_x >= mapping.beam_size || src_y < 0 || src_y >= mapping.beam_size) continue;

            // Generate all the azimuth/range elements for this point in the
            // full size map
            std::function<void(unsigned, unsigned)> fullsize_sample = [&f](unsigned azimuth, unsigned range) {
                f(azimuth, range);
            };
            mapping.sample(beam_count, src_x, src_y, fullsize_sample);
        }
}

void ScaledIndexMapping::map_max_sample(const PolarScan<double>& scan, const FullsizeIndexMapping& mapping)
{
    // ciclo sui punti della nuova matrice. per il primo prenderò il massimo tra i primi sedici etc..
    for (unsigned y = 0; y < height; ++y)
        for (unsigned x = 0; x < width; ++x)
        {
            bool first = true;
            double maxval;
            unsigned maxval_az = missing;
            unsigned maxval_range = missing;

            // Load each sample with a value from a 4x4 window on the original image
            for(unsigned sy = 0; sy < fullsize_pixels_per_scaled_pixel; ++sy)
                for(unsigned sx = 0; sx < fullsize_pixels_per_scaled_pixel; ++sx)
                {
                    // Use the full size mapping to get the volume value at this point
                    int src_x = x * fullsize_pixels_per_scaled_pixel + sx + image_offset;
                    int src_y = y * fullsize_pixels_per_scaled_pixel + sy + image_offset;
                    if (src_x < 0 || src_x >= mapping.width || src_y < 0 || src_y >= mapping.height) continue;
                    unsigned range = mapping.map_range(src_y, src_x);
                    if (range == missing) continue;
                    unsigned az = mapping.map_azimuth(src_y, src_x);
                    double sample = scan.get(az, range);

                    // Find the source point with the maximum sample
                    if (first || sample > maxval)
                    {
                        maxval = sample;
                        maxval_az = az;
                        maxval_range = range;
                        first = false;
                    }
                }

            // If nothing has been found, skip this sample
            if (first) continue;

            // Store the indices for this mapping
            map_azimuth(y, x) = maxval_az;
            map_range(y, x) = maxval_range;
        }
}

}
