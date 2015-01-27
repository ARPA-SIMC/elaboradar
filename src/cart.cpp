#include "cart.h"
#include <cmath>

namespace elaboradar {

CoordinateMapping::CoordinateMapping(unsigned beam_size)
    : beam_size(beam_size),
      map_azimuth(beam_size * 2, beam_size * 2),
      map_range(beam_size * 2, beam_size * 2)
{
    for (unsigned y = 0; y < beam_size * 2; ++y)
        for (unsigned x = 0; x < beam_size * 2; ++x)
        {
            // x and y centered on the map center
            double absx = x + 0.5 - beam_size;
            double absy = y + 0.5 - beam_size;

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


const unsigned IndexMapping::missing;

IndexMapping::IndexMapping(unsigned beam_size)
    : beam_size(beam_size),
      map_azimuth(Matrix2D<unsigned>::Constant(beam_size * 2, beam_size * 2, missing)),
      map_range(Matrix2D<unsigned>::Constant(beam_size * 2, beam_size * 2, missing))
{
}

void IndexMapping::map_max_sample(const PolarScan<double>& scan)
{
    CoordinateMapping raw_mapping(scan.beam_size);
    map_max_sample(scan, raw_mapping);
}

void IndexMapping::map_max_sample(const PolarScan<double>& scan, const CoordinateMapping& mapping)
{
    for (unsigned y = 0; y < beam_size * 2; ++y)
        for (unsigned x = 0; x < beam_size * 2; ++x)
        {
            // Map cartesian coordinates to angles and distances
            unsigned range = floor(mapping.map_range(y, x));
            if (range >= scan.beam_size) continue;

            // Exact angle
            double az = mapping.map_azimuth(y, x);

            // Iterate indices 0.45° before and after
            int az_min = floor((az - .45) * scan.beam_count / 360);
            int az_max = ceil((az + .45) * scan.beam_count / 360);
            if (az_min < 0)
            {
                az_min += scan.beam_count;
                az_max += scan.beam_count;
            }

            // Look for the position of the maximum value of the scan in
            // this range
            bool first = true;
            double maxval;
            unsigned maxval_idx;
            for (unsigned iaz = az_min; iaz < (unsigned)az_max; ++iaz)
            {
                unsigned azidx = iaz % scan.beam_count;
                unsigned char sample = scan.get(azidx, range);

                if (first || sample > maxval)
                {
                    maxval = sample;
                    maxval_idx = azidx;
                    first = false;
                }
            }

            // If nothing has been found, skip this sample
            if (first) continue;

            // map_azimuth(y, x) = round(az * scan.beam_count / 360);
            map_azimuth(y, x) = maxval_idx;
            map_range(y, x) = range;
        }
}

}
