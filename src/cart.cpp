#include "cart.h"
#include <cmath>

namespace elaboradar {

const unsigned CartFullRes::missing;


CartFullRes::CartFullRes(const PolarScan<double>& scan, bool ignore_data)
    : beam_size(scan.beam_size),
      map_azimuth(Matrix2D<unsigned>::Constant(beam_size * 2, beam_size * 2, missing)),
      map_range(Matrix2D<unsigned>::Constant(beam_size * 2, beam_size * 2, missing))
{
    for (unsigned y = 0; y < beam_size * 2; ++y)
        for (unsigned x = 0; x < beam_size * 2; ++x)
        {
            // x and y centered on the map center
            double absx = x + 0.5 - beam_size;
            double absy = y + 0.5 - beam_size;
            unsigned short range = floor(hypot(absx, absy));

            // If we are out of range of the base scan, leave missing in the
            // mapping matrices
            if (!ignore_data && range >= scan.beam_size) continue;

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

            if (ignore_data)
            {
                map_azimuth(y, x) = (unsigned short)floor(az * scan.beam_count / 360) % scan.beam_count;
                map_range(y, x) = range;
                continue;
            }

            // Iterate pixels 0.45° before and after
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
