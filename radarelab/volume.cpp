#include "volume.h"
#include "logging.h"
#include <stdexcept>
#include <algorithm>

static const double ker = 8494666.666667; // c'è qualcosa in geo_par.h

using namespace std;

namespace radarelab {

PolarScanBase::PolarScanBase(unsigned beam_count, unsigned beam_size)
    : beam_count(beam_count), beam_size(beam_size),
      azimuths_real(beam_count), elevations_real(beam_count)
{
}

PolarScanBase::PolarScanBase(const PolarScanBase& s)
    : beam_count(s.beam_count), beam_size(s.beam_size),
      azimuths_real(s.azimuths_real),
      elevation(s.elevation),
      elevations_real(s.elevations_real),
      cell_size(s.cell_size)
{
}

double PolarScanBase::height(unsigned rg, double beam_half_width)
{
    return sample_height(elevation + beam_half_width, (double)rg*cell_size) / 1000.; // km
}

double PolarScanBase::diff_height(unsigned rg_start, unsigned rg_end)
{
    return fabs(height(rg_end) - height(rg_start));
}

double PolarScanBase::sample_height(unsigned cell_idx) const
{
    return sample_height(elevation, (double)cell_idx * cell_size);
}

double PolarScanBase::sample_height(double elevation, double range, double equiv_earth_radius)
{
    return sqrt(
            range * range
            + equiv_earth_radius * equiv_earth_radius
            + 2. * equiv_earth_radius * range * sin(elevation * M_PI / 180.)
           ) - equiv_earth_radius; //meters
}

double PolarScanBase::sample_height(double elevation, double range)
{
    return sample_height(elevation, range, ker);
}

void VolumeStats::print(FILE* out)
{
    fprintf(out, "Nel    Zeros     Ones   Others      Sum\n");
    for (size_t iel =0; iel<count_zeros.size(); ++iel){
        fprintf(out, "%4zu %8u %8u %8u %8u\n",iel,count_zeros[iel],count_ones[iel],count_others[iel],sum_others[iel]);
    }
}

template class PolarScan<double>;
template class Volume<double>;

namespace volume {
template class Scans<double>;
}

}
