#ifndef ARCHIVIATORE_VOLUME_CLEANER_H
#define ARCHIVIATORE_VOLUME_CLEANER_H

#include <volume.h>

namespace cumbac {
namespace volume {

struct Cleaner
{
    const unsigned min_segment_length =  4;
    const unsigned max_segment_length = 40;

    const double Z_missing;
    const double W_threshold;
    const double V_missing;
    const double bin_wind_magic_number;

    Cleaner(double Z_missing, double W_threshold, double V_missing, double bin_wind_magic_number)
        : Z_missing(Z_missing), W_threshold(W_threshold), V_missing(V_missing), bin_wind_magic_number(bin_wind_magic_number) {}

    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v) const;

    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V);
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, double bin_wind_magic_number);
};


//  void print_config(FILE* out) const;

//  void clean_beams(Beams<TB>& b, unsigned beam_size, std::vector<bool>& corrected) const;

#if 0
	void print_config(FILE* out) const
	{
	    fprintf(out, "bin_wind_magic_number: %u, min_segment_length: %u, max_segment_length: %u\n",
	      (unsigned)bin_wind_magic_number, min_segment_length, max_segment_length);
	}
#endif

}
}

#endif
