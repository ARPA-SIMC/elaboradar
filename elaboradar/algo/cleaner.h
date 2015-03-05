#ifndef ELABORADAR_ALGO_CLEANER_H
#define ELABORADAR_ALGO_CLEANER_H

#include <elaboradar/volume.h>

namespace elaboradar {
namespace algo {

struct Cleaner
{
    const unsigned min_segment_length =  4;
    const unsigned max_segment_length = 40;

    const double Z_missing;
    const double W_threshold;
    const double V_missing;
    const double bin_wind_magic_number;
	const double sd_threshold = 2;

    Cleaner(double Z_missing, double W_threshold, double V_missing, double bin_wind_magic_number)
        : Z_missing(Z_missing), W_threshold(W_threshold), V_missing(V_missing), bin_wind_magic_number(bin_wind_magic_number) {}

    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v) const;
    std::vector<bool> Cleaner::clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& SD,int iray) const;

    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V);
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, double bin_wind_magic_number);
};

}
}

#endif
