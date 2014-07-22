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

    Cleaner(const Variable<double>& var_Z, const Variable<double>& var_W, const Variable<double>& var_V)
        : Z_missing(var_Z.undetect), W_threshold(var_W.undetect), V_missing(var_V.nodata),
          bin_wind_magic_number(var_V.undetect )
    {
  }

    void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V);

//  void print_config(FILE* out) const;

//  void clean_beams(Beams<TB>& b, unsigned beam_size, std::vector<bool>& corrected) const;

#if 0
	BeamCleaner(TB V, TB Z, TB W)
	{
	  bin_wind_magic_number =  V;
	  Z_missing             =  Z;
	  W_threshold           =  W;
	  V_missing		=  V;
	}
	BeamCleaner(TB V, TB Z, TB W, TB VM)
	{
	  bin_wind_magic_number =  V;
	  Z_missing             =  Z;
	  W_threshold           =  W;
	  V_missing		= VM;
//std::cout<<bin_wind_magic_number<<" "<<Z_missing<<" "<<W_threshold<<" "<<V_missing<<std::endl;
	}
#endif

#if 0
	void print_config(FILE* out) const
	{
	    fprintf(out, "bin_wind_magic_number: %u, min_segment_length: %u, max_segment_length: %u\n",
	      (unsigned)bin_wind_magic_number, min_segment_length, max_segment_length);
	}
#endif

    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v) const;
};


}
}

#endif
