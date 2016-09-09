#ifndef RADARELAB_ALGO_VPR_H
#define RADARELAB_ALGO_VPR_H

#include <radarelab/volume.h>
#include <radarelab/algo/dbz.h>
#include <vector>
#include <array>

namespace radarelab {
namespace algo {

constexpr unsigned VPR_NMAXLAYER = 70;
constexpr float VPR_MISSING = -9999.;

struct VPR
{
    std::array<float, VPR_NMAXLAYER> val;
    std::array<long int, VPR_NMAXLAYER> area;

    VPR(const VPR&) = default;
    VPR(VPR&&) = default;
    VPR& operator=(const VPR&) = default;
    VPR& operator=(VPR&&) = default;
    VPR()
    {
        for (unsigned i = 0; i < size(); ++i)
        {
            val[i] = VPR_MISSING;
            area[i] = 0;
        }
    }

    size_t size() const { return VPR_NMAXLAYER; }
};


struct Livmin
{
    /// Index in VPR of the minimum level
    unsigned idx = 0;

    /// Value of the minimum level
    int livmin = 0;

    /// True if the minimum level has been found
    bool found = false;

    Livmin(const VPR& vpr);
};


struct InstantaneousVPR
{
    const Volume<double>& volume;
    const Volume<unsigned char>& qual;  ///< quality flags for the volume
    /// 1 for points where VPR can be computed, else 0
    Volume<unsigned char>& flag_vpr;
    int az_min;
    int az_max;
    DBZ dbz;
    long int cv = 0;
    long int ct = 0;
    VPR vpr;
    bool success = false;

#if 0
    /**
     *  crea vpr istantaneo
     *  @brief funzione che calcola il profilo istantaneo  secondo il metodo di Germann e Joss (2003)
     *  @details   calcola il VPR istantaneo secondo il metodo di Germann e Joss (2003) 
     *  Per il calcolo si considerano i punti con Z>THR_VPR, qualità>QMIN_VPR, BeamBlocking<20 percento e clutter free all'interno del volume scelto.
     *  Il profilo è poi soggetto a quality check e eventualmente viene rigettato (return(1)) 
     *  @param[out] cv volume precipitante
     *  @param[out] ct volume totale 
     *  @param[out] vpr1 vettore vpr istantaneo
     *  @param[out] area_vpr vettore area di ogni stato 
     *  @return  0 se ok 1 se fallisce
     */ 
    int func_vpr(long int *cv, long int *ct, std::vector<float>& vpr1, std::vector<long int>& area_vpr);
#endif

    InstantaneousVPR(const Volume<double>& volume, const Volume<unsigned char>& qual, Volume<unsigned char>& flag_vpr, int az_min, int az_max);
    void compute();
};

}
}

#endif
