#ifndef RADARELAB_ALGO_VPR_H
#define RADARELAB_ALGO_VPR_H

#include <radarelab/volume.h>
#include <radarelab/algo/dbz.h>
#include <vector>

namespace radarelab {
namespace algo {

struct InstantaneousVPR
{
    const Volume<double>& volume;
    const Volume<unsigned char>& qual;  ///< quality flags for the volume
    /// 1 for points where VPR can be computed, else 0
    Volume<unsigned char>& flag_vpr;
    int az_min;
    int az_max;
    DBZ dbz;


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
    int compute(long int *cv, long int *ct, std::vector<float>& vpr1, std::vector<long int>& area_vpr);
};

}
}

#endif
