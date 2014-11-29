#ifndef ELABORADAR_ALGO_TOP_H
#define ELABORADAR_ALGO_TOP_H

#include <elaboradar/volume.h>

namespace elaboradar {
namespace algo {

template<typename T>
void compute_top(const Volume<T>& volume, T threshold, Matrix2D<unsigned char>& top)
{
    top.resize(400, volume.max_beam_size());
    top.fill(0);
    for (unsigned l=0; l < volume.size(); ++l)
    {
        const auto& scan = volume[l];
        for (int i=0; i < NUM_AZ_X_PPI; ++i)
        {
            const double elevaz = scan.elevations_real(i); //--- elev reale in gradi
            for (unsigned k = 0; k < scan.beam_size; ++k)
                if (scan.get(i, k) > threshold)
                    //top in ettometri
                    top(i, k) = (unsigned char)(PolarScanBase::sample_height(
                                elevaz, (k + 0.5) * scan.cell_size) / 100.);
        }
    }
}

}
}

#endif

