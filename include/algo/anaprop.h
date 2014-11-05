#ifndef ELABORADAR_ALGO_ANAPROP_H
#define ELABORADAR_ALGO_ANAPROP_H

#include <volume.h>
#include <volume/elev_fin.h>

namespace elaboradar {
namespace algo {

namespace anaprop {
// Matrici per statistiche
struct GridStats
{
    // dim azim griglia per stat anap
    const unsigned step_stat_az = 25;
    // dim range griglia per stat anap
    const unsigned step_stat_range = 40;

    // Number of cells in the azimut direction
    unsigned size_az = 0;
    // Number of cells in the beam direction
    unsigned size_beam = 0;

    // statistica anaprop
    unsigned* stat_anap = 0;
    // contatore punti dentro ogni box per statistica
    unsigned* stat_tot = 0;
    // statistica beam blocking
    unsigned* stat_bloc = 0;
    // statistica cambio elevazione rispetto mappa statica
    unsigned* stat_elev = 0;

    GridStats();
    ~GridStats();

    void init(const Volume<double>& volume);

    inline unsigned idx(unsigned az, unsigned beam) const
    {
        return az / step_stat_az * size_beam + beam / step_stat_range;
    }

    void incr_anap(unsigned az, unsigned beam) { stat_anap[idx(az, beam)]++; }
    void incr_tot(unsigned az, unsigned beam) { stat_tot[idx(az, beam)]++; }
    void incr_elev(unsigned az, unsigned beam) { stat_elev[idx(az, beam)]++; }
    void incr_bloc(unsigned az, unsigned beam, unsigned amount) { stat_bloc[idx(az, beam)]++; }

    unsigned count(unsigned az, unsigned beam) const
    {
        return stat_tot[idx(az, beam)];
    }

    unsigned char perc_anap(unsigned az, unsigned beam) const
    {
        return stat_anap[idx(az, beam)] * 100 / stat_tot[idx(az, beam)];
    }

    unsigned char perc_elev(unsigned az, unsigned beam) const
    {
        return stat_elev[idx(az, beam)] * 100 / stat_tot[idx(az, beam)];
    }

    unsigned char perc_bloc(unsigned az, unsigned beam) const
    {
        return stat_bloc[idx(az, beam)] * 100 / stat_tot[idx(az, beam)];
    }
};
}

template<class T>
class Anaprop
{
public:
    static const int ANAP_OK        = 0;
    static const int ANAP_YES       = 1;
    static const int ANAP_NODAT     = 2;
    static const int ANAP_NOCONTROL = 3;

    log4c_category_t* logging_category;

    // Anaprop configuration

    /**
     * Threshold value to use on the texture parameter, used to discriminate
     * meteorological echoes from non-meteorological echoes
     */
    double conf_texture_threshold = 3;
    bool do_quality = false;
    bool do_beamblocking = false;
    bool do_bloccorr = false;


    // Output data from the anaprop algorithm

    anaprop::GridStats grid_stats;
    volume::ElevFin<T> elev_fin;
    Matrix2D<unsigned char> dato_corrotto; // uscita controllo anaprop in coordinate azimut range
    Matrix2D<unsigned short> quota; /*quota fascio in prop standard e elev reali in coordinate azimut range*/

    Anaprop();

    // Initialise basic structures
    void init(const Volume<T>& volume);

    // Initialise basic structures and copy first_level_static to elev_fin
    void init_elev_fin_static(const Volume<T>& volume, const PolarScan<unsigned char>& first_level_static);

    // Initialise basic structures and performe anomalous propagation using sd and testing vertical gradient
    void remove(
            Volume<T>& volume,
            PolarScan<unsigned char>& beam_blocking,
            const PolarScan<unsigned char>& first_level,
            const PolarScan<unsigned char>& first_level_static,
            const Volume<double>& sd);

    // Initialise basic structures and performe anomalous propagation testing vertical gradient
    void remove_without_SD(
            Volume<T>& volume,
            PolarScan<unsigned char>& beam_blocking,
            const PolarScan<unsigned char>& first_level,
            const PolarScan<unsigned char>& first_level_static,
            const Volume<double>& sd);
};

}
}

#endif

