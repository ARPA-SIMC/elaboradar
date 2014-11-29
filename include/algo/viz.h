#ifndef ARCHIVIATORE_VIZ_CLASS_H
#define ARCHIVIATORE_VIZ_CLASS_H

#include <elaboradar/logging.h>
#include <elaboradar/matrix.h>

namespace elaboradar {
struct CylindricalVolume;

namespace algo {

struct CalcoloVIZ
{
    log4c_category_t* logging_category;

    const CylindricalVolume& cil;

    const unsigned x_size;
    const unsigned z_size;
    const double htbb;
    const double hbbb;
    const double t_ground;

    Matrix2D<unsigned char> conv_VIZ;
    Matrix2D<unsigned char> stratiform;

    CalcoloVIZ(const CylindricalVolume& cil, double htbb, double hbbb, double t_ground);

    /**
     *  classifica tramite Vertical Integrated Z
     *  @brief funzione  che classifica secondo il metodo VIZ
     *  @details calcolo per ogni pixel polare l'integrale verticale esclusa la fascia della bright band
     */
    void classifico_VIZ();
};

}
}

#endif
