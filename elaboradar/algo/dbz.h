#ifndef ELABORADAR_ALGO_DBZ_H
#define ELABORADAR_ALGO_DBZ_H

#include <elaboradar/volume.h>
#include <elaboradar/elev_fin.h>

namespace elaboradar {
namespace algo {

class DBZ
{
public:
    log4c_category_t* logging_category;

    double base_cell_size;
    //coeff a e b relazione Z-R
    double aMP, bMP;   /*  coeff a e b relazione Z-R  */

    DBZ();

    void setup(int month, double base_cell_size);

    /**
     *  ingresso:dbz in quel punto e attenuazione fin lì
     *  Doviak,Zrnic,1984 for rain as reported in cost 717 final document
     *
     *  @brief   funzione che calcola l'att enuazione totale 
     *  @details Ricevuto in ingresso il dato di Z in byte e l'attenuazione complessiva sul raggio fino al punto in considerazione, calcola l'attenuazione totale 
     *  @param[in]  DBZbyte valore in byte della Z nel pixel
     *  @param[in]  PIA attenuazione totale fino al punto
     *  @return att_tot attenuazione incrementata del contributo nel pixel corrente
     */
    double attenuation(unsigned char DBZbyte, double  PIA);

    // TODO: find better names for these:
    // RtoDBZ calcolato su aMP e bMP
    double RtoDBZ(double rain) const;
    double DBZtoR(double dbz) const;
    double DBZ_snow(double dbz) const;
    double DBZ_conv(double dbz) const;
    double RtoDBZ_class(double R) const;
    double DBZ_to_mp_func(double sample) const;
};

}
}

#endif

