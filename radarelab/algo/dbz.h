/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_ALGO_DBZ_H
#define RADARELAB_ALGO_DBZ_H

#include <radarelab/volume.h>
#include <radarelab/elev_fin.h>

namespace radarelab {
namespace algo {
/**
 * Class to manage reflectivity functions (simply attenuation correction, conversion between Z, dBZ, R)
 */
class DBZ
{
public:
    log4c_category_t* logging_category;			///< logging category label

    double base_cell_size;				///< cella size dimension
    double aMP, bMP;   					////< Marshall-Palmer coefficient for Z-R relationship

/// Constriuctor
    DBZ();


/**
 * @brief Seasonal setup function
 * @param [in] month - month 
 * @param [in] base_cell_size - cell size dimension [m]
 */
    void setup(int month, double base_cell_size);

    /**
     *
     *  @brief   funzione che calcola l'attenuazione totale 
     *  @details Ricevuto in ingresso il dato di Z in byte e l'attenuazione complessiva sul raggio fino al punto in considerazione, calcola l'attenuazione totale 
     *  ingresso:dbz in quel punto e attenuazione fin lì
     *  Doviak,Zrnic,1984 for rain as reported in cost 717 final document
     *  @param[in]  DBZbyte valore in byte della Z nel pixel
     *  @param[in]  PIA attenuazione totale fino al punto
     *  @return att_tot attenuazione incrementata del contributo nel pixel corrente
     */
    double attenuation(unsigned char DBZbyte, double  PIA);
    double attenuation(double DBZvalue, double  PIA);   /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/

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

