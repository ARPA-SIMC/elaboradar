/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_ALGO_DBZ_H
#define RADARELAB_ALGO_DBZ_H

#include <cmath>

namespace radarelab {
template<typename T> class Volume;

namespace algo {

/**
 * Class to manage reflectivity functions (simply attenuation correction, conversion between Z, dBZ, R)
 */
class DBZ
{
protected:
    void init(int month, double base_cell_size);

public:
    double base_cell_size;              ///< cella size dimension
    double aMP, bMP;                    ////< Marshall-Palmer coefficient for Z-R relationship

    DBZ(const Volume<double>& volume);

    /**
     * @brief Seasonal setup function
     * @param [in] month - month 
     * @param [in] base_cell_size - cell size dimension [m]
     */
    DBZ(int month, double base_cell_size);

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

    /**
     * @function
     * Compute the corrected value (in dB) given the original value (in dB) and the
     * beam blocking percentage (from 0 to 100) for that value
     * @param [in] val_db - uncorrected dBZ value
     * @param [in] beamblocking - percentage of beam blocking
     * @return corrected dBZ value
     */
    static constexpr inline double beam_blocking_correction(double val_db, double beamblocking)
    {
       return val_db - 10 * log10(1. - beamblocking / 100.);
    }

    /**
     * @brief funzione che converte Z unsigned char in  DBZ
     * @param[in] DBZbyte riflettivita' in byte
     * @param [in] gain  - first conversion factor 
     * @param [in] offset - second conversion factor 
     * @return equivalente in DBZ
     */
    static constexpr inline double BYTEtoDB(unsigned char DBZbyte, double gain=80./255., double offset=-20.)
    {
        return (DBZbyte * gain + offset);
    }

    /**
     * @brief funzione che converte  dB in valore intero tra 0 e 255
     * @param[in] DB dBZ in ingresso
     * @param [in] gain  - first conversion factor 
     * @param [in] offset - second conversion factor 
     * @return converted value in byte
     */
    static inline unsigned char DBtoBYTE(double DB, double gain=80./255., double offset=-20.)
    {
        int byt = round((DB - offset) / gain);
        if (byt <= 0)
            return 0;
        else if (byt <= 255)
            return (unsigned char)(byt);
        else 
            return 255;
    }

    /**
     * @brief funzione che converte byte in Z
     * @param[in] byte valore da convertire in z espresso tra 0 e 255
     * @param [in] gain  - first conversion factor 
     * @param [in] offset - second conversion factor 
     * @return Z value (linear not dBZ)
     */
    static double BYTEtoZ(unsigned char byte);

    /**
     * @brief funzione che converte dBZ in Z
     * @param[in] DBZ valore da convertire in z 
     * @return Z value (linear not dBZ)
     */
    static constexpr inline double DBZtoZ(double DBZ) { return pow(10., DBZ * 0.1); }

    /**
     * @brief funzione che converte Z in dBZ
     * @param[in] Z valore da convertire in dB
     * @return dBZ value
     */
    static constexpr inline double ZtoDBZ(double Z) { return 10 * log10(Z); }

    /**
     * @brief funzione che converte dbZ  in R usando a e b variabili
     * @param[in] dbz  dB da convertire
     * @param[in] aMP  a della relazione z-r
     * @param[in] bMP  b della relazione z-r
     * @return R [mmh-1]
     */
    static constexpr inline double DBZtoR( double dbz, double aMP, double bMP) { return pow(pow(10., dbz / 10.) / aMP,1. / bMP); }

    /**
     * @brief funzione che converte R in  dbZ usando a e b variabili
     * @param[in]  rain  tasso di pioggia
     * @param[in] aMP  a della relazione z-r
     * @param[in] bMP  b della relazione z-r
     * @return dBZ 
     */
    static constexpr inline double RtoDBZ( double rain, double aMP, double bMP) { return 10. * (log10(aMP * pow(rain, bMP))); }

    /**
     * @brief funzione che converte R in  Z usando a e b variabili
     * @param[in]  rain  tasso di pioggia
     * @param[in] aMP  a della relazione z-r
     * @param[in] bMP  b della relazione z-r
     * @return Z value (linear not dBZ)
     */
    static constexpr inline double RtoZ( double rain, double aMP, double bMP) { return aMP * pow(rain, bMP); }

    /**
     * @brief funzione che converte Z  in R usando a e b variabili
     * @param[in] z      da convertire
     * @param[in] aMP  a della relazione z-r
     * @param[in] bMP  b della relazione z-r
     * @return R value [mmh-1]
     */
    static constexpr inline double ZtoR( double z, double aMP, double bMP) { return  pow(z / aMP, 1. / bMP); }

    /**
     * @brief funzione che converte byte in R usando a e b variabili
     * @param[in] byte byte da convertire
     * @param[in] aMP  a della relazione z-r
     * @param[in] bMP  b della relazione z-r
     * @param [in] gain  - first conversion factor 
     * @param [in] offset - second conversion factor 
     * @return R [mmh-1]
     */
    static constexpr inline double BYTE_to_mp_func(unsigned char byte, double aMP, double bMP, double gain=80./255., double offset=-20.)
    {
       return pow(pow(10., 0.1 * (byte * gain + offset)) / aMP, 1. / bMP);
    }
};

}
}

#endif

