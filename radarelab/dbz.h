/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_ALGO_UTILS_H
#define RADARELAB_ALGO_UTILS_H

#include <math.h>

namespace radarelab {
namespace algo {

/**
 * @file algo/utils.h
 */

/**
*  @brief funzione che converte dBZ in Z
*  @param[in] DBZ valore da convertire in z 
*  @return Z value (linear not dBZ)
*/
inline double DBZtoZ(double DBZ){ return pow(10.,DBZ*0.1);  } 

/**
*  @brief funzione che converte Z in dBZ
*  @param[in] Z valore da convertire in dB
*  @return dBZ value
*/
inline double ZtoDBZ(double Z) { return 10 * log10(Z); }

/**
*  @brief funzione che converte byte in R usando a e b variabili
*  @param[in] byte byte da convertire
*  @param[in] aMP  a della relazione z-r
*  @param[in] bMP  b della relazione z-r
*  @param [in] gain  - first conversion factor 
*  @param [in] offset - second conversion factor 
*  @return R [mmh-1]
*/
inline double BYTE_to_mp_func( unsigned char byte,  double aMP,  double bMP, double gain=80./255., double offset=-20.){ 
   return pow(pow(10.,0.1*(byte*gain+offset))/aMP,1./bMP);
}
/**
*  @brief funzione che converte R in  dbZ usando a e b variabili
*  @param[in]  rain  tasso di pioggia
*  @param[in] aMP  a della relazione z-r
*  @param[in] bMP  b della relazione z-r
*  @return dBZ 
*/
inline double RtoDBZ( double rain, double aMP, double bMP){ return 10.*( log10(aMP*pow(rain,bMP)) ) ; }

/**
*  @brief funzione che converte dbZ  in R usando a e b variabili
*  @param[in] dbz  dB da convertire
*  @param[in] aMP  a della relazione z-r
*  @param[in] bMP  b della relazione z-r
*  @return R [mmh-1]
*/
inline double DBZtoR( double dbz, double aMP, double bMP) { return  pow(pow(10., dbz/10.)/aMP,1./bMP); }

/**
*  @brief funzione che converte R in  Z usando a e b variabili
*  @param[in]  rain  tasso di pioggia
*  @param[in] aMP  a della relazione z-r
*  @param[in] bMP  b della relazione z-r
*  @return Z value (linear not dBZ)
*/
inline double RtoZ( double rain, double aMP, double bMP){ return aMP*pow(rain,bMP);} 
     
/**
*  @brief funzione che converte Z  in R usando a e b variabili
*  @param[in] z      da convertire
*  @param[in] aMP  a della relazione z-r
*  @param[in] bMP  b della relazione z-r
*  @return R value [mmh-1]
*/
inline double ZtoR( double z, double aMP, double bMP) { return  pow(z/aMP,1./bMP); }

}
}

#endif

