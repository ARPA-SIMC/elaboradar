/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_ALGO_CLEANER_H
#define RADARELAB_ALGO_CLEANER_H

#include <radarelab/volume.h>
#include <radarelab/loader.h>

#include <fstream>
#include <string>

using namespace Eigen;
using namespace std;
namespace radarelab {
namespace algo {
/**
 * Struttura per cleaner dati grezzi sulla base dei valori di V, W e la deviazione standard di Z.
 */
struct Cleaner
{
    const unsigned min_segment_length =  2;		///< lunghezza minima segmento in celle 
    const unsigned max_segment_length = 40;		///< lunghezza massima segmento in celle se piÃ¹ lungo pulisce in ogni caso

    const double Z_missing;				///< Valore dato mancante DBZH.
    const double W_threshold;				///< Soglia per WRAD. 
    const double V_missing;				///< Dato mancante per VRAD. 
    const double bin_wind_magic_number;			///< valore magico per dati in formato SP20
    const double sd_threshold = 2;			///< Soglia per devizione standard DBZH.

/// Constructor
    Cleaner(double Z_missing, double W_threshold, double V_missing, double bin_wind_magic_number)
        : Z_missing(Z_missing), W_threshold(W_threshold), V_missing(V_missing), bin_wind_magic_number(bin_wind_magic_number) {}
/**
 * Funzione per ripulire raggio.Utilizza (sigmaV, V)
 * @param [in]	beam_z	- raggio DBZH
 * @param [in]	beam_w	- raggio WRAD
 * @param [in]	beam_v	- raggio VRAD
 * @return raggio di valori boleani per correzione (true-> Da Correggere)
 */
    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v,int i) const;
/**
 * Funzione per ripulire raggio.Utilizza (sigmaV, V) 
 * Analoga a  clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v,int i), ma 
 * restituisce un vettore di unsigned char 
 * @param [in]	beam_z	- raggio DBZH
 * @param [in]	beam_w	- raggio WRAD
 * @param [in]	beam_v	- raggio VRAD
 * @return raggio di valori unsigned char per correzione (1 -> Da Correggere)
 */
    std::vector<unsigned char> eval_clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v,int i) const;

/**
 * Funzione per ripulire raggio. Utilizza (sigmaV, V, dev.std Z, dev.std. ZDR)
 * @param [in]	beam_z	- raggio DBZH
 * @param [in]	beam_w	- raggio WRAD
 * @param [in]	beam_v	- raggio VRAD
 * @param [in]	beam_sd	- raggio deviazione standard DBZH
 * @param [in]	beam_sdzdr- raggio deviazione standard ZDR
 * @param [in]	scan_z	- per debug
 * @param [in]	scan_w	- per debug
 * @param [in]	scan_v	- per debug
 * @param [in]	SD	- per debug
 * @param [in]	iray	- index of the ray per debug
 * @return raggio di flag per correzione
 */
    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdzdr, PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& SD,int iray) const;

/**
 * Funzione per ripulire raggio. Utilizza (sigmaV, V, dev.std Z)
 * @param [in]	beam_z	- raggio DBZH
 * @param [in]	beam_w	- raggio WRAD
 * @param [in]	beam_v	- raggio VRAD
 * @param [in]	beam_sd	- raggio deviazione standard DBZH
 * @param [in]	scan_z	- per debug
 * @param [in]	scan_w	- per debug
 * @param [in]	scan_v	- per debug
 * @param [in]	SD	- per debug
 * @param [in]	iray	- index of the ray per debug
 * @return raggio di flag per correzione
 */
    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& SD,int iray) const;
    std::vector<unsigned char> eval_clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, int iray) const;
  std::vector<unsigned char> eval_classID_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_zdr, const Eigen::VectorXd& beam_rohv, const Eigen::VectorXd& beam_sqi, const Eigen::VectorXd& beam_snr, const Eigen::VectorXd& beam_zvd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, const Eigen::VectorXd& beam_zdr_sd, int iray, const string radar, double v_ny) const;

/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	iel	- indice elevazione solo per debug
 */
  std::vector<unsigned char> eval_classID_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, int iray, const string radar, double v_ny) const;
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V,unsigned iel=0, bool set_undetect=false);
  static void evaluateCleanID(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V,PolarScan<unsigned char>& scan_cleanID, unsigned iel=0);
/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	bin_wind_magic_number	- soglia vento per dati in formato SP20
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, double bin_wind_magic_number,unsigned iel=0, bool set_undetect=false);

/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	scan_ZDR - volume di ZDR
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, PolarScan<double>& scan_ZDR,unsigned iel=0, bool set_undetect=false);

    static void evaluateCleanID(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V,PolarScan<unsigned char>& scan_cleanID, double bin_wind_magic_number, unsigned iel=0);

  static void evaluateClassID(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, PolarScan<double>& scan_zdr, PolarScan<double>& scan_rohv, PolarScan<double>& scan_sqi, PolarScan<double>& scan_snr, PolarScan<double>& scan_zvd, PolarScan<unsigned char>& scan_cleanID,double bin_wind_magic_number, const string radar, unsigned iel=0);
/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	scan_zdr - volume di ZDR
 * @param [in]	bin_wind_magic_number	- soglia vento per dati in formato SP20
 * @param [in]	iel	- indice elevazione solo per debug
 */
  static void evaluateClassID(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, PolarScan<unsigned char>& scan_cleanID, double bin_wind_magic_number, const string radar, unsigned iel=0);
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, PolarScan<double>& scan_ZDR, double bin_wind_magic_number,unsigned iel=0, bool set_undetect=false);

/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	loader 	- Struttura di caricamento dati 
 * @param [in]	bin_wind_magic_number	- soglia vento per dati in formato SP20
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean( radarelab::volume::Loader load_structure, double bin_wind_magic_number,unsigned iel=0, bool set_undetect=false);
  
/*!
 * trapezoidal probability function
 */
  double trap(double x1, double x2, double x3, double x4, double val, double x5=-9999.) const;

  /*function reading matrix from txt file*/
  vector<string> read_matrix_from_txt(string fin) const;
  
};

}
}

#endif
