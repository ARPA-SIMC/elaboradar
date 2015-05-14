#ifndef ELABORADAR_ALGO_CLEANER_H
#define ELABORADAR_ALGO_CLEANER_H

#include <elaboradar/volume.h>
#include <elaboradar/loader.h>

namespace elaboradar {
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
 * @return raggio di flag per correzione
 */
    std::vector<bool> clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v,int i) const;

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

/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V,unsigned iel=0);
/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	scan_ZDR - volume di ZDR
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, PolarScan<double>& scan_ZDR,unsigned iel=0);
/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	bin_wind_magic_number	- soglia vento per dati in formato SP20
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, double bin_wind_magic_number,unsigned iel=0);
/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	scan_Z	- volume di DBZH
 * @param [in]	scan_W	- volume di WRAD
 * @param [in]	scan_V	- volume di V
 * @param [in]	scan_zdr - volume di ZDR
 * @param [in]	bin_wind_magic_number	- soglia vento per dati in formato SP20
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean(PolarScan<double>& scan_Z, PolarScan<double>& scan_W, PolarScan<double>& scan_V, PolarScan<double>& scan_ZDR, double bin_wind_magic_number,unsigned iel=0);

/**
 * Funzione che crea l'oggetto cleaner, lo inizializza, pulisce i dati e modifica il PolarScan di DBZH.
 * @param [in,out]	loader 	- Struttura di caricamento dati 
 * @param [in]	bin_wind_magic_number	- soglia vento per dati in formato SP20
 * @param [in]	iel	- indice elevazione solo per debug
 */
    static void clean( elaboradar::volume::Loader load_structure, double bin_wind_magic_number,unsigned iel=0);
};

}
}

#endif
