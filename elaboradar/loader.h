#ifndef SRC_LOADER_H
#define SRC_LOADER_H

/**
 *  @file
 *  @ingroup progetto_cum_bac
 *  @brief Codice per il caricamento di volumi in elaboradar
*/
#include <elaboradar/volume.h>
#include <string>
#include <map>
#include <vector>

namespace elaboradar {
namespace volume {

/**
 *  Struttura che contiene mappa per caricamento dati
 */
struct Loader
{
    std::map<std::string, Scans<double>*> to_load;  ///< Map used to specify quantity to be loaded

};


}
}

#endif 

