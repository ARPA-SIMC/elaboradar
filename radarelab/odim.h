#ifndef RADARELAB_ODIM_H
#define RADARELAB_ODIM_H
/**
 *  @file
 *  @ingroup radarelab
 *  @brief Codice per il caricamento di volumi ODIM in radarelab
*/


#include <radarelab/volume.h>
#include <string>
#include <map>
#include <vector>
#include <radarelab/loader.h>

namespace radarelab {
namespace volume {

/**
 *  Struttura che eredita da Loader e definisce i metodi per accedere ai dati ODIM
 */
struct ODIMLoader : Loader
{
    /**
     * Define a request - Fill to_load attribute  
     * @param [in] name   - Quantity requested
     * @param [in] volume - Scans where the data will be loaded
     */
    void request_quantity(const std::string& name, Scans<double>* volume);

    /**
     * Load method 
     * @param [in] pathname - full path for data file
     */
    void load(const std::string& pathname);

    /**
     *  Get list of available elevation
     * 
     */
    std::vector<double> get_nominal_elevations() { return Available_Elevations; }

    std::vector<double> Available_Elevations;

};

/*
union bit64
{
	double fp64;
	int64_t int64;
};
*/

class ODIMStorer
{
public:
    std::vector<Volume<double>*> to_store_fp;
    std::vector<Volume<int>*> to_store_int;
    std::vector<Volume<unsigned char>*> to_store_uchar;
    void store_quantity_fp(Volume<double>* vol_fp) {to_store_fp.push_back(vol_fp);}
    void store_quantity_int(Volume<int>* vol_int) {to_store_int.push_back(vol_int);}
    void store_quality_uchar(Volume<unsigned char>* vol_uchar) {to_store_uchar.push_back(vol_uchar);}
    void store(const std::string& pathname);
    void storeQuality(const std::string& pathname, const std::string& task, bool RemoveQualityFields=true);
};

} // volume
} // radarelab

#endif

