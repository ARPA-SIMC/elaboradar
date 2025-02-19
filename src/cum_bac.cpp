#include "cum_bac.h"
#include <radarelab/logging.h>
#include <radarelab/utils.h>
#include <radarelab/image.h>
#include <radarelab/sp20.h>
#include <radarelab/odim.h>
#include <radarelab/algo/cleaner.h>
#include <radarelab/algo/anaprop.h>
#include <radarelab/algo/azimuth_resample.h>
#include <radarelab/algo/dbz.h>
#include <radarelab/algo/vpr.h>
#include "site.h"
#include <radarelab/algo/steiner.h>
#include <radarelab/algo/viz.h>
#include <radarelab/algo/elabora_volume.h>
#include "cartproducts.h"
#include <radarelab/algo/top.h>
#include <radarelab/cylindrical.h>
#include <radarelab/interpola_vpr.h>
#include <radarelab/cart.h>
#include <radarlib/radar.hpp>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <unistd.h>
#include "setwork.h"
#include <sstream>

//#ifdef __cplusplus
//extern "C" {
//#endif
//#include <func_Z_R.h>
//#ifdef __cplusplus
//}
//#endif

#include <func_Q3d.h>

#include <qual_par.h>
#include <radarelab/par_class.h>

#ifdef NEL
#undef NEL
#endif

// Soglie algoritmi
#define OVERBLOCKING 51 /* minimo BB non accettato*/
#define SOGLIA_TOP 45.0 // soglia per trovare top
#define MISSING 0 /*valore mancante*/

//Definizioni geometriche
#define AMPLITUDE 0.9 /* esternalizzo?*/ // ampiezza fascio radar

// anaprop
#define LIMITE_ANAP 240/* LIMITE in numero bins per cambiare controllo anaprop*/

#define DTOR  M_PI/180. /* esternalizzo?*/ //fattore conversione gradi-radianti
#define CONV_RAD 360./4096.*DTOR  // fattore conversione unità angolare radar-radianti

using namespace std;
using namespace radarelab;
using namespace radarelab::algo;

namespace elaboradar {

namespace {
/**
 *  @brief funzione che legge la quota del centro fascio e del limite inferiore del fascio da file
 *  @details legge la quota del centro fascio e del limite inferiore del fascio da file e li memorizza nei vettori hray_inf e hray
 *  @return 
 */
struct HRay : public Matrix2D<double>
{
    static const int NSCAN = 6;

    // distanza temporale radiosondaggio
    double dtrs;

    template<typename T>
    HRay(const Volume<T>& vol) : Matrix2D<double>(vol.size(), vol.max_beam_size())
    {
        const double radius = 6378137.0;
        const double kea = 4. / 3. * radius;

        for (unsigned iel = 0; iel < vol.size(); ++iel)
        {
            const double elev = vol.scan(iel).elevation;
            const double cell_size = vol.scan(iel).cell_size;

            for (unsigned ibin = 0; ibin < cols(); ++ibin)
            {
                double range = (ibin + 0.5) * cell_size;
                (*this)(iel, ibin) = ::sqrt(::pow(range, 2.) + ::pow(kea, 2.) + 2 * kea * range * ::sin(DTOR * elev)) - kea;
            }
        }
    }

    void load_hray(Assets& assets)
    {
        // quota centro fascio in funzione della distanza e elevazione
        dtrs = assets.read_file_hray([this](unsigned el, unsigned bin, double value) {
            if (el >= rows() || bin >= cols()) return;
            (*this)(el, bin) = value;
        });
    }
    void load_hray_inf(Assets& assets)
    {
        // quota limite inferiore fascio in funzione della distanza e elevazione
        dtrs = assets.read_file_hray_inf([this](unsigned el, unsigned bin, double value) {
            if (el >= rows() || bin >= cols()) return;
            (*this)(el, bin) = value;
        });
    }
};
}

CUM_BAC::CUM_BAC(Volume<double>& volume, const Config& cfg, const Site& site, bool medium, unsigned max_bin)
    : MyMAX_BIN(max_bin), cfg(cfg), site(site), assets(cfg),
      do_medium(medium), volume(volume), SD_Z6(volume.beam_count), cil(volume, 0, RES_HOR_CIL, RES_VERT_CIL),
      dbz(volume), flag_vpr(volume, 0),
      first_level(NUM_AZ_X_PPI, MyMAX_BIN), first_level_static(NUM_AZ_X_PPI, MyMAX_BIN),
      bb_first_level(NUM_AZ_X_PPI, 1024), beam_blocking(NUM_AZ_X_PPI, 1024),
      anaprop(volume), dem(NUM_AZ_X_PPI, 1024),
      qual(volume, 0),
      top(volume.beam_count, volume.max_beam_size())
{
    logging_category = log4c_category_get("radar.cum_bac");
    assets.configure(site, volume.load_info->acq_date);

    //definisco stringa data in modo predefinito
    time_t Time = volume.load_info->acq_date;
    struct tm* tempo = gmtime(&Time);
    // scrivo la variabile char date con la data in formato aaaammgghhmm
    strftime(date, 20, "%Y%m%d%H%M", tempo);

    // ------definisco i coeff MP in base alla stagione( mese) che servono per calcolo VPR e attenuazione--------------
    algo::compute_top(volume, SOGLIA_TOP, top);
}

CUM_BAC::~CUM_BAC()
{
    if (calcolo_vpr) delete calcolo_vpr;
}

void CUM_BAC::want_vpr()
{
    calcolo_vpr = new CalcoloVPR(*this);
}

void CUM_BAC::read_sp20_volume(Volume<double>& volume, const Site& site, const char* nome_file, bool do_clean, bool do_medium)
{
    using namespace radarelab::volume;
    LOG_CATEGORY("radar.io");
    LOG_INFO("Reading %s for site %s", nome_file, site.name.c_str());

    SP20Loader loader;

    Scans<double> z_volume;
    Scans<double> w_volume;
    Scans<double> v_volume;
    loader.vol_z = &z_volume;
    loader.vol_w = &w_volume;
    loader.vol_v = &v_volume;
    loader.load(nome_file);

    // Normalise the scan elevations to match the elevations requested in Site
    auto elev_array = site.get_elev_array(do_medium);
    z_volume.normalize_elevations(elev_array);
    w_volume.normalize_elevations(elev_array);
    v_volume.normalize_elevations(elev_array);

    if (do_clean)
    {
        for (unsigned i = 0; i < z_volume.size(); ++i)  {
            double bin_wind_magic_number = site.get_bin_wind_magic_number(v_volume.load_info->acq_date)
                * v_volume.at(i).gain + v_volume.at(i).offset;
            algo::Cleaner::clean(z_volume.at(i), w_volume.at(i), v_volume.at(i), bin_wind_magic_number);
        }
    }

    algo::azimuthresample::MaxOfClosest<double> resampler;
    resampler.resample_volume(z_volume, volume, 1.0);
// Copy all radar site information to volume data
    volume.radarSite = site.radarSite;

}

  void CUM_BAC::read_odim_volume(Volume<double>& volume, const Site& site, const char* nome_file, char* fuzzypath, bool do_clean, bool do_medium, bool set_undetect)
{
    using namespace radarelab::volume;
    LOG_CATEGORY("radar.io");
    namespace odim = OdimH5v21;
    LOG_INFO("Reading %s for site %s", nome_file, site.name.c_str());

    volume::ODIMLoader loader;

    Scans<double> dbzh_volume;
    Scans<double> th_volume;
    Scans<double> v_volume;
    Scans<double> w_volume;
    Scans<double> zdr_volume;
    Scans<double> rhohv_volume;
    Scans<double> sqi_volume;
    Scans<double> snr_volume;
    //Scans<unsigned char> full_volume_cleanID;
    //Scans<double> full_volume_diffprob;
    //full_volume_diffprob.quantity="Diffprob";

    string radar_name = site.name.c_str(); // da elaboradar/src/site.cpp : 'SPC' o 'GAT'
    bool init_sqi = false;
    
    loader.request_quantity(odim::PRODUCT_QUANTITY_DBZH, &dbzh_volume);
    loader.request_quantity(odim::PRODUCT_QUANTITY_TH, &th_volume);

    if (do_clean)
    {
        loader.request_quantity(odim::PRODUCT_QUANTITY_VRAD, &v_volume);
        loader.request_quantity(odim::PRODUCT_QUANTITY_WRAD, &w_volume);
        loader.request_quantity(odim::PRODUCT_QUANTITY_ZDR, &zdr_volume);
	loader.request_quantity(odim::PRODUCT_QUANTITY_RHOHV,&rhohv_volume);
	loader.request_quantity(odim::PRODUCT_QUANTITY_SQI,&sqi_volume);
	loader.request_quantity(odim::PRODUCT_QUANTITY_SNR,&snr_volume);
    }
    loader.load(nome_file);

    // FIXME: are they really empty? isn't make_scan called on all of them?
    if (dbzh_volume.empty() && th_volume.empty())
    {
        LOG_ERROR("neither DBZH nor TH were found in %s", nome_file);
        throw runtime_error("neither DBZH nor TH were found");
    }

    // Normalise the scan elevations to match the elevations requested in Site
    auto elev_array = site.get_elev_array(do_medium);
    for (auto i: loader.to_load){
      if(!i.second->empty() ) i.second->normalize_elevations(elev_array);
    }
    Scans<double>* z_volume;
    if (!dbzh_volume.empty()) {
        LOG_WARN(" DBZH found");
       z_volume = &dbzh_volume;
    }
    else {
        LOG_WARN("no DBZH found: using TH");
        z_volume = &th_volume;
    }

    //cout<<"z max ="<<z_volume->at(0).maxCoeff()<<endl;
    if (do_clean && !w_volume.empty() && !v_volume.empty())
    {
      if(sqi_volume.empty()) init_sqi = true;
      //cout<<"init_sqi"<<init_sqi<<endl;
      if (zdr_volume.empty())
      {
	//caso ZDR e grandezze polarimetriche assenti -> lascio versione master al 16/2/2022
	volume::Scans<unsigned char> full_volume_cleanID;
        //for (unsigned i = 0; i < 1; ++i){
        for (unsigned i = 0; i < z_volume->size(); ++i){
//            radarelab::algo::Cleaner::clean(z_volume->at(i), w_volume.at(i), v_volume.at(i),i);
	  full_volume_cleanID.append_scan(z_volume->at(i).beam_count,z_volume->at(i).beam_size,z_volume->at(i).elevation, z_volume->at(i).cell_size, 0);
            radarelab::algo::Cleaner::evaluateCleanID(z_volume->at(i), w_volume.at(i), v_volume.at(i),full_volume_cleanID.at(i),i);
            for (unsigned ibeam=0;ibeam<z_volume->at(i).beam_count; ibeam++)
		for (unsigned j=0; j<z_volume->at(i).beam_size; j++){
		  if (full_volume_cleanID.at(i)(ibeam,j) != 0) z_volume->at(i)(ibeam,j)=z_volume->at(i).undetect;
	        }
	}
      } else {
	cout<<"applico logica fuzzy"<<endl;
	//caso ZDR e grandezze polarimetriche presenti--> uso fuzzy logic del branch elaboradar_updating al 16/2/2022
	//volume::Scans<unsigned char> full_volume_cleanID;
	//volume::Scans<double>full_volume_diffprob;
	volume::Scans<unsigned char> full_volume_cleanID;

	unsigned last = z_volume->size() -1;
        for (unsigned i = 0; i < z_volume->size(); ++i){
	  //cout<<"it="<<i<<endl;
	    //volume::Scans<unsigned char> full_volume_cleanID;
	    volume::Scans<double>full_volume_diffprob;
	    full_volume_cleanID.append_scan(z_volume->at(i).beam_count,z_volume->at(i).beam_size,z_volume->at(i).elevation, z_volume->at(i).cell_size);
	    full_volume_diffprob.append_scan(z_volume->at(i).beam_count,z_volume->at(i).beam_size,z_volume->at(i).elevation, z_volume->at(i).cell_size);
	    //full_volume_cleanID.at(0).setZero();

	    volume::Scans<double> Texture;
	    //cout<<"init_sqi = "<<init_sqi<<endl;
	    if(init_sqi){
	      sqi_volume.append_scan(z_volume->at(i).beam_count,z_volume->at(i).beam_size,z_volume->at(i).elevation, z_volume->at(i).cell_size);
	      sqi_volume.at(i).setZero();
	    }

            //calcolo texture di dbzh sulla verticale tra prima elevazione e seconda elevazione:
	    if(i< last){
                volume::Scans<double> Input,Input2;        
		Input.push_back(z_volume->at(i));
	        Input2.push_back(z_volume->at(i+1));
	        radarelab::volume::textureVD(Input, Input2, Texture, true);
	        Texture.at(0).nodata=65535.;
	        Texture.at(0).undetect=0.;
		//cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
	    }
	    else{
	        Texture.clear();
		//cout<<"it="<<i<<", Texture size = "<<Texture.size()<<endl;
		Texture.append_scan(z_volume->at(i).beam_count,z_volume->at(i).beam_size,z_volume->at(i).elevation, z_volume->at(i).cell_size);
		Texture.at(0).setZero();
		Texture.at(0).nodata=65535.;
	        Texture.at(0).undetect=0.;
		//cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
	    }
	    
	    //algo::Cleaner::evaluateClassID(z_volume->at(i), w_volume.at(i), v_volume.at(i), zdr_volume.at(i), rhohv_volume.at(i), sqi_volume.at(i), snr_volume.at(i), Texture.at(0), full_volume_cleanID.at(i), v_volume.at(i).undetect , radar_name, i);
	    //modifico il force_bool a true (ultimo parametro)
	    
	    algo::Cleaner::evaluateClassID(z_volume->at(i), w_volume.at(i), v_volume.at(i), zdr_volume.at(i), rhohv_volume.at(i), sqi_volume.at(i), snr_volume.at(i), Texture.at(0), full_volume_cleanID.at(i), full_volume_diffprob.at(0), v_volume.at(i).undetect , radar_name, fuzzypath, i, true);

            double new_value=z_volume->at(0).nodata;
	    //provo a commentare per test di uso solo nodata dei punti clssificati come non meteo
	    if (set_undetect) new_value=z_volume->at(i).undetect;

	    for (unsigned ii = 0; ii < z_volume->at(i).beam_count; ++ii)
                for (unsigned ib = 0; ib < z_volume->at(i).beam_size; ++ib) {
		  
     	          if(full_volume_cleanID.at(i)(ii,ib) ) z_volume->at(i)(ii,ib)= new_value;

        	}	      
	    }
	    
	    // commento doppia ripulitura tramite clean() della versione master al 16/2/2022
            //algo::Cleaner::clean(z_volume->at(i), w_volume.at(i), v_volume.at(i),zdr_volume.at(i),i,true);
            //algo::Cleaner::clean(z_volume->at(i), w_volume.at(i), v_volume.at(i),zdr_volume.at(i),i+100,true);
        }
      }

    //cout<<"arrivo al ponghino"<<endl;

    algo::azimuthresample::MaxOfClosest<double> resampler;
    resampler.resample_volume(*z_volume, volume, 1.0);
    cout<<"resampler fatto!!!"<<endl;

    /*
    printf("fbeam ϑ%f α%f", this->volume.scan(0)[0].teta, this->volume.scan(0)[0].alfa);
    for (unsigned i = 0; i < 20; ++i)
        printf(" %d", (int)this->volume.scan(0).get_raw(0, i));
    printf("\n");
    */

    /*
    int numRaggi»···»···= scan->getNumRays();
    NUM_AZ_X_PPI

    NEL

        se due scan per stessa elecvazione, prendo il primo

        guardare se il passo di griglia è 0.9 o dare errore
        sennò prendere il beam che ha l'angolo piú vicino

        fill_bin in sp20lib

        leggere DBZH o TH (fare poi DBtoBYTE)
        */

    /*
    struct VOL_POL volume.scan(NEL)[NUM_AZ_X_PPI];
    T_MDB_data_header   old_data_header;

    //--------lettura volume------
    int tipo_dati_richiesti = INDEX_Z;
    int ier = read_dbp_SP20((char*)nome_file,volume.vol_pol,&old_data_header,
                            tipo_dati_richiesti,volume.nbeam_elev);

    if (ier != OK)
        LOG_ERROR("Reading %s returned error code %d", nome_file, ier);

    //  ----- Test sul volume test_volume.......  --------
    if (!test_volume(file_type))
    {
        LOG_ERROR("test_volume failed");
        return false;
    }
    */

    // TODO: look for the equivalent of declutter_rsp and check its consistency
    // like in test_volume
}


void CUM_BAC::declutter_anaprop()
{
    //-------------leggo mappa statica ovvero first_level (funzione leggo_first_level)------------
    leggo_first_level();

    //-------------se definita qualita' leggo dem e altezza fascio (mi servono per calcolare qualità)
    if (do_quality)
    {
        assets.load_dem(dem);
        dem.resize_beams_and_propagate_last_bin(volume.max_beam_size());
    }

    //------------se definito DECLUTTER , non rimuovo anap e riscrivo  volume polare facedndo declutter solo con mappa statica.... ancora valido?

    if (do_declutter)
    {
        for(unsigned i=0; i<NUM_AZ_X_PPI; i++)
            for(unsigned k=0; k<volume[0].beam_size; k++)
            {
                //---assegno el_inf a mappa statica
                unsigned el_inf = first_level_static(i, k);
                //---ricopio valori a mappa statica sotto
                for(unsigned l=0; l<=el_inf; l++)
                {
                    // Enrico: cerca di non leggere/scrivere fuori dal volume effettivo
                    if (k >= volume[l].beam_size) continue;
                    if (k < volume[el_inf].beam_size)
                        volume[l].set(i, k, volume[el_inf].get(i, k));
                    else
                        volume[l].set(i, k, MISSING_DB);
                    //------------se definito BEAM BLOCKING e non definito BLOCNOCORR (OPZIONE PER non correggere il beam blocking a livello di mappa statica PUR SAPENDO QUANT'È)
                    if (do_beamblocking && do_bloccorr)
                        volume[l].set(i, k, algo::DBZ::beam_blocking_correction(volume[l].get(i, k), beam_blocking(i, k)));
                }
            }

        anaprop.init_elev_fin_static(volume, first_level_static);
        LOG_INFO("declutter_anaprop completed with static declutter");
    }

    //------------se non definito DECLUTTER inizio rimozione propagazione anomala al livello mappa dinamica e elaborazioni accessorie
    else if (do_anaprop)
    {
        /* 26-5-2004 : se sono alla 1 o successive elevazioni
           e range > 60 km cambio le soglie, in modo
           da evitare di riconoscere come anaprop una pioggia shallow
           Il criterio diventa: - se la differenza tra Z all'elevazione più bassa della
           corrente e la Z corrente è <10 dbZ allora
           rendo inefficaci i limiti di riconoscimento anaprop. */

        //--------ciclo sugli azimut e bins per trovare punti con propagazione anomala----------------

        textureSD(volume,SD_Z6,6000., false);

        // test to define the more appropriate value for textture_threshold for rainy and clutter data
        std::vector <double>  above_0  (4,0);
        std::vector <double>  above_15 (4,0);
        std::vector <double>  above_30 (4,0);
        std::vector <double>  above_40 (4,0);
        for( unsigned el =0; el <4; el++)
        for (unsigned iaz=0; iaz<volume[el].beam_count; iaz++)
            for (unsigned k=80; k < volume[el].beam_size; k ++)
                if (volume[el].get(iaz,k) > 40.){
                        above_0[el]++;
                        above_15[el]++;
                        above_30[el]++;
                        above_40[el]++;
                } else if (volume[el].get(iaz,k) > 30.){
                        above_0[el]++;
                        above_15[el]++;
                        above_30[el]++;
                } else if (volume[el].get(iaz,k) > 15.){
                        above_0[el]++;
                        above_15[el]++;
                } else if (volume[el].get(iaz,k) > 0.){
                        above_0[el]++;
                }

        anaprop.do_quality = do_quality;
        anaprop.do_beamblocking = do_beamblocking;
        anaprop.do_bloccorr = do_bloccorr;
        if ( above_15[2]/above_15[0] >= 0.025){
           if (above_0[1]/above_0[0] >= 0.6 && above_30[2]/above_15[2] <0.15 && above_0[1] >=50000){
              anaprop.conf_texture_threshold = 5.;
        LOG_WARN("TEXTURE THRESHOLD USED %4.1f -- 0. %6d %6d %6d %6d -- 15. %6d %6d %6d %6d -- 30. %6d %6d %6d %6d -- 40. %6d %6d %6d %6d", anaprop.conf_texture_threshold, (int)above_0[0], (int)above_0[1], (int)above_0[2], (int)above_0[3], (int)above_15[0], (int)above_15[1], (int)above_15[2], (int)above_15[3], (int)above_30[0], (int)above_30[1], (int)above_30[2], (int)above_30[3], (int)above_40[0], (int)above_40[1], (int)above_40[2], (int)above_40[3] );
              anaprop.remove(volume, beam_blocking, first_level, first_level_static, SD_Z6);
           } else { 
        //      anaprop.conf_texture_threshold = 5.;
              anaprop.remove_without_SD(volume, beam_blocking, first_level, first_level_static,SD_Z6);
        LOG_WARN("THUNDERSTORM           %4.1f -- 0. %6d %6d %6d %6d -- 15. %6d %6d %6d %6d -- 30. %6d %6d %6d %6d -- 40. %6d %6d %6d %6d", -9.9, (int)above_0[0], (int)above_0[1], (int)above_0[2], (int)above_0[3], (int)above_15[0], (int)above_15[1], (int)above_15[2], (int)above_15[3], (int)above_30[0], (int)above_30[1], (int)above_30[2], (int)above_30[3], (int)above_40[0], (int)above_40[1], (int)above_40[2], (int)above_40[3] );
           }
        } else {
              anaprop.remove(volume, beam_blocking, first_level, first_level_static, SD_Z6);
        LOG_WARN("TEXTURE THRESHOLD USED %4.1f -- 0. %6d %6d %6d %6d -- 15. %6d %6d %6d %6d -- 30. %6d %6d %6d %6d -- 40. %6d %6d %6d %6d", anaprop.conf_texture_threshold, (int)above_0[0], (int)above_0[1], (int)above_0[2], (int)above_0[3], (int)above_15[0], (int)above_15[1], (int)above_15[2], (int)above_15[3], (int)above_30[0], (int)above_30[1], (int)above_30[2], (int)above_30[3], (int)above_40[0], (int)above_40[1], (int)above_40[2], (int)above_40[3] );
	
        }
        LOG_INFO("declutter_anaprop completed with anaprop");
        ScrivoStatistica(anaprop.grid_stats);
    }
    else
    {
        LOG_WARN("declutter_anaprop completed without doing anything");
    }

    //---------------------------- Code to plot data from polarMatrix
    /*   Image <unsigned char> toBePlotted (volume[0].beam_size, volume[0].beam_count);
       for(unsigned i=0; i<volume[0].beam_count; i++)
          for(unsigned k=0 ; k<volume[0].beam_size; k++){
             toBePlotted(i,k)= DBtoBYTE(volume[0].get(i, k));
        }
            radarelab::write_image(toBePlotted, "/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/Polarplot.png", "PNG");*/
        LOG_INFO("elabora_Dato completata");
}

void CUM_BAC::leggo_first_level()
{
    if (do_readStaticMap)
    {
        // Leggo mappa statica
      assets.load_first_level(first_level_static);
        // Allargo per coprire la dimensione del volume
        if (first_level_static.cols() < volume.max_beam_size())
            first_level_static.resize_beams_and_propagate_last_bin(volume.max_beam_size());

        // copio mappa statica su matrice first_level
        first_level = first_level_static;
        LOG_INFO("Letta mappa statica");
    }

    if (do_beamblocking)
    {
        // Leggo file elevazioni per BB
        assets.load_first_level_bb_el(bb_first_level);
        // Allargo per coprire la dimensione del volume
        if (bb_first_level.cols() < volume.max_beam_size())
            bb_first_level.resize_beams_and_propagate_last_bin(volume.max_beam_size());

        // Leggo file valore di BB
        assets.load_first_level_bb_bloc(beam_blocking);
        // Allargo per coprire la dimensione del volume
        if (beam_blocking.cols() < volume.max_beam_size())
            beam_blocking.resize_beams_and_propagate_last_bin(volume.max_beam_size());

        /* Se elevazione clutter statico < elevazione BB, prendi elevazione BB,
           altrimeti prendi elevazione clutter statico e metti a 0 il valore di BB*/
        for(unsigned i=0; i < first_level.rows(); ++i)
            for (unsigned j=0; j < first_level.cols(); ++j)
            {
                if (do_bloccorr)
                {
                    if (first_level_static(i, j)<=bb_first_level(i, j))
                        first_level(i, j)=bb_first_level(i, j);
                    else
                    {
                        beam_blocking(i, j)=0;
                        first_level(i, j)=first_level_static(i, j);
                    }
                } else {
                    if (first_level_static(i, j)>bb_first_level(i, j))
                        beam_blocking(i, j)=0;
                    if (first_level_static(i, j)<bb_first_level(i, j))
                        beam_blocking(i, j)=OVERBLOCKING;
                }
            }
    }

    /*-------------------------------
      patch per espandere il clutter
      -------------------------------*/
    if(do_medium){
        PolarScan<unsigned char> first_level_tmp(first_level);
        LOG_INFO(" Dentro patch %u ",MyMAX_BIN);
        for (unsigned i=NUM_AZ_X_PPI; i<800; i++)
            for (unsigned j=0; j<MyMAX_BIN; j++)
                for (unsigned k=i-1; k<i+2; k++)
                    if(first_level(i%NUM_AZ_X_PPI, j) < first_level_tmp(k%NUM_AZ_X_PPI, j))
                        first_level(i%NUM_AZ_X_PPI, j)=first_level_tmp(k%NUM_AZ_X_PPI, j);
        LOG_INFO(" fine patch %u ",MyMAX_BIN);
    }
}

void CUM_BAC::ScrivoStatistica(const algo::anaprop::GridStats& grid_stats)
{
    //Definizioni per statistica anap
    static const int DIM1_ST = 16;
    static const int DIM2_ST = 13;
    /*--- numero minimo di celle presenti in un
      settore per la statistica            ---*/
    static const int N_MIN_BIN = 500;

    int az,ran;
    unsigned char statistica[DIM1_ST][DIM2_ST];
    unsigned char statistica_bl[DIM1_ST][DIM2_ST];
    unsigned char statistica_el[DIM1_ST][DIM2_ST];

    memset(statistica,255,DIM1_ST*DIM2_ST);
    memset(statistica_bl,255,DIM1_ST*DIM2_ST);
    memset(statistica_el,255,DIM1_ST*DIM2_ST);

    for(az=0; az<DIM1_ST; az++)
        for(ran=0; ran<DIM2_ST; ran++){
            if (grid_stats.count(az, ran) >= N_MIN_BIN)
            {
                statistica[az][ran] = grid_stats.perc_anap(az, ran);
                statistica_bl[az][ran] = grid_stats.perc_bloc(az, ran);
                statistica_el[az][ran] = grid_stats.perc_elev(az, ran);
            }
        }
    File f_stat;

    if (f_stat.open_from_env("ANAP_STAT_FILE", "a"))
    {
        fwrite(date,12,1,f_stat);
        fwrite(statistica,DIM1_ST*DIM2_ST,1,f_stat);
    }

    if (f_stat.open_from_env("BLOC_STAT_FILE", "a"))
    {
        fwrite(date,12,1,f_stat);
        fwrite(statistica_bl,DIM1_ST*DIM2_ST,1,f_stat);
    }

    if (f_stat.open_from_env("ELEV_STAT_FILE", "a"))
    {
        fwrite(date,12,1,f_stat);
        fwrite(statistica_el,DIM1_ST*DIM2_ST,1,f_stat);
    }
    return ;
}

/*
   comstart caratterizzo_volume
   idx calcola qualita' volume polare
   calcola qualita' volume polare
   NB il calcolo è fatto considerando q=0 al di sotto della mappa dinamica.
   per ora drrs=dist nche nel caso di Gattatico, mentre dtrs è letto da file
   si puo' scegliere tra qualita' rispetto a Z e rispetto a R, in realtà per ora sono uguali.

   float quo:  quota bin
   float el:   elevazione
   float rst:  raggio equivalente in condizioni standard
   float dZ:   correzione vpr
   float sdevZ:  stand. dev. correzione vpr

//--- nb. non ho il valore di bb sotto_bb_first_level
comend
*/
void CUM_BAC::caratterizzo_volume()
{
    LOG_DEBUG("start caratterizzo_volume");

    HRay hray_inf(volume); /*quota limite inferiore fascio in funzione della distanza e elevazione*/
    hray_inf.load_hray_inf(assets);

    cout<<"sono in caratterizzo volume"<<endl;

    // path integrated attenuation
    double PIA;
    // dimensione verticale bin calcolata tramite approcio geo-ottico
    float dh=1.;
    // distanza radiosondaggio,
    float dhst=1.;
    // tempo dal radiosondaggio
    float drrs=1.;
    // distanza dal radar
    float dist=1.;
    // beam blocking
    unsigned char bb=0;
    // indice clutter da anaprop
    unsigned char cl=0;

    //----------ciclo su NSCAN(=6), cioè sul numero di elevazioni (nominali) per le quali ho calcolato il beam blocking
    /* a questo punto servono: bb, cl,  PIA, dtrs e drrs radiosond, quota, hsup e hinf beam-----------------*/

    for (unsigned l=0; l<volume.size(); l++)/*ciclo elevazioni*/// VERIFICARE CHE VADA TUTTO OK
    {
        const auto& scan = volume[l];
        for (int i=0; i<NUM_AZ_X_PPI; i++)/*ciclo azimuth*/
        {
            const double elevaz = scan.elevations_real(i); //--- elev reale in gradi

            //--assegno PIA=0 lungo il raggio NB: il ciclo nn va cambiato in ordine di indici!
            PIA=0.;

            for (unsigned k=0; k<scan.beam_size; k++)/*ciclo range*/
            {
                double sample = scan.get(i, k);

                //---------distanza in m dal radar (250*k+125 x il corto..)
                dist = k * scan.cell_size + scan.cell_size / 2.;/*distanza radar */

                //-----distanza dal radiosondaggio (per GAT si finge che sia colocato ..), perchè? (verificare che serva )
                drrs=dist;
                /* if (!(strcmp(sito,"GAT")) ) {  */
                /*     drrs=dist; */
                /* } */
                /* if (!(strcmp(sito,"SPC")) ) {  */
                /*     drrs=dist; */
                /* } */


                //assegno la PIA (path integrated attenuation) nel punto e POI la incremento  (è funzione dell'attenuazione precedente e del valore nel punto)
                PIA=dbz.attenuation(DBZ::DBtoBYTE(sample),PIA);

                //------calcolo il dhst ciè l'altezza dal bin in condizioni standard utilizzando la funzione quota_f e le elevazioni reali
                dhst = PolarScanBase::sample_height(elevaz + 0.45, dist)
                     - PolarScanBase::sample_height(elevaz - 0.45, dist);

                //----qui si fa un po' di mischione: finchè ho il dato dal programma di beam blocking uso il dh con propagazione da radiosondaggio, alle elevazioni superiori assegno dh=dhst  e calcolo quota come se fosse prop. standard, però uso le elevazioni nominali

                if (l<volume.size() -1   ) {
                    // differenza tra limite sup e inf lobo centrale secondo appoccio geo-ott
                    dh = hray_inf(l + 1, k) - hray_inf(l, k);
                }
                else {
                    // non ho le altezze oltre nscan-1 pero' suppongo che a tali elevazioni la prop. si possa considerare standard
                    dh = dhst;
                }

                if (l < anaprop.elev_fin(i, k)) {
                    cl=algo::ANAP_YES;
                    bb=BBMAX;
                } else if (l == anaprop.elev_fin(i, k)) {
                    cl=anaprop.dato_corrotto(i, k);  /*cl al livello della mappa dinamica*/
                    bb=beam_blocking(i, k);  /*bb al livello della mappa dinamica *///sarebbe da ricontrollare perchè con la copia sopra non è più così
                } else if (l > anaprop.elev_fin(i, k)) {
                    cl=0;       /*per come viene scelta la mappa dinamica si suppone che al livello superiore cl=0 e bb=0*/
                    bb=0;   // sarebbe if (l-bb_first_level(i, k) >0  bb=0;  sopra all'elevazione per cui bb<soglia il bb sia =0 dato che sono contigue o più però condiz. inclusa
                }

                //------dato che non ho il valore di beam blocking sotto i livelli che ricevo in ingresso ada progrmma beam blocking e
                //--------dato che sotto elev_fin rimuovo i dati come fosse anaprop ( in realtà c'è da considerare che qui ho pure bb>50%)
                //--------------assegno qualità zero sotto il livello di elev_fin (si può discutere...), potrei usare first_level_static confrontare e in caso sia sotto porre cl=1
                if (l < anaprop.elev_fin(i, k)) {
                    qual.scan(l).set(i, k, 0);
                    cl=2;
                } else {
                    //--bisogna ragionare di nuovo su definizione di qualità con clutter se si copia il dato sopra.--

                    //--calcolo la qualità--
                    // FIXME: qui tronca: meglio un round?
                    qual.scan(l).set(i, k, (unsigned char)(func_q_Z(cl,bb,dist,drrs,hray_inf.dtrs,dh,dhst,PIA)*100));
                }

                if (qual.scan(l).get(i, k) ==0) qual.scan(l).set(i, k, 1);//????a che serve???
                if (calcolo_vpr)
                {
                    /* sezione PREPARAZIONE DATI VPR*/
                    if(cl==0 && bb<BBMAX_VPR )   /*pongo le condizioni per individuare l'area visibile per calcolo VPR, riduco il bb ammesso (BBMAX_VPR=20)*/ //riveder.....?????
                        flag_vpr.scan(l).set(i, k, 1);
                }
            }
        }
    }

    LOG_DEBUG("End caratterizzo_volume");
    return;
}

void CalcoloVPR::classifica_rain()
{
    LOG_CATEGORY("radar.class");
    int hmax=-9999;

    /* ;---------------------------------- */
    /* ;          FASE 0 :                  */
    /* ;---------------------------------- */
    // DEFINISCO QUOTE DELLA BASE E DEL TOP DELLA BRIGHT BAND USANDO IL DATO quota del picco  DEL PRECEDENTE RUN O, SE NON PRESENTE LA QUOTA DELLO ZERO DA MODELLO

    // Lettura quota massimo da VPR  calcolo base e top bright band
    LOG_INFO("data= %s",cum_bac.date);
    // calcolo il gap
    gap = cum_bac.assets.read_profile_gap();
    //-- se gap < memory leggo hmax da VPR
    if (gap<=MEMORY){
        hmax = cum_bac.assets.read_vpr_hmax();
            //---suppongo una semiampiezza massima della bright band di 600 m e definisco htopbb e hbasebb come hmassimo +600 m (che da clima ci sta) e hmassimo -600 m
    }

    if (hmax >= 0)
    {
        hbbb=(hmax-600.)/1000.;
        htbb=(hmax+600.)/1000.;
    } else {
        //-- se gap > memory o se non ho trovato il file
        // Lettura 0 termico da modello, e calcolo base e top bright band
        LOG_INFO("leggo 0termico per class da file %s",getenv("FILE_ZERO_TERMICO"));
        //  leggo informazioni di temperatura da modello*/
        float zeroterm;//zerotermico
        if (cum_bac.assets.read_0term(zeroterm))
        {
            //-- considerato che lo shift medio tra il picco e lo zero è tra 200 e 300 m, che il modello può avere un errore, definisco cautelativamente htbb come quota zero + 400 m e hbbb come  quota zero -700 m  .
            htbb=zeroterm/1000. + 0.4; // se non ho trovato il vpr allora uso un range più ristretto, potrebbe essere caso convettivo
            hbbb=zeroterm/1000. - 1.0;
        } else {
            LOG_ERROR("non ho trovato il file dello zero termico");
            LOG_INFO("attenzione, non ho trovat zero termico ne da vpr ne da radiosondaggio");
            htbb=0.; // discutibile così faccio tutto con VIZ
            hbbb=0.;
        }
    }

    // se hbasebb è <0 metto 0
    if (hbbb<0.) hbbb=0.;

    LOG_INFO("calcolati livelli sopra e sotto bright band hbbb=%f  htbb=%f",hbbb,htbb);

    const CylindricalVolume& cil = cum_bac.cil;

    // ricampionamento del volume in coordinate cilindriche
    LOG_DEBUG ("Matrice cilindrica Naz %3d Nrange %4d Nheight %4d", cil.slices.size(), cil.x_size, cil.z_size);
    //-------------------------------------------------------------------------------------------------------------------------
    // faccio la classificazione col metodo Vertical Integrated Reflectivity
    algo::CalcoloVIZ viz(cil, htbb, hbbb, t_ground);
    viz.classifico_VIZ();

    //classificazione con STEINER
    //  if (hmax > 2000.) {// per evitare contaminazioni della bright band, si puo' tunare
    // if (hbbb > 500.) {// per evitare contaminazioni della bright band, si puo' tunare

    algo::CalcoloSteiner steiner(cum_bac.volume, cum_bac.anaprop.elev_fin, cil.x_size);
    steiner.calcolo_background();
    steiner.classifico_STEINER();
    //  }
    merge_metodi(steiner, viz);
    return ;
}

void CalcoloVPR::merge_metodi(const algo::CalcoloSteiner& steiner, const algo::CalcoloVIZ& viz)
{
    for (unsigned j = 0; j < NUM_AZ_X_PPI; ++j)
        for (unsigned k = 0; k < steiner.max_bin; ++k)
            if (cum_bac.anaprop.quota(j, k) < hbbb*1000. && steiner.conv_STEINER(j, k) == viz.conv_VIZ(j, k) && steiner.conv_STEINER(j, k) > 0)
                conv(j,k) = steiner.conv_STEINER(j, k);
            else
                if (steiner.conv_STEINER(j, k) == viz.conv_VIZ(j, k) && steiner.conv_STEINER(j, k) > 0 && viz.stratiform(j, k) < 1)
                    conv(j,k) = viz.conv_VIZ(j, k);
}

//----------ALGORITMO
/*  combina il profilo verticale corrente con quello precedente tramite il metodo di Germann (2003)
    a) calcolo gap tra ultimo profilo e istante corrente
    b) se gap < MEMORY leggo profilo
    c) se gap> MEMORY peso 0 il profilo storico e cerco oltre la data (per casi vecchi)
    d) faccio func_vpr
    f) cerco il profilo con cui combinare (->proprio, se gap<MEMORY ->dell'altro radar se gap_res<MEMORY e profile_heating_res=WARM)
    g) Combino livelli con peso sottostante
    Dati cv e ct, volume totale e volume precipitante il peso del vpr istantaneo è calcolato come segue:
    c0=2*cv;
    peso=(float)ct/(c0+ct)
    long int c0,cv,ct; costanti di combinazione (v. ref.)
    h) trovo livello minimo, se livello minimo profilo combinato più alto del precedente calcolo la diff media e sommo al vecchio
    e) ricalcolo livello minimo
    float vpr0.val[NMAXLAYER],vpr1.val[NMAXLAYER],vpr.val[NMAXLAYER]; profilo precedente, ultimo e combinato
    float alfat,noval; peso, nodata
    FILE *file;
    int mode,ilay;  modalità calcolo profilo (0=combinazione, 1=istantaneo),indice di strato
*/
int CalcoloVPR::combina_profili(const InstantaneousVPR& inst_vpr)
{
    LOG_CATEGORY("radar.vpr");

    LOG_DEBUG (" modalita %d", MOD_VPR);
    VPR vpr0;
    bool combinante; // combinante: variabile che contiene presenza vpr alternativo
    if (MOD_VPR == 0)
    {
        /* MOD_VPR=0: VPR combinato */
        combinante = cum_bac.assets.find_vpr0(cum_bac.dbz, vpr0, gap);

        for (unsigned i=0; i<vpr0.size(); i++)
            LOG_DEBUG (" Profilo vecchio - livello %2d valore %6.2f",i,vpr0.val[i]);
        //----a fine calcolo sul sito in esame stampo il valore del gap
        LOG_INFO("gap %li",gap);
    } else {
        /* MOD_VPR=1: VPR istantaneo */
        combinante = false;
    }

    if (combinante)
    {
        if (inst_vpr.success)
        {
            vpr = combine_profiles(vpr0, inst_vpr.vpr, inst_vpr.cv, inst_vpr.ct);
        } else {
            // se il calcolo dell'istantaneo non è andato bene , ricopio l'altro vpr e la sua area
            vpr = vpr0;
        }
    } else {
        if (inst_vpr.success)
        {
            // se il calcolo dell'istantaneo  è andato bene ricopio il profilo
            vpr = inst_vpr.vpr;
        } else {
            //-----se è andata male la ricerca dell'altro e anche il calcolo dell'istantaneo esco
            return 1;
        }
    }

    //------------- trovo livello minimo -------
    Livmin livmin(vpr);
    LOG_INFO(" livmin %i", livmin.livmin);

    if (livmin.idx >= vpr.size() - 1 || !livmin.found)
        return (1);

    this->livmin = livmin.livmin;


    //-----scrivo il profilo e la sua area-----
    cum_bac.assets.write_vpr0(vpr);
    for (unsigned i=0; i<vpr.size(); i++) LOG_DEBUG (" Profilo nuovo - livello %2d valore %6.2f",i,vpr.val[i]);

    return(0);
}

int CalcoloVPR::profile_heating(bool has_inst_vpr)
#include <radarelab/vpr_par.h>
{
    LOG_CATEGORY("radar.vpr");
    //---leggo ultimo file contenente riscaldamento , se non esiste impongo heating=0 (verificare comando)
    int heating = cum_bac.assets.read_vpr_heating();

    //--una volta letto il file, se il calcolo del vpr è andato bene incremento di uno heating sottraendo però la differenza di date (in quarti d'ora)-1 tra gli ultimi due profili
    //--lo faccio perchè potrei avere heating più alto del dovuto se ho avuto un interruzione del flusso dei dati
    //-- heating ha un valore massimo pari  WARM dopodichè diventa heating = MEMORY e così resta finchè non sono passati MEMORY istanti non aggiornati ( va bene?)
    //---se il profil non  è stato aggiornato invece decremento la variabile riscaldamento di gap con un minimo pari a 0
    if (!has_inst_vpr){
        heating=heating-gap; /*se il profilo non è stato aggiornato, ho raffreddamento, in caso arrivi sotto WARM riparto da 0, cioè serve riscaldamento  */
    }
    else  {
        heating = heating - (gap-1) + 1; /*se il profilo è stato aggiornato, ho riscaldamento , in caso arrivi sopra WARM riparto da MEMORY  */
        if (heating>=WARM) heating=MEMORY;  /* se heating raggiunge WARM allora lo pongo uguale a MEMORY     */
    }
    if (heating<0) heating=0;

    //----stampo heating su file
    cum_bac.assets.write_vpr_heating(heating);

    //----stampo log vpr
    LOG_INFO("gap %li heating %i",gap,heating);

    return(heating);
}


int CalcoloVPR::stampa_vpr()
{
    float vpr_dbz;

    File file(logging_category);
    file.open_from_env("VPR_ARCH", "wt", "ultimo vpr in dBZ per il plot");

    fprintf(file," QUOTA   DBZ    AREA PRECI(KM^2/1000)\n" );
    for (unsigned ilay=0;  ilay<NMAXLAYER; ilay++){
        if (vpr.val[ilay]> 0.001 ) {
            vpr_dbz=cum_bac.dbz.RtoDBZ(vpr.val[ilay]);
            fprintf(file," %i %10.3f %li\n", ilay*TCK_VPR+TCK_VPR/2, vpr_dbz, vpr.area[ilay]);
        }
        else
            fprintf(file," %i %10.3f %li\n", ilay*TCK_VPR+TCK_VPR/2, NODATAVPR, vpr.area[ilay]);
    }

    return 0;
}
/*=======================================================================================*/
/*
   comstart corr_vpr
   idx corregge i dati tramite profilo verticale
   corregge i dati tramite profilo verticale a partire da quelli con valore maggiore di THR_CORR(v vpr_par.h): riporto il dato alla quota del livello 'liquido' tramite "traslazione" cioè aggiungo al valore in dBZ la differenza tra il valore del VPR alla quota 'liquida' e quello alla quota della misura, anche in caso di neve,purchè esista il livello liquido nel profilo. In questo caso pero' flaggo il bin tramte la variabile neve[][]. In caso il profilo intero sia di neve allora riporto al valore di Z al suolo (o al livello rappresentativo) perchè non ho un valore di riferimento 'liquido'.
   La correzione avviene solo se heating>warm

   int ilref : livello del suolo o della quota rappresentativa di esso nel VPR (dove considero buoni i dati a partire dati da REPR_LEV)
   int ilray : livello a cui ho il dato
   int ilay2 : livello sopra o sotto ilray a seconda che il fascio stia sopra o sotto il centro di ilray, serve per interpolare il dato vpr su retta ed evitare correzioni a gradino
   int heating,warm: indicano quanto è caldo il profilo, e la soglia di riscaldamento
   int snow : indica se il profilo è di neve (snow=1)
   int neve[NUM_AZ_X_PPI][MAX_BIN]: 1 se c'è neve, identificata se quota dem> hvprmax+300m
   float corr: correzione in dB
   float vpr_liq: valore di R nel VPR alla quota 'liquida' ricavato dalla funzione analyse_VPR

   ilref= (dem[i][k]>REPR_LEV)?(floor(dem[i][k]/TCK_VPR)):( floor(REPR_LEV/TCK_VPR)); in pratica livello dem se > REPR_LEV, livello di REPR_LEV altrimenti.
   ilray=floor((hbin[i][k])/TCK_VPR)
   corr=RtoDBZ(vpr_liq)-RtoDBZ(vpr_hray)
   comend
*/
int CalcoloVPR::corr_vpr()
    //* ====correzione profilo====================================*/

#include <radarelab/vpr_par.h>

{
    LOG_CATEGORY("radar.vpr");

    int ilray,ilref,ilay2,ier_ana,snow,strat;
    float corr,vpr_liq,vpr_hray,hbin,hliq;

    /*inizializzazione variabili */
    snow=0;
    //vpr al livello liquido liv liquido e liv max
    vpr_liq=NODATAVPR;
    hliq=NODATAVPR;
    hvprmax=INODATA;

    // analisi vpr

    ier_max=trovo_hvprmax(&hvprmax);
    ier_ana=analyse_VPR(&vpr_liq,&snow,&hliq);
    LOG_INFO("ier_analisi %i",ier_ana) ;

    /* se analisi dice che non è il caso di correggere non correggo (NB in questo caso non riempio la matrice di neve)*/
    if (ier_ana) return 1;

    LOG_INFO("altezza bright band %i",hvprmax);
    LOG_INFO("CORREGGO VPR");

    //correzione vpr
    for (unsigned i=0; i<NUM_AZ_X_PPI; i++){
        for (unsigned k=0; k<cum_bac.volume[0].beam_size; k++){
            corr=0.;
            /* trovo elevazione reale e quota bin*/
            //elevaz=(float)(volume_at_elev_preci(i, k).teta_true)*CONV_RAD;
            hbin=(float)cum_bac.anaprop.quota(i, k);

            /* se dall'analisi risulta che nevica assegno neve ovunque*/
            if (snow) neve(i, k)=1;
            strat=1;
            if (cum_bac.do_class)
            {
                if (conv(i,k) >= CONV_VAL){
                    strat=0;
                }
            }
            //--- impongo una soglia per la correzione pari a 0 dBZ
            if (cum_bac.volume[0].get(i, k) > THR_CORR && hbin > hliq  && strat){

                //---trovo lo strato del pixel, se maggiore o uguale a NMAXLAYER lo retrocedo di 2, se minore di livmn lo pongo uguale a livmin
                ilray=(hbin>=livmin)?(floor(hbin/TCK_VPR)):(floor(livmin/TCK_VPR));//discutibile :livello del fascio se minore di livmin posto=livmin
                if (ilray>= NMAXLAYER ) ilray=NMAXLAYER-2;//livello del fascio se >= NMAXLAYER posto =NMAXLAYER-2

                //---trovo ilay2 strato con cui mediare per calcolare il vpr a una quota intermedia tra 2 livelli, se l'altezza del bin è sopra metà strato prendo quello sopra altrimenti quello sotto
                if ((int)hbin%TCK_VPR > TCK_VPR/2) ilay2=ilray+1;
                else ilay2=ilray-1;
                if (ilay2< floor(livmin/TCK_VPR)) ilay2=floor(livmin/TCK_VPR);

                //trovo ilref: livello di riferimento per ricostruire il valore vpr al suolo nel caso di neve.
                // in caso di profilo di pioggia mi riporto sempre al valore del livello liquido e questo può essere un punto critico.. vedere come modificarlo.

                ilref=(cum_bac.dem(i, k)>livmin)?(floor(cum_bac.dem(i, k)/TCK_VPR)):(floor(livmin/TCK_VPR));//livello di riferimento; se livello dem>livmin = livello dem altrimenti livmin


                if (vpr.val[ilref] > 0 && vpr.val[ilray] > 0 ){ /*devo avere dati validi nel VPR alle quote considerate!*/
                    //-- calcolo il valore del profilo alla quota di interesse
                    vpr_hray=vpr.val[ilray]+((vpr.val[ilray]-vpr.val[ilay2])/(ilray*TCK_VPR-TCK_VPR/2-ilay2*TCK_VPR))*(hbin-ilray*TCK_VPR-TCK_VPR/2); /*per rendere la correzione continua non a gradini */
                    //--identifico le aree dove nevica stando alla quota teorica dello zero termico

                    if (cum_bac.dem(i, k)> hvprmax+HALF_BB-TCK_VPR || snow){ /*classifico neve*/
                        neve(i, k)=1;

                    }

                    //--se nevica la correzione consiste solo nel riportare il valore del vpr al suolo: PROPOSTA: qui si potrebbe generare una mappa di intensità di neve ma deve essere rivisto tutto


                    //if(snow) //A rimosso, faccio una cosa diversa
                    if(neve(i, k)){

                        //faccio la regressione lineare dei punti del profilo sopra il punto del dem
                        //calcolo il valore al livello del dem e lo sostituisco a vpr.val[ilref] nella correzione
                        // faccio linearizzazione in maniera becera:
                        //vpr.val[ilref]=(vpr.val[ilref+7]-vpr.val[ilref+2])/(5)*(ilref-(ilref+2))+vpr.val[ilref+2];

                        //passaggio=BYTEtoR(volume.vol_pol,aMP_SNOW,bMP_SNOW)

                        //volpol[0][i][k]=RtoBYTE(passaggio)

                        corr=cum_bac.dbz.RtoDBZ(vpr.val[ilref])-cum_bac.dbz.RtoDBZ(vpr_hray);

                        cum_bac.volume[0].set(i, k, cum_bac.dbz.DBZ_snow(cum_bac.volume[0].get(i, k)));
                    }
                    else{
                        // -- altrimenti correggo comunque a livello liquido :
                        corr = cum_bac.dbz.RtoDBZ_class(vpr_liq) - cum_bac.dbz.RtoDBZ_class(vpr_hray);/*riporto comunque al valore liquido anche se sono sopra la bright band*/
                    }
                    // --  controllo qualità su valore correzione
                    if (corr>MAX_CORR) corr=MAX_CORR; /*soglia sulla massima correzione*/
                    if (hbin<hvprmax && corr>0.) corr=0; /*evito effetti incrementi non giustificati*/

                    //controllo qualità su valore corretto e correzione
                    double corrected = cum_bac.volume[0].get(i, k) + corr;
                    if (corrected > MAXVAL_DB) // se dato corretto va fuori scala assegno valore massimo
                        cum_bac.volume[0].set(i, k, MAXVAL_DB);
                    else if ( corrected < MINVAL_DB) // se dato corretto va a fodoscala assegno valore di fondo scala
                        cum_bac.volume[0].set(i, k, MINVAL_DB);
                    else
                        cum_bac.volume[0].set(i, k, corrected);  // correggo

                    corr_polar(i, k)=(unsigned char)(corr)+128;

                    //inserisco un ponghino per rifare la neve con aMP e bMP modificati // DA SCOMMENTARE SE DECIDO DI FARLO

                    //if (neve[i][k]) volume.scan(0).get_raw(i, k)=DBtoBYTE(RtoDBZ( BYTE_to_mp_func(volume.scan(0).get_raw(i, k),aMP_SNOW,bMP_SNOW),aMP_class,bMP_class )) ;


                }
            }
        }
    }
    return(0);
}

int CalcoloVPR::trovo_hvprmax(int *hmax)
{
    int i,imax,istart,foundlivmax;
    float h0start,peak,soglia;


    if (t_ground != NODATAVPR)
    {
        LOG_DEBUG("trovo hvprmax  a partire da 400 m sotto lo zero dell'adiabatica secca");
        h0start=t_ground/9.8*1000 ;
        istart=h0start/TCK_VPR -2;
        if (istart< livmin/TCK_VPR) istart=livmin/TCK_VPR;
        LOG_DEBUG("t_ground h0start istart %f %f %i",t_ground,h0start,istart);
    }
    else {
        LOG_DEBUG("trovo hvprmax  a partire da livmin");
        istart=livmin/TCK_VPR+1;
    }


    /* trovo hvprmax e il suo livello a partire dal livello istart */

    //--inizializzazione
    foundlivmax=0;
    peak=NODATAVPR;
    *hmax=INODATA;
    // Enrico vprmax=NODATAVPR;
    imax=INODATA;
    soglia = DBZ::DBZtoR(THR_VPR,200,1.6); // CAMBIATO, ERRORE, PRIMA ERA RtoDBZ!!!!VERIFICARE CHE IL NUMERO PARAMETRI FUNZIONE SIA CORRETTO

    //--se vpr al livello corrente e 4 layer sopra> soglia, calcolo picco
        LOG_DEBUG(" istart %d low %6.2f  up %6.2f  soglia %6.2f  peak %6.2f  imax %d", istart, vpr.val[istart] , vpr.val[istart+4], soglia, peak, imax); 
    if (vpr.val[istart] >soglia && vpr.val[istart+4] > soglia){
        peak=10*log10(vpr.val[istart]/vpr.val[istart+4]);//inizializzo il picco
        LOG_DEBUG("peak1 = %f",peak);
    }
    //----se picco > MINIMO il punto è ok
    if(peak> MIN_PEAK_VPR){
        imax=istart;
        // Enrico vprmax=vpr.val[imax];
        LOG_DEBUG("il primo punto soddisfa le condizioni di picco");
    }
    for (i=istart+1;i<NMAXLAYER-4;i++) //la ricerca è un po' diversa dall'originale.. trovo il picco + alto con valore  rispetto a 4 sopra > soglia
    {
        if (vpr.val[i] <soglia || vpr.val[i+4] < soglia) break;
        peak=10*log10(vpr.val[i]/vpr.val[i+4]);
        if (vpr.val[i]>vpr.val[i-1]  && peak> MIN_PEAK_VPR ) // se vpr(i) maggiore del massimo e picco sufficientemente alto
        {
            imax=i;
            // Enrico vprmax=vpr.val[imax];
        }
        LOG_DEBUG(" low %6.2f  up %6.2f  soglia %6.2f  peak %6.2f  imax %d", vpr.val[i] , vpr.val[i+4], soglia, peak, imax); 
    }

    if ( imax  > INODATA ){
        foundlivmax=1;
        peak=10*log10(vpr.val[imax]/vpr.val[imax+4]);
        *hmax=imax*TCK_VPR+TCK_VPR/2;
        LOG_DEBUG("trovato ilaymax %i %i",*hmax,imax);
        LOG_DEBUG(" picco in dbR %f",peak);
    }

    LOG_DEBUG("exit status trovo_hvprmax %i",foundlivmax);
    return (foundlivmax);
}

/*
   1)hvprmax - semiampiezza BB > liv 100 m => la bright band sta sopra il suolo => interpolo il profilo per trovare il livello liquido
   se interpolazione fallisce NON CORREGGO (scelta cautelativa, potrei avere caso convettivo
   o pochi punti o molto disomogeneo)
   2)hvprmax - semiampiezza BB < liv 100 m  => la bright contiene o è sotto il livello 100 m  oppure ho un massimo 'spurio'
   quindi calcolo la Tground, temperatura media nelle stazioni più vicine al radar, se non trovo T torno al caso 1.
   2A) Tground>T_MAX_ML ----> prob. caso convettivo o max preci passa sopra radar, interpolo il profilo e calcolo livello liquido
   se interpolazione fallisce NON CORREGGO
   2B) T_MIN_ML<T0<T_MAX_ML : prob. bright band al suolo, interpolo il profilo per trovare il livello liquido
   se interpolazione fallisce NON CORREGGO
   2C) Tground<T_MIN_ML ----> prob. profilo neve NON interpolo e non calcolo vpr al livello liquido perchè non ho livello liquido

   comend
*/
int CalcoloVPR::analyse_VPR(float *vpr_liq,int *snow,float *hliq)
    /*=======analisi profilo============ */
{
    int ier=1,ier_ana=0,liv0;
    char date[]="000000000000";
    struct tm *tempo;
    time_t Time;

    // ------------inizializzazione delle variabili ----------

    //strcpy(date,"000000000000");

    int tipo_profilo=-1;
    float v600sottobb=NODATAVPR;
    float v1000=NODATAVPR;
    float v1500=NODATAVPR;
    float vliq=NODATAVPR;
    float vhliquid=NODATAVPR;
    float vprmax=NODATAVPR;
    //*togliere gli ultimi tre*/;

    //ier_max=trovo_hvprmax(&hvprmax);


    if (t_ground == NODATAVPR) //1  se non ho nè T nè il massimo esco altrimenti tipo_profilo=0
    {
        LOG_WARN("non ho T,...");

        if ( ! ier_max ) {
            LOG_ERROR(" non ho trovato hvprmax, nè T, esco");
            return 1;
        }
        tipo_profilo=0;
    }
    else
    {

        if (t_ground >= T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100.){ //1  se T > T_MAX_ML e non ho il massimo esco
            if ( ! ier_max ) {
                LOG_ERROR(" temperatura alta e non ho trovato hvprmax, esco");
                return 1;
            }
            tipo_profilo=0;
        }


        // if (t_ground >= T_MAX_SN+0.65*(float)(livmin+TCK_VPR/2)/100  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100. )
        if (t_ground >= T_MAX_SN  && t_ground < T_MAX_ML+0.65*(float)(livmin+TCK_VPR/2)/100. )
        {

            if (  ier_max ) {
                LOG_INFO(" temperatura da scioglimento e massimo in quota");
                tipo_profilo=2;
            }
            else{
                LOG_ERROR(" temperatura da scioglimento ma superiore a temperatura max neve e non ho trovato hvprmax, esco");
                return 1;
            }
            // solo una scritta per descrivere cos'è accaduto
            liv0=livmin+HALF_BB;
            if (hvprmax > liv0) LOG_INFO(" il livello %i è sotto la Bright band, ma T bassa  interpolo",livmin);
            else LOG_INFO(" il livello %i potrebbe essere dentro la Bright Band, interpolo",livmin);

        }

        //if (t_ground >= T_MIN_ML  && t_ground < T_MAX_SN+0.65*(float)(livmin+TCK_VPR/2)/100.)
        if (t_ground < T_MAX_SN)
        {
            if ( ier_max ){
                LOG_INFO(" temperatura da neve o scioglimento e massimo in quota");
                tipo_profilo=2;
            }
            else {
                LOG_INFO(" temperatura da neve o scioglimento e massimo non trovato,neve , non interpolo");
                tipo_profilo=3;
                hvprmax=0;
            }
        }

    }

    // InterpolaVPR_NR iv;
    InterpolaVPR_GSL iv;

    switch
        (tipo_profilo)
        {
            case 0:
            case 1:
            case 2:
                ier=iv.interpola_VPR(vpr.val.data(), hvprmax, livmin);
                if (ier){
                    LOG_INFO(" interpolazione fallita");
                    switch (tipo_profilo)
                    {

                        /*Questo fallisce se la bright band non è al suolo (per evitare correzioni dannose in casi poco omogenei)*/
                        case 0:
                        case 1:
                            ier_ana=1;
                            *vpr_liq=NODATAVPR;
                            *hliq=NODATAVPR;
                            break;
                        case 2:
                            *vpr_liq=vpr.val[(hvprmax+1000)/TCK_VPR]*2.15;/*21 aprile 2008*/
                            *hliq=0;
                            break;
                    }
                }
                else{
                    LOG_INFO(" interpolazione eseguita con successo");
                    //
                    // stampa del profilo interpolato
                    const char* vpr_arch = getenv("VPR_ARCH");
                    if (!vpr_arch) throw runtime_error("VPR_ARCH is not defined");
                    string fname(vpr_arch);
                    fname += "_int";
                    File file(logging_category);
                    file.open(fname, "wt", "vpr interpolato");
                    for (unsigned i = 0; i < NMAXLAYER; ++i)
                        fprintf(file," %f \n", cum_bac.dbz.RtoDBZ(iv.vpr_int[i]));

                    /*calcolo valore di riferimento di vpr_liq per l'acqua liquida nell'ipotesi che a[2]=quota_bright_band e a[2]-1.5*a[3]=quota acqua liquida*/
                    if (tipo_profilo == 2 ) {
                        *hliq=(iv.E-2.1*iv.G)*1000.;
                        //lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);
                        if (*hliq<0)
                            *hliq=0;  /*con casi di bright band bassa.. cerco di correggere il più possibile*/
                        *vpr_liq=vpr.val[(hvprmax+1000)/TCK_VPR]*2.15;
                    }
                    else {
                        *hliq=(iv.E-2.1*iv.G)*1000.;
                        //lineargauss(a[2]-2.1*a[3], a, vpr_liq, dyda, ndata);
                        if ( *hliq > livmin) {
                            *vpr_liq=vpr.val[(int)(*hliq/TCK_VPR)]; // ... SE HO IL VALORE VPR USO QUELLO.
                        }
                        else // altrimenti tengo il valore vpr neve + 6 dB* e metto tipo_profilo=2
                        {
                            if (*hliq<0) *hliq=0;
                            tipo_profilo=2;
                            //*vpr_liq=vpr.val[(hvprmax+1000)/TCK_VPR]*2.15;
			    // iv.C = risultato dell'interpolazione (interpolated vpr)
                            *vpr_liq=iv.C;
                        }
                    }
                }
                break;
            case 3:
                *snow=1;
                *vpr_liq=NODATAVPR;
                *hliq=NODATAVPR;
                break;
        }
    LOG_INFO("TIPO_PROFILO= %i vpr_liq %f hliq %f", tipo_profilo, *vpr_liq,*hliq );


    /* parte di stampa test vpr*/

    /* nome data */
    //definisco stringa data in modo predefinito
    Time = cum_bac.volume.load_info->acq_date;
    tempo = gmtime(&Time);
    sprintf(date,"%04d%02d%02d%02d%02d",tempo->tm_year+1900, tempo->tm_mon+1,
            tempo->tm_mday,tempo->tm_hour, tempo->tm_min);
    if (! ier ) {
        if(*hliq > livmin +200 )
            vhliquid=cum_bac.dbz.RtoDBZ(vpr.val[(int)(*hliq)/TCK_VPR]);
        vliq=cum_bac.dbz.RtoDBZ(*vpr_liq);
    }
    if (ier_max) {
        if ( hvprmax-600 >= livmin )
            v600sottobb=cum_bac.dbz.RtoDBZ(vpr.val[(hvprmax-600)/TCK_VPR]);
        if ((hvprmax+1000)/TCK_VPR < NMAXLAYER )
            v1000=cum_bac.dbz.RtoDBZ(vpr.val[(hvprmax+1000)/TCK_VPR]);
        if ((hvprmax+1500)/TCK_VPR < NMAXLAYER )
            v1500=cum_bac.dbz.RtoDBZ(vpr.val[(hvprmax+1500)/TCK_VPR]);
        vprmax=cum_bac.dbz.RtoDBZ(vpr.val[(hvprmax/TCK_VPR)]);
    }

    if (FILE* test_vpr=fopen(getenv("TEST_VPR"),"a+"))
    {
        fprintf(test_vpr,"%s %i %i -1 %f %f  %f %f %f %f %f %f %f %f  %f %f %f  %f \n",date,hvprmax,tipo_profilo,iv.chisqfin,*hliq,vliq,vhliquid,v600sottobb,v1000+6,v1500+6,vprmax,iv.rmsefin,iv.B,iv.E,iv.G,iv.C,iv.F);
        fclose(test_vpr);
    }

    // fine parte di stampa test vpr

    //---SCRIVO ALTEZZA MASSIMO PER CLASSIFICAZIONE AL RUN SUCCESSIVO

    cum_bac.assets.write_vpr_hmax(hvprmax);
    LOG_INFO("fatta scrittura hmax vpr = %d",hvprmax);

    return (ier_ana);
}

void CUM_BAC::conversione_convettiva()
{
    for (unsigned i=0; i<NUM_AZ_X_PPI; i++)
        for (unsigned k=0; k<volume[0].beam_size; k++)
            if (calcolo_vpr->conv(i,k) > 0)
                volume[0].set(i, k, dbz.DBZ_conv(volume[0].get(i, k)));
}

void CUM_BAC::vpr_class()
{
    //--------------se definita la qualita procedo con il calcolo qualita e del VPR (perchè prendo solo i punti con qual > soglia?)-----------------//
    if (do_quality)
    {
        //-------------calcolo qualita' e trovo il top

        caratterizzo_volume();
        /* //---------trovo il top (a X dbZ) */
        /* printf ("trovo top \n") ; */
        /* ier=trovo_top(); */


        //--------------se definita CLASS procedo con  classificazione -----------------//
        if (do_class)
	 {
            calcolo_vpr->classifica_rain();
	  }

        //--------------se definito VPR procedo con calcolo VPR -----------------//

        if (calcolo_vpr)
            calcolo_vpr->esegui_tutto();
    }

    if (do_class)
        conversione_convettiva();
}

void CUM_BAC::generate_maps(CartProducts& products)
{
    // Generate products and write them out
    LOG_INFO("Scrittura File Precipitazione 1X1\n");
    if (do_zlr_media)
    {
        std::function<unsigned char(const vector<double>&)> convert = [this](const vector<double>& samples) {
            // Samples are in contains dB (logaritmic values), so mediating
            // them is not the correct operation, and we need to convert
            // them to Z (linear) values to average them.
            // TODO: there may be more efficient way to mediate logaritmic
            // values.
            double sum = 0;
            for (const auto& s: samples)
                sum += DBZ::DBZtoZ(s);
            unsigned char res = DBZ::DBtoBYTE(DBZ::ZtoDBZ(sum / samples.size()));
            // il max serve perchè il valore di MISSING è 0
            if (res == 0) return (unsigned char)1;
            return res;
        };
        products.scaled.to_cart_average(volume[0], convert, products.z_out);
    } else {
        std::function<unsigned char(unsigned, unsigned)> assign_cart =
            [this](unsigned azimuth, unsigned range) {
                // il max serve perchè il valore di MISSING è 0
                unsigned char sample = DBZ::DBtoBYTE(volume[0].get(azimuth, range));
                return max(sample, (unsigned char)1);
            };
        products.scaled.to_cart(assign_cart, products.z_out);
    }

    std::function<unsigned char(unsigned, unsigned)> assign_cart =
         [this](unsigned azimuth, unsigned range) {
         // il max serve perchè il valore di MISSING è 0
         unsigned char sample = DBZ::DBtoBYTE(volume[0].get(azimuth, range));
         return max(sample, (unsigned char)1);
    };
    products.fullres.to_cart(assign_cart, products.z_fr);

    products.scaled.to_cart(top, products.top_1x1);
    if (do_quality)
    {
        const auto& elev_fin = anaprop.elev_fin;
        const auto& quota = anaprop.quota;

        std::function<unsigned char(unsigned, unsigned)> assign_qual =
          [this, &elev_fin](unsigned azimuth, unsigned range) {
              const auto& el = elev_fin(azimuth, range);
              if (range >= volume[el].beam_size)
                  return (unsigned char)0;
              return qual.scan(el).get(azimuth, range);
          };
        products.scaled.to_cart(assign_qual, products.qual_Z_1x1);
        
        std::function<unsigned char(unsigned, unsigned)> assign_quota =
            [&quota](unsigned azimuth, unsigned range) {
                return 128 + round(quota(azimuth, range) / 100.0);
            };
        products.scaled.to_cart(assign_quota, products.quota_1x1);
        products.scaled.to_cart(anaprop.dato_corrotto, products.dato_corr_1x1);

        std::function<unsigned char(unsigned, unsigned)> assign_elev_fin = [&elev_fin](unsigned azimuth, unsigned range) {
                return elev_fin(azimuth, range);
        };
        products.scaled.to_cart(assign_elev_fin, products.elev_fin_1x1);
        
        products.scaled.to_cart(beam_blocking, products.beam_blocking_1x1);
    }

    if (calcolo_vpr)
    {
        const auto& neve = calcolo_vpr->neve;
        std::function<unsigned char(unsigned, unsigned)> assign = [&neve](unsigned azimuth, unsigned range) {
            return neve(azimuth, range) ? 0 : 1;
        };
        products.scaled.to_cart(assign, products.neve_1x1);
        
        products.scaled.to_cart(calcolo_vpr->corr_polar, products.corr_1x1);

        if (do_class)
            products.scaled.to_cart(calcolo_vpr->conv, products.conv_1x1);
    }

    if (do_devel)
    {
        SD_Z6 *= 10.;

        SingleCart SC_SD(SD_Z6.max_beam_size());
        for (unsigned int i=0; i<SD_Z6.size(); i++){
           SC_SD.creo_cart(SD_Z6, i);
           std::ostringstream oss;
           oss<<"SD_"<<i;
           SC_SD.write_out(assets,oss.str());
        }
    }

//    return true;
}


CalcoloVPR::CalcoloVPR(CUM_BAC& cum_bac)
    : cum_bac(cum_bac),
      inst_vpr(cum_bac.volume, cum_bac.qual, cum_bac.flag_vpr, cum_bac.site.vpr_iaz_min, cum_bac.site.vpr_iaz_max),
      conv(NUM_AZ_X_PPI, cum_bac.volume.max_beam_size(),0),
      corr_polar(NUM_AZ_X_PPI, cum_bac.volume.max_beam_size()),
      neve(NUM_AZ_X_PPI, cum_bac.volume.max_beam_size())
{
    logging_category = log4c_category_get("radar.vpr");
    MyMAX_BIN=cum_bac.MyMAX_BIN;
    htbb=-9999.; hbbb=-9999.;
    t_ground=NODATAVPR;

    /*
    for (int k=0; k<NUM_AZ_X_PPI*MyMAX_BIN;k++ ){
      lista_conv[k][0]=-999;
      lista_conv[k][1]=-999;
    }
    */

    t_ground = cum_bac.assets.read_t_ground();
}

CalcoloVPR::~CalcoloVPR()
{
}

void CalcoloVPR::esegui_tutto()
{
    LOG_INFO("processo file dati: %s", cum_bac.volume.load_info->filename.c_str());
    printf ("calcolo VPR \n") ;

    //VPR  // ------------inizializzo hvprmax ---------------

    hvprmax=INODATA;

    //VPR  // ------------chiamo combina profili con parametri sito, sito alternativo ---------------

    inst_vpr.compute(); // ho fatto func_vpr, il profilo istantaneo
    LOG_INFO("fatta func vpr %s", inst_vpr.success ? "ok" : "errore");
    for (unsigned i=0; i<inst_vpr.vpr.size(); i++) LOG_DEBUG (" Profilo istantaneo - livello %2d valore %6.2f",i,inst_vpr.vpr.val[i]);

    int ier_comb;               ///< flag d'errore su combinazione vpr

    //  ier_comb=combina_profili(sito,argv[4]);
    ier_comb=combina_profili(inst_vpr);
    LOG_INFO ("exit status calcolo VPR istantaneo: %s", inst_vpr.success ? "ok": "fallito") ; // debug
    LOG_INFO("exit status combinaprofili: (1--fallito 0--ok) %i ",ier_comb) ; // debug


    //VPR  // ------------chiamo profile_heating che calcola riscaldamento profilo ---------------

    heating=profile_heating(inst_vpr.success);
    printf ("heating %i \n", heating);
    LOG_INFO("ier_comb %i", ier_comb);

    if (inst_vpr.success)
        cum_bac.assets.write_last_vpr();

    //VPR  // ------------se combina profili ok e profilo caldo correggo --------------

    if (!ier_comb && heating >= WARM){

        int ier=corr_vpr();
        LOG_INFO("exit status correggo vpr: (1--fallito 0--ok) %i",ier) ; // debug


        //VPR // ------------se la correzione è andata bene e il profilo è 'fresco' stampo profilo con data-------

        if ( ! ier && inst_vpr.success)
            ier_stampa_vpr=stampa_vpr();
    }
}

namespace {
struct CartData
{
    Image<double> azimut;
    Image<double> range;

    CartData(int max_bin=512)
        : azimut(max_bin), range(max_bin)
    {
        for(int i=0; i<max_bin; i++)
            for(int j=0; j<max_bin; j++)
            {
                range(i, j) = hypot(i+.5,j+.5) ;
                azimut(i, j) = 90. - atan((j+.5)/(i+.5)) * M_1_PI*180.;
            }
    }
};
}

SingleCart::SingleCart(unsigned max_bin)
    : max_bin(max_bin),
      cart(max_bin*2) 
{
}

void SingleCart::creo_cart(const Volume <double>& volume, unsigned el_index)
{
    LOG_CATEGORY("radar.singlecart");

    //matrici per ricampionamento cartesiano
    //int x,y,irange,az,iaz,az_min,az_max,cont;
    int x,y,iaz,az_min,az_max;
    float az;
    CartData cd(max_bin);

    for(unsigned i=0; i<max_bin *2; i++)
        for(unsigned j=0; j<max_bin *2; j++)
            cart(i, j) = MISSING;

    LOG_INFO("Creo_cart - %u", max_bin);

    for(unsigned quad=0; quad<4; quad++)
        for(unsigned i=0; i<max_bin; i++)
            for(unsigned j=0; j<max_bin; j++)
            {
                unsigned irange = (unsigned)round(cd.range(i, j));
                if (irange >= max_bin)
                    continue;
                switch(quad)
                {
                    case 0:
                        x = max_bin + i;
                        y = max_bin - j;
                        az = cd.azimut(i, j);
                        break;
                    case 1:
                        x = max_bin + j;
                        y = max_bin + i;
                        az = cd.azimut(i, j) + 90.;
                        break;
                    case 2:
                        x = max_bin - i;
                        y = max_bin + j;
                        az = cd.azimut(i, j) + 180.;
                        break;
                    case 3:
                        x = max_bin - j;
                        y = max_bin - i;
                        az = cd.azimut(i, j)+270.;
                        break;
                }

                az_min = (int)((az - .45)/.9);
                az_max = ceil((az + .45)/.9);


                if(az_min < 0)
                {
                    az_min = az_min + NUM_AZ_X_PPI;
                    az_max = az_max + NUM_AZ_X_PPI;
                }
                for(iaz = az_min; iaz<az_max; iaz++){
                    // Enrico: cerca di non leggere fuori dal volume effettivo
                    unsigned char sample = 0;
                    if (irange < volume[el_index].beam_size)
                        sample = max((unsigned char) (volume[el_index].get(iaz%NUM_AZ_X_PPI, irange)), (unsigned char)1);   // il max serve perchè il valore di MISSING è 0
                    if(cart(y, x) <= sample)  cart(y, x) = sample;
                }
            }
}

void SingleCart::write_out(Assets& assets, const std::string tagname, const std::string format)
{
    if (getenv("DIR_DEBUG") == NULL) return;
    assets.write_gdal_image(cart, "DIR_DEBUG", tagname.c_str(), format.c_str());
}

}

char *PrendiOra()
{
    time_t clock;
    struct tm *tempo;

    clock=time(0);

    tempo=gmtime(&clock);

    return   asctime(tempo);
}

void prendo_tempo()
{
    static time_t time_tot = 0,time1 = 0,time2 = 0;
    LOG_CATEGORY("radar.timing");

    if(time1 == 0){
        time1=time(&time1);
        time_tot = time1;
    }
    time2 = time(&time2);

    LOG_INFO("tempo parziale %ld ---- totale %ld", time2-time1, time2-time_tot);
    time1=time2;
    return;
}
