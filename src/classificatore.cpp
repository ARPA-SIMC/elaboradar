/**
 *  @file
 *  @defgroup HydroClass 
 *  @brief Hydrometeor Calssifier
 *  @details Main program for HydroMeteor classification based on radar polarimetry. 
 *      Classification is based on the HCA from Park et al. (2009).
 * 	Input data are Z, Zdr, rhohv and phidp and are provided through
 * 	odim volume data file.
 */

#include <iostream>
#include <radarelab/volume.h>
#include <radarelab/odim.h>

#include "classifier.h"
#include <radarelab/algo/elabora_volume.h>

using namespace radarelab;
using namespace std;

int main(int argc,char* argv[])
{
//const elaboradar::Site& sito(elaboradar::Site::get("SPC"));
	
	//Volume<double> volume;
	//volume::ODIMLoader loader(sito, false, 1024);

	volume::classifier classificatore(argv[1]);

	cout<<"riempito classificatore"<<endl;
	classificatore.compute_derived_volumes();
	cout<<"calcolati i volumi derivati"<<endl;
	classificatore.HCA_Park_2009();
//	classificatore.print_ppi_class();

	volume::ODIMStorer storer;
	storer.store_quantity_int((Volume<int>*)(&classificatore.vol_hca));
	storer.store(argv[1]);
/*

	cout<<"parte interessante"<<endl;
	Volume<double> dummy;

	Volume<double> fil_rho;
	Volume<double> fil_zdr;
	Volume<double> fil_z;
	Volume<double> fil_phi;
	cout<<"filtro"<<endl;
	bool undet=false;
	filter(classificatore.vol_rhohv,fil_rho,2000.,0.,undet);
	filter(classificatore.vol_zdr,fil_zdr,2000.,0.,undet);
	filter(classificatore.vol_z,fil_z,1000.,0.,undet);
	filter(classificatore.vol_phidp_6km,fil_phi,0.,0.,undet);
	volume::ODIMStorer store2(sito, false, 1024);
	cout<<"salvo"<<endl;
	store2.store_quantity_fp(&fil_rho);
	store2.store_quantity_fp(&fil_zdr);
	store2.store_quantity_fp(&fil_z);
	store2.store_quantity_fp(&fil_phi);
	cout<<"scrivo"<<endl;
	store2.store("vol_rad.h5");	
	filter(classificatore.vol_rhohv,fil_rho,0.,3.,undet);
	filter(classificatore.vol_zdr,fil_zdr,0.,3.,undet);
	filter(classificatore.vol_z,fil_z,0.,3.,undet);
	volume::ODIMStorer store3(sito, false, 1024);
	cout<<"salvo"<<endl;
	store3.store_quantity_fp(&fil_rho);
	store3.store_quantity_fp(&fil_zdr);
	store3.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store3.store("vol_az.h5");
	filter(classificatore.vol_rhohv,fil_rho,2000.,3.,undet);
	filter(classificatore.vol_zdr,fil_zdr,2000.,3.,undet);
	filter(classificatore.vol_z,fil_z,1000.,3.,undet);
	volume::ODIMStorer store4(sito, false, 1024);
	cout<<"salvo"<<endl;
	store4.store_quantity_fp(&fil_rho);
	store4.store_quantity_fp(&fil_zdr);
	store4.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store4.store("vol_azrad.h5");



	filter(classificatore.vol_rhohv,dummy,2000.,0.,undet);
	filter(dummy,fil_rho,2000.,0.,undet);
	filter(classificatore.vol_zdr,dummy,2000.,0.,undet);
	filter(dummy,fil_zdr,2000.,0.,undet);
	filter(classificatore.vol_z,fil_z,1000.,0.,undet);
	volume::ODIMStorer store22(sito, false, 1024);
	cout<<"salvo"<<endl;
	store22.store_quantity_fp(&fil_rho);
	store22.store_quantity_fp(&fil_zdr);
	store22.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store22.store("vol_rd_2.h5");
	
	filter(classificatore.vol_rhohv,dummy,0.,3.,undet);
	filter(dummy,fil_rho,0.,5.,undet);
	filter(classificatore.vol_zdr,dummy,0.,3.,undet);
	filter(dummy,fil_zdr,0.,5.,undet);
	filter(classificatore.vol_z,fil_z,0.,3.,undet);
	volume::ODIMStorer store33(sito, false, 1024);
	cout<<"salvo"<<endl;
	store33.store_quantity_fp(&fil_rho);
	store33.store_quantity_fp(&fil_zdr);
	store33.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store33.store("vol_az_2.h5");

	filter(classificatore.vol_rhohv,dummy,2000.,3.,undet);
	filter(dummy,fil_rho,2000.,5.,undet);
	filter(classificatore.vol_zdr,dummy,2000.,3.,undet);
	filter(dummy,fil_zdr,2000.,5.,undet);
	filter(classificatore.vol_z,fil_z,1000.,3.,undet);
	volume::ODIMStorer store44(sito, false, 1024);
	cout<<"salvo"<<endl;
	store44.store_quantity_fp(&fil_rho);
	store44.store_quantity_fp(&fil_zdr);
	store44.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store44.store("vol_azrd_2.h5");

	Volume<double> vol_sdphidp;
	Volume<double> vol_sdz;
	textureSD(classificatore.vol_phidp,vol_sdphidp,2000.,0.,true);
	textureSD(classificatore.vol_z,vol_sdz,1000.,0.,true);
	volume::ODIMStorer storeSD(sito,false,1024);
	storeSD.store_quantity_fp(&vol_sdphidp);
	storeSD.store_quantity_fp(&vol_sdz);
	storeSD.store("vol_sd.h5");
*/	

	cout<<endl<<"Fine"<<endl;
	
}
