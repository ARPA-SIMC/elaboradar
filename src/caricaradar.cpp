#include<iostream>
#include"volume.h"
#include"volume/odim.h"
#include"volume/loader.h"
#include"site.h"
#include "volume/resample.h"

#include"classifier.h"
#include"algo/elabora_volume.h"

using namespace elaboradar;
using namespace std;

int main(int argc,char* argv[])
{
	const Site& sito(Site::get("SPC"));
	
	//Volume<double> volume;
	//volume::ODIMLoader loader(sito, false, 1024);

	volume::classifier classificatore(argv[1],sito);
/*
	cout<<"riempito classificatore"<<endl;
	classificatore.compute_derived_volumes();
	cout<<"calcolati i volumi derivati"<<endl;
	classificatore.HCA_Park_2009();
	classificatore.print_ppi_class();

	volume::ODIMStorer storer(sito, false, 1024);
	storer.store_quantity_int((Volume<int>*)(&classificatore.vol_hca));
	storer.store(argv[1]);
*/

	cout<<"parte interessante"<<endl;
	Volume<double> fil_rho;
	Volume<double> fil_zdr;
	Volume<double> fil_z;
	cout<<"filtro"<<endl;
	bool undet=false;
//	textureSD(classificatore.vol_rhohv,fil_rho,0.,0.,undet);
//	textureSD(classificatore.vol_zdr,fil_zdr,0.,0.,undet);
	textureSD(classificatore.vol_z,fil_z,1000.,0.,undet);
	volume::ODIMStorer store2(sito, false, 1024);
	cout<<"salvo"<<endl;
//	store2.store_quantity_fp(&fil_rho);
//	store2.store_quantity_fp(&fil_zdr);
	store2.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store2.store("vol_rad.h5");
	
//	textureSD(classificatore.vol_rhohv,fil_rho,0.,3.,undet);
//	textureSD(classificatore.vol_zdr,fil_zdr,0.,3.,undet);
	textureSD(classificatore.vol_z,fil_z,0.,3.,undet);
	volume::ODIMStorer store3(sito, false, 1024);
	cout<<"salvo"<<endl;
//	store3.store_quantity_fp(&fil_rho);
//	store3.store_quantity_fp(&fil_zdr);
	store3.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store3.store("vol_az.h5");

//	textureSD(classificatore.vol_rhohv,fil_rho,0.,3.,undet);
//	textureSD(classificatore.vol_zdr,fil_zdr,0.,3.,undet);
	textureSD(classificatore.vol_z,fil_z,1000.,3.,undet);
	volume::ODIMStorer store4(sito, false, 1024);
	cout<<"salvo"<<endl;
//	store4.store_quantity_fp(&fil_rho);
//	store4.store_quantity_fp(&fil_zdr);
	store4.store_quantity_fp(&fil_z);
	cout<<"scrivo"<<endl;
	store4.store("vol_azrad.h5");

	cout<<endl<<"Fine"<<endl;
	
}
