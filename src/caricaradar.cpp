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
	cout<<"riempito classificatore"<<endl;
	classificatore.compute_derived_volumes();
	cout<<"calcolati i volumi derivati"<<endl;
	classificatore.HCA_Park_2009();
	classificatore.print_ppi_class();

	volume::ODIMStorer storer(sito, false, 1024);

// CHECK: It is safe to pass pointers to volumes?
//	Volume<volume::bit64>* puntatore;
//	puntatore = new Volume<volume::bit64>[1];
//	puntatore = &classificatore.vol_z;
//	storer.store_quantity(puntatore);
	storer.store_quantity_int((Volume<int>*)(&classificatore.vol_hca));
	storer.store_quantity_fp(&classificatore.vol_z);
	cout<<storer.to_store_fp.size()<<endl;
	cout<<storer.to_store_int.size()<<endl;
	cout<<classificatore.vol_z.scan(0).get(0,0)<<" "<<storer.to_store_fp[0]->scan(0).get(0,0)<<endl;
	cout<<classificatore.vol_z.scan(0).get(20,200)<<" "<<storer.to_store_fp[0]->scan(0).get(20,200)<<endl;
	cout<<classificatore.vol_hca.scan(0).get(0,0)<<" "<<storer.to_store_int[0]->scan(0).get(0,0)<<endl;
	cout<<classificatore.vol_hca.scan(0).get(20,200)<<" "<<storer.to_store_int[0]->scan(0).get(20,200)<<endl;

	storer.store("volume_scritto.h5");

	cout<<endl<<"Fine"<<endl;
	
}
