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
	const Site& sito(Site::get("GAT"));
	//site.name="GAT";
	
	Volume<double> volume;
	volume::ODIMLoader loader(sito, false, 1024);

	volume::Scans<double> full_volume;
	loader.request_quantity(argv[2],&full_volume);
	loader.load(argv[1]);
	volume::volume_resample<double>(full_volume, loader.azimuth_maps, volume, volume::merger_closest<double>);
	cout<<"a me si"<<volume.scan(volume.size()-1).beam_size<<endl;
	cout<<"elevations "<<volume.size()<<endl;

	Volume<double> Z_fil;
	Volume<double> filtro;
//	filter(volume,Z_fil,1000.,false);
//	cout<<"ritornato"<<endl;
	Volume<double> Z_SD;
	Volume<double> rms;
//	textureSD(volume,Z_SD,1000.,false);
	Volume<double> Z_slope;
//	moving_average_slope(volume,Z_slope,1000.,false);

//	lin2dB(volume,Z_fil);
//	lin2dB(volume);

	//for(unsigned idx=0;idx<1000;idx++)
	//{
//		textureSD(volume,Z_fil,1000.,false);
//		textureSD(volume,Z_SD,1000.,false);
	//}




//	cout<<endl;
//	cout<<Z_fil.size()<<endl;
//	for(unsigned rg=0;rg<50;rg++)
//		cout<<fixed<<volume[0](0,rg)<<"\t"<<Z_fil[0](0,rg)<<endl;//"\t"<<filtro[0](70,rg)<<"\t"<<Z_SD[0](70,rg)<<"\t"<<rms[0](70,rg)<<"\t"<<Z_slope[0](70,rg)<<endl;
//	cout<<endl;


	volume::classifier classificatore(argv[1],sito);
	cout<<"riempito classificatore"<<endl;
	classificatore.compute_derived_volumes();
	cout<<"calcolati i volumi derivati"<<endl;

	classificatore.HCA_Park_2009();
	cout<<endl<<"Fine"<<endl;

}
