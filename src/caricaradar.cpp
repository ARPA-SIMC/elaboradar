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
	//loader.vol_z = &full_volume;

	loader.request_quantity(argv[2],&full_volume);

	//if(loader.load(argv[1])) std::cout<<"tutto bene"<<std::endl;
	//else std::cout<<"tutto male"<<std::endl;
	loader.load(argv[1]);

//	cout<<volume.beam_count<<" "<<full_volume[0].beam_count<<endl;

	volume::volume_resample<double>(full_volume, loader.azimuth_maps, volume, volume::merger_closest<double>);
	cout<<"a me si"<<volume.scan(volume.size()-1).beam_size<<endl;
	cout<<"elevations "<<volume.size()<<endl;

	Volume<double> Z_fil;
	Volume<double> filtro;
//	filter(volume,Z_fil,1000.);
//	filtro.filter(volume,1000.);
//	cout<<"ritornato"<<endl;
	Volume<double> Z_SD;
	Volume<double> rms;
//	textureSD(volume,Z_SD,1000.);
//	rms.textureSD(volume,1000.);
	Volume<double> Z_slope;
//	moving_average_slope(volume,Z_slope,1000.);

	for(unsigned idx=0;idx<1000;idx++)
	{
		//filter(volume,Z_fil,1000.);
		//filtro.filter(volume,1000.);
		//textureSD(volume,Z_SD,1000.);
		rms.textureSD(volume,1000.);
	}




//	cout<<endl;
//	cout<<Z_fil.size()<<endl;
//	for(unsigned rg=0;rg<50;rg++)
//		cout<<fixed<<volume[0](70,rg)<<"\t"<<Z_fil[0](70,rg)<<"\t"<<filtro[0](70,rg)<<"\t"<<Z_SD[0](70,rg)<<"\t"<<rms[0](70,rg)<<"\t"<<Z_slope[0](70,rg)<<endl;

//	cout<<endl;
/*
	for(unsigned el=0;el<volume.size();el++)
	{
		cout<<"el "<<el<<" elev "<<volume.scan(el).elevation
		<<" beams "<<volume.scan(el).beam_count<<" "<<volume.scan(el).rows()
		<<" size "<<volume.scan(el).beam_size<<" "<<volume.scan(el).cols()<<endl; 
	}
	for(unsigned el=0;el<full_volume.size();el++)
	{
		cout<<"el "<<el<<" elev "<<full_volume[el].elevation
		<<" beams "<<full_volume[el].beam_count<<" "<<full_volume[el].rows()
		<<" size "<<full_volume[el].beam_size<<" "<<full_volume[el].cols()<<endl; 
	}

	for(unsigned az=0;az<10;az++)
	{
		for(unsigned rg=0;rg<10;rg++) cout<<volume.scan(2).get(az,rg)<<"\t";
		cout<<endl;
	}
*/
//	volume::classifier classificatore(argv[1],sito);
//	cout<<"riempito classificatore"<<endl;
//	classificatore.compute_derived_volumes();
//	cout<<"calcolati i volumi derivati"<<endl;
/*
	for(unsigned az=0;az<10;az++)
	{
		for(unsigned rg=0;rg<10;rg++)
		{
			cout<<classificatore.vol_phidp_2km.scan(0).get(az,rg)<<"\t";
		}
		cout<<endl;
	}
	cout<<endl;
	for(unsigned az=0;az<10;az++)
	{
		for(unsigned rg=0;rg<10;rg++)
		{
			cout<<classificatore.vol_lkdp_2km.scan(0).get(az,rg)<<"\t";
		}
		cout<<endl;
	}
	cout<<endl;
*/

//	classificatore.HCA_Park_2009();

	cout<<endl<<"Fine"<<endl;

}
