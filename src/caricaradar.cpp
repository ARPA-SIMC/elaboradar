#include<iostream>
#include"volume.h"
#include"volume/odim.h"
#include"volume/loader.h"
#include"site.h"
#include "volume/resample.h"

#include"classifier.h"

using namespace cumbac;
using namespace std;

int main(int argc,char* argv[])
{
	const Site& sito(Site::get("GAT"));
	//site.name="GAT";
	
	Volume<double> volume;
	volume::ODIMLoader loader(sito, false, false, 1024);
	loader.coherent_loader=true;

	volume::LoadInfo load_info;
	loader.load_info = &load_info;

	volume::Scans<double> full_volume;
	loader.vol_z = &full_volume;

	if(loader.load(argv[1],argv[2])) std::cout<<"tutto bene"<<std::endl;
	else std::cout<<"tutto male"<<std::endl;

	volume::volume_resample<double>(full_volume, loader.azimuth_maps, volume, volume::merger_closest<double>);
	
	cout<<"elevations "<<volume.size()<<endl;
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

	volume::classifier classificatore(argv[1],sito);
	classificatore.compute_derived_volumes();

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

	classificatore.HCA_Park_2009();

	cout<<endl<<"Fine"<<endl;

}
