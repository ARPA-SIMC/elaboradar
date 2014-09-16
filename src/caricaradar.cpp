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
	cout<<endl<<"Fine"<<endl;
}
