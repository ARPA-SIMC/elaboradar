
#include <iostream>
#include <cstring>
#include <radarelab/volume.h>
#include <radarelab/odim.h>
#include <radarlib/radar.hpp>
#include <sstream>
#include <radarelab/image.h>
#include <radarelab/algo/azimuth_resample.h>
#include <radarelab/algo/cleaner.h>

#include <radarelab/algo/elabora_volume.h>

using namespace radarelab;
using namespace std;

using namespace volume;
namespace odim = OdimH5v21;

int main(int argc,char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <h5-volume-input> <h5-volume-output> [--Use_undetect] \n", argv[0]);
		exit(1);
	}
	
	std::string pathname = argv[1];	

	printf("il nome del mio file è %s\n", pathname.c_str());

	volume::ODIMLoader loader_all;

	volume::Scans<double> full_volume_z;
	volume::Scans<double> full_volume_zdr;
	volume::Scans<double> full_volume_vrad;
	volume::Scans<double> full_volume_wrad;
	volume::Scans<unsigned char> full_volume_cleanID;
	std::string task;

	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_WRAD,&full_volume_wrad);

	loader_all.load(argv[1]);
    if ( !full_volume_wrad.empty() && !full_volume_vrad.empty())
    {
      if (full_volume_zdr.empty())
      {
	unsigned last = full_volume_z.size() -1; 
        for (unsigned i = 0; i < full_volume_z.size(); ++i){
	    full_volume_cleanID.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
            radarelab::algo::Cleaner::evaluateCleanID(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i),full_volume_cleanID.at(i),i);
	    task="Cleaner base";
	    double new_value=full_volume_z.at(last).nodata;
	    if (argc == 4 && ! sscanf(argv[3], " --Use_undetect") ) new_value=full_volume_z.at(last).undetect;
	    for (unsigned ii = 0; ii < full_volume_z.at(i).beam_count; ++ii)
        	for (unsigned ib = 0; ib < full_volume_z.at(i).beam_size; ++ib) {
	//	    printf(" %4d %4d %4d %4d %5.2f %5.2f %5.2f %5.2f  ---> ", i,ii,ib, full_volume_cleanID.at(i)(ii,ib), full_volume_z.at(i)(ii,ib),
	//			    	full_volume_z.at(last).nodata, full_volume_z.at(last).gain, full_volume_z.at(last).offset);
     		    if(full_volume_cleanID.at(i)(ii,ib) ) 
			    full_volume_z.at(i)(ii,ib)= new_value;
	//	    printf(" %6.2f \n", full_volume_z.at(i)(ii,ib));
        	}
	}
      } else {
        for (unsigned i = 0; i < full_volume_z.size(); ++i){
            algo::Cleaner::clean(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i),full_volume_zdr.at(i),i,true);
            algo::Cleaner::clean(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i),full_volume_zdr.at(i),i+100,true);
        }
      }
    }
 

      std::cout<<"Finito Cleaner, salvo risultati"<<std::endl;
	volume::ODIMStorer storer;
	storer.replace_quantity((Volume<double>*)(&full_volume_z));
	storer.store_quantity_uchar((Volume<unsigned char>*)(&full_volume_cleanID));
	storer.store(argv[2]);
	cout<<endl<<"Fine"<<endl;
}
