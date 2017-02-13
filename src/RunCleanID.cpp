
#include <iostream>
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
printf("Chiamo cleaner senza zdr\n");
        //for (unsigned i = 0; i < 1; ++i){
        for (unsigned i = 0; i < full_volume_z.size(); ++i){
printf("Creo scan per output cleaner\n");
	    full_volume_cleanID.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
printf("Ora chiamo evaluateCleanID\n");
            radarelab::algo::Cleaner::evaluateCleanID(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i),full_volume_cleanID.at(i),i);
	    task="Cleaner base";
        }
      } else {
printf("Chiamo cleaner con zdr\n");
        for (unsigned i = 0; i < full_volume_z.size(); ++i){
            algo::Cleaner::clean(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i),full_volume_zdr.at(i),i);
            algo::Cleaner::clean(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i),full_volume_zdr.at(i),i+100);
        }
      }
    }


	volume::ODIMStorer storer;
	storer.store_quality_uchar((Volume<unsigned char>*)(&full_volume_cleanID));
	storer.storeQuality(argv[2],task );
	cout<<endl<<"Fine"<<endl;
}
