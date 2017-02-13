
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

//	printf("il nome del mio file è %s\n", pathname.c_str());

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
//printf("Chiamo cleaner senza zdr\n");
        //for (unsigned i = 0; i < 1; ++i){
        for (unsigned i = 0; i < full_volume_z.size(); ++i){
	    full_volume_cleanID.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
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
    vector <std::string> Sweep(full_volume_cleanID.at(0).beam_count, "" );
    std::string my_time;         // formato "YYYY-MM-DD hh:mm:s
    my_time = Radar::timeutils::absoluteToString(full_volume_z.load_info->acq_date);
    for (unsigned iel = 0; iel< full_volume_cleanID.size(); iel++){
      //VectorXd  conteggi ;
      auto Weather = (full_volume_cleanID.at(iel).array() == 0 && full_volume_z.at(iel).array() >= -30.).rowwise().count() ;
      auto Clutter = (full_volume_cleanID.at(iel).array() == 1 ).rowwise().count() ;
      auto Interf  = (full_volume_cleanID.at(iel).array() >= 2 && full_volume_cleanID.at(iel).array() <= 4 ).rowwise().count() ;
      auto Noise   = (full_volume_cleanID.at(iel).array() == 5 ).rowwise().count() ;
      for (unsigned iray=0; iray < full_volume_cleanID.at(iel).beam_count; iray ++){
        char Ray[50];
	sprintf(Ray,", %5.1f, %5d, %5d, %5d, %5d",full_volume_z.at(iel).elevation, Weather(iray), Clutter(iray), Interf(iray), Noise(iray));
        Sweep[iray] += Ray;	
//	printf("%s,%5.1f,%6.1f, %5d, %5d, %5d, %5d\n",my_time.c_str(), full_volume_z.at(iel).elevation, full_volume_z.at(iel).azimuths_real(iray), 
//		Weather(iray), Clutter(iray), Interf(iray), Noise(iray));
      }
    }
	printf("DateTime, Azimuth, Elev_01, Weather_E01,Clutter_E01, Iterf_E01,Noise_E01, Elev_02, Weather_E02,Clutter_E02, Iterf_E02,Noise_E02, Elev_03, Weather_E03,Clutter_E03, Iterf_E03,Noise_E03, Elev_04, Weather_E04,Clutter_E04, Iterf_E04,Noise_E04, Elev_05, Weather_E05,Clutter_E05, Iterf_E05,Noise_E05, Elev_06, Weather_E06,Clutter_E06, Iterf_E06,Noise_E06\n");
    for (unsigned iray=0; iray < Sweep.size(); iray ++){
	printf("%s,%6.1f %s\n",my_time.c_str(), full_volume_z.at(0).azimuths_real(iray),Sweep[iray].c_str()); 
    }
}
