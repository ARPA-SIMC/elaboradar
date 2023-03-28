
#include <iostream>
#include <radarelab/volume.h>
#include <radarelab/odim.h>
#include <radarlib/radar.hpp>
#include <sstream>
#include <radarelab/image.h>
#include <radarelab/algo/azimuth_resample.h>
#include <radarelab/algo/cleaner.h>
#include <radarelab/algo/elabora_volume.h>
#include <tclap/CmdLine.h>


using namespace radarelab;
using namespace std;

using namespace volume;
namespace odim = OdimH5v21;

int main(int argc,char* argv[])
{
	TCLAP::CmdLine cmd("stat_CleanID ", ' ', "0.1" );

	TCLAP::UnlabeledValueArg<std::string> cmd_vol_input("h5_volume_input", "hdf5 volume input", true, "NULL", "h5-volume-output");
	cmd.add(cmd_vol_input);

	TCLAP::UnlabeledValueArg<std::string> cmd_vol_output("h5_volume_output", "hdf5 volume output", true, "NULL", "h5-volume-output");
	cmd.add(cmd_vol_output);

	TCLAP::UnlabeledValueArg<std::string> cmd_radar("radar_name", "radar name", true, "NULL", "radar-name");
	cmd.add(cmd_radar);

	TCLAP::ValueArg<std::string> cmd_fuzzy_path("F", "fuzzy-path", "Set path of fuzzy logic files", false, FUZZY_PATH "/dati", "path");
	cmd.add(cmd_fuzzy_path);

	TCLAP::SwitchArg cmd_use_undetect("U", "use-undetect", "Use undetect TODO", false);
	cmd.add(cmd_use_undetect);

	cmd.parse(argc,argv);

	const char* fuzzy_path = cmd_fuzzy_path.getValue().c_str();
	volume::ODIMLoader loader_all;

	volume::Scans<double> full_volume_z;
	volume::Scans<double> full_volume_zdr;
	volume::Scans<double> full_volume_vrad;
	volume::Scans<double> full_volume_wrad;
	volume::Scans<double> full_volume_rohv;
	volume::Scans<double> full_volume_sqi;
	volume::Scans<double> full_volume_snr;
	volume::Scans<unsigned char> full_volume_cleanID;
	volume::Scans<double> full_volume_diffprob;
	full_volume_cleanID.quantity="ClassID";
	full_volume_diffprob.quantity="Diffprob";
	bool is_zdr=true;
	string radar_name = cmd_radar.getValue();
	bool init_sqi = false;
	std::string task;

	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_WRAD,&full_volume_wrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_RHOHV,&full_volume_rohv);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SNR,&full_volume_snr);
	if(radar_name=="GAT")
	  init_sqi=true;
	else
	  loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);

	loader_all.load(cmd_vol_input.getValue());

    if ( !full_volume_wrad.empty() && !full_volume_vrad.empty())
    {
      unsigned last = full_volume_z.size() -1;
      if (full_volume_zdr.empty())
      {
	is_zdr = false;
      }
        //for (unsigned i = 0; i < 1; ++i){
      for (unsigned i = 0; i < full_volume_z.size(); ++i){
	full_volume_cleanID.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size, 0);
	  full_volume_diffprob.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
	  if(init_sqi){
	      full_volume_sqi.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
	      full_volume_sqi.at(i).setZero();
	  }
	  volume::Scans<double> Texture;

              //calcolo texture V:
	      if(i< last){
                volume::Scans<double> Input,Input2;
		Input.push_back(full_volume_z.at(i));
	        Input2.push_back(full_volume_z.at(i+1));
	        radarelab::volume::textureVD(Input, Input2, Texture, true);
	        Texture.at(0).nodata=65535.;
	        Texture.at(0).undetect=0.;
		Texture.at(0).gain=200./65535.;
	        Texture.at(0).offset=-100.;
	        //Z_VD.push_back(Texture.at(0));
		cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
		//cout<<"z min"<<std::min(full_volume_z.at(i),100000)<<"z max"<<std::max(full_volume_z.at(i),0)<<endl;
	      }
	      else{
	        Texture.clear();
		//cout<<"it="<<i<<", Texture size = "<<Texture.size()<<endl;
		//M_start = Matrix2D::Zero(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size);
		Texture.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
		Texture.at(0).setZero();
		Texture.at(0).nodata=65535.;
	        Texture.at(0).undetect=0.;
		Texture.at(0).gain=200./65535.;
	        Texture.at(0).offset=-100.;
		cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
		//cout<<"Texture is zero?"<<Texture.at(0)(30,50)<<endl;
	      }

	      if(is_zdr){
	        radarelab::algo::Cleaner::evaluateClassID(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i), full_volume_zdr.at(i), full_volume_rohv.at(i), full_volume_sqi.at(i), full_volume_snr.at(i), Texture.at(0), full_volume_cleanID.at(i), full_volume_diffprob.at(i), full_volume_vrad.at(i).undetect , radar_name, fuzzy_path, i, true);

	      }else{
		radarelab::algo::Cleaner::evaluateClassID(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i), full_volume_cleanID.at(i), full_volume_vrad.at(i).undetect, radar_name, fuzzy_path, i);
	      }
	      task="Cleaner base";
	      double new_value=full_volume_z.at(last).nodata;
	      if (cmd_use_undetect.getValue()) new_value=full_volume_z.at(last).undetect;
	      cout<<"novalue"<<new_value<<endl;
	      for (unsigned ii = 0; ii < full_volume_z.at(i).beam_count; ++ii)
                for (unsigned ib = 0; ib < full_volume_z.at(i).beam_size; ++ib) {
		  
     	          if(full_volume_cleanID.at(i)(ii,ib) ) 
		    full_volume_z.at(i)(ii,ib)= new_value;
		  //cout<<"full_clean_ID(i)(ii,ib)= "<<full_volume_cleanID.at(i)(ii,ib)<<endl;
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
      //auto Noise   = (full_volume_cleanID.at(iel).array() == 5 ).rowwise().count() ;
      for (unsigned iray=0; iray < full_volume_cleanID.at(iel).beam_count; iray ++){
        char Ray[50];
	//sprintf(Ray,", %5.1f, %5d, %5d, %5d, %5d",full_volume_z.at(iel).elevation, Weather(iray), Clutter(iray), Interf(iray), Noise(iray));
	sprintf(Ray,", %5.1f, %5d, %5d, %5d",full_volume_z.at(iel).elevation, Weather(iray), Clutter(iray), Interf(iray));
        Sweep[iray] += Ray;	
//	printf("%s,%5.1f,%6.1f, %5d, %5d, %5d, %5d\n",my_time.c_str(), full_volume_z.at(iel).elevation, full_volume_z.at(iel).azimuths_real(iray), 
//		Weather(iray), Clutter(iray), Interf(iray), Noise(iray));
      }
    }
	printf("DateTime, Azimuth, Elev_01, Weather_E01,Clutter_E01, Iterf_E01, Elev_02, Weather_E02,Clutter_E02, Iterf_E02, Elev_03, Weather_E03,Clutter_E03, Iterf_E03, Elev_04, Weather_E04,Clutter_E04, Iterf_E04, Elev_05, Weather_E05,Clutter_E05, Iterf_E05, Elev_06, Weather_E06,Clutter_E06, Iterf_E06\n");
    for (unsigned iray=0; iray < Sweep.size(); iray ++){
	printf("%s,%6.1f %s\n",my_time.c_str(), full_volume_z.at(0).azimuths_real(iray),Sweep[iray].c_str()); 
    }
}
