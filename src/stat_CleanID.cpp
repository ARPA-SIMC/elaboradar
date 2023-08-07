
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
	TCLAP::CmdLine cmd("stat_CleanID ", ' ', "1.0" );

	TCLAP::UnlabeledValueArg<std::string> cmd_vol_input("h5_volume_input", "hdf5 volume input", true, "NULL", "h5-volume-input");
	cmd.add(cmd_vol_input);

	TCLAP::ValueArg<std::string> cmd_fuzzy_path("F", "fuzzy-path", "Set path of fuzzy logic files", false, FUZZY_PATH, "path");
	cmd.add(cmd_fuzzy_path);

	TCLAP::ValueArg<int> how_many_elev("N", "how-many-elev", "Number of sweeps to be analised", false, 6, "Number of sweeps to be analised");
	cmd.add(how_many_elev);

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
	bool init_sqi = false;
	std::string task;

	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_WRAD,&full_volume_wrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_RHOHV,&full_volume_rohv);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SNR,&full_volume_snr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);

	loader_all.load(cmd_vol_input.getValue());
	string radar_name = full_volume_z.load_info->source_name;
	radar_name = radar_name.substr(2);
	transform(radar_name.begin(), radar_name.end(), radar_name.begin(),::toupper);

    int N_ELEV = full_volume_z.size() < how_many_elev.getValue() ? full_volume_z.size() : how_many_elev.getValue();
    if ( !full_volume_wrad.empty() && !full_volume_vrad.empty())
    {
      unsigned last = full_volume_z.size() -1;
      if (full_volume_zdr.empty())
      {
	is_zdr = false;
      }
      if (full_volume_sqi.empty()){
	init_sqi = true;
      }
      for (unsigned i = 0; i < N_ELEV ; ++i){
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
		//cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
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
		//cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
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
	      //cout<<"novalue"<<new_value<<endl;
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
    for (unsigned iel = 0; iel< N_ELEV; iel++){
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
      }
    }

    // stampo tabella statistica output
    std::string TabHeader="DateTime, Azimuth";
    for (unsigned i=1; i<=N_ELEV; i++){
	    char Header_Elev[50];
	    sprintf(Header_Elev,", Elev_%2.2d, Weather_E%2.2d,Clutter_E%2.2d, Iterf_E%2.2d",i,i,i,i);
	    TabHeader += Header_Elev;
    }
    printf("%s\n",TabHeader.c_str());
    for (unsigned iray=0; iray < Sweep.size(); iray ++){
	printf("%s,%6.1f %s\n",my_time.c_str(), full_volume_z.at(0).azimuths_real(iray),Sweep[iray].c_str()); 
    }
}
