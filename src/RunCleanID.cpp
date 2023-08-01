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
#include <tclap/CmdLine.h>

using namespace radarelab;
using namespace std;

using namespace volume;
namespace odim = OdimH5v21;

int main(int argc,char* argv[])
{
	TCLAP::CmdLine cmd("RunCleanID ", ' ', "0.1" );

	TCLAP::UnlabeledValueArg<std::string> cmd_vol_input("h5_volume_input", "input hdf5 volume", true, "NULL", "h5-volume-output");
	cmd.add(cmd_vol_input);

	TCLAP::UnlabeledValueArg<std::string> cmd_vol_output("h5_volume_output", "post-processed overwritten input hdf5 volume", true, "NULL", "h5-volume-output");
	cmd.add(cmd_vol_output);

	TCLAP::ValueArg<std::string> cmd_fuzzy_path("F", "fuzzy-path", "Optional: Set path of fuzzy logic files and clutter maps. \n Default: /usr/share/elaboradar ", false, FUZZY_PATH, "path");
	cmd.add(cmd_fuzzy_path);

	TCLAP::SwitchArg cmd_use_undetect("U", "use-undetect", "Optional: Use undetect value (-31.15 dBZ) as DBZH replacing value for pixels classified as non meteorological echo. \nIf not passed, nodata value is used instead (99.95 dBZ)", false);
	cmd.add(cmd_use_undetect);

	cmd.parse(argc,argv);

	const char* fuzzy_path = cmd_fuzzy_path.getValue().c_str();
	cout<<"fuzzypath="<<fuzzy_path<<" "<<endl;

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
	//volume::Scans<double> Z_VD;
	std::string task;
	bool is_zdr=true;
	bool init_sqi = false;

	volume::Scans<double> SDZ6;
	SDZ6.quantity="SDZ6";
        
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_WRAD,&full_volume_wrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_RHOHV,&full_volume_rohv);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SNR,&full_volume_snr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);

	loader_all.load(cmd_vol_input.getValue());
	string radar_name = full_volume_z.load_info->source_name;
	radar_name = radar_name.substr(2);
	transform(radar_name.begin(), radar_name.end(), radar_name.begin(),::toupper);
	
	//unica funzione fuzzy
	if( !full_volume_wrad.empty() && !full_volume_vrad.empty()){

	  unsigned last = full_volume_z.size() -1;
	  cout<<"last="<<last<<endl;
	  
	  if (full_volume_zdr.empty()){
	    is_zdr = false;
	  }
	  if(full_volume_sqi.empty()){
	    init_sqi=true;
	  }
	  cout<<"full volume zdr size = "<<full_volume_zdr.size()<<" and z size "<<full_volume_z.size()<<endl;
	  cout<<"is zdr="<<is_zdr<<endl;
	  for (unsigned i=0; i<full_volume_z.size();++i){//1 anziche full_volume_z.size()
	  //for (unsigned i=0; i<1;++i){
	    cout<<"elev="<<full_volume_z.at(i).elevation<<endl;
	    full_volume_cleanID.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size, 0);
	      full_volume_diffprob.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
	      if(init_sqi){
		full_volume_sqi.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
		full_volume_sqi.at(i).setZero();
	      }
	      //SDZ6.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);

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

		//full_volume_diffprob.at(i).gain=100.0;

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

	      volume::Scans<double> sdz6, z_cur;
	      z_cur.push_back(full_volume_z.at(i));
	      textureSD(z_cur,sdz6,6000., false);
	      sdz6.at(0).nodata=65535.;
	      sdz6.at(0).undetect=0.;
	      SDZ6.push_back(sdz6.at(0));
	    	      
	      }
	  

	}

      std::cout<<"Finito Cleaner, salvo risultati"<<std::endl;
	volume::ODIMStorer storer;
	storer.replace_quantity((Volume<double>*)(&full_volume_z));
	cout<<"replaced quantity"<<endl;
	storer.store_quantity_fp((Volume<double>*)(&SDZ6));
	cout<<"stored quantity"<<endl;
	storer.store_quantity_uchar((Volume<unsigned char>*)(&full_volume_cleanID));
	cout<<"stored_quantity_uchar"<<endl;
	storer.store_quantity_fp((Volume<double>*)(&full_volume_diffprob));
	cout<<"stored quantity"<<endl;
	storer.store(cmd_vol_output.getValue());
	cout<<endl<<"Fine"<<endl;
}
