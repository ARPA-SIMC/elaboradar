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
	//string radar_name = "SPC";
	string radar_name = argv[3];
	bool init_sqi = false;
        
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_WRAD,&full_volume_wrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_RHOHV,&full_volume_rohv);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SNR,&full_volume_snr);
	//try{
	//loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);
	//}
	//catch(const std::exception &exc){
	//init_sqi = true;
	//cout<<"impossible loading SQI: initialize to 0"<<endl;
	//}
	if(radar_name=="GAT")
	  init_sqi=true;
	else
	  loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);

	loader_all.load(argv[1]);
        cout<<argv[1]<<endl;
	//unica funzione fuzzy
	if( !full_volume_wrad.empty() && !full_volume_vrad.empty()){

	  unsigned last = full_volume_z.size() -1;
	  
	  if (full_volume_zdr.empty()){
	    //inizializzo matrice di zeri
	    //cout<<"ZDR empty"<<endl;
	    //for (unsigned i=0; i<full_volume_z.size();++i){
	    //full_volume_zdr.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
	    //full_volume_zdr.at(i).row(i).setZero(); //0 sarebbe undetected, puoi scegliere altro valore
	    //}
	    is_zdr = false;
	  }
	  cout<<"full volume zdr size = "<<full_volume_zdr.size()<<" and z size "<<full_volume_z.size()<<endl;
	  cout<<"is zdr="<<is_zdr<<endl;
	  for (unsigned i=0; i<full_volume_z.size();++i){//1 anziche full_volume_z.size()
	  //for (unsigned i=0; i<1;++i){
	      full_volume_cleanID.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
	      full_volume_diffprob.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
	      if(init_sqi){
		full_volume_sqi.append_scan(full_volume_z.at(i).beam_count,full_volume_z.at(i).beam_size,full_volume_z.at(i).elevation, full_volume_z.at(i).cell_size);
		full_volume_sqi.at(i).setZero();
	      }

	      volume::Scans<double> Texture;

              //calcolo texture V:
	      if(i< last){
                volume::Scans<double> Input,Input2;        ;
		Input.push_back(full_volume_z.at(i));
	        Input2.push_back(full_volume_z.at(i+1));
	        radarelab::volume::textureVD(Input, Input2, Texture, true);
	        Texture.at(0).nodata=65535.;
	        Texture.at(0).undetect=0.;
		Texture.at(0).gain=200./65535.;
	        Texture.at(0).offset=-100.;
	        //Z_VD.push_back(Texture.at(0));
		cout<<"it="<<i<<", Texture size = "<<Texture.size()<<" "<<Texture.at(0).size()<<endl;
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
	        radarelab::algo::Cleaner::evaluateClassID(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i), full_volume_zdr.at(i), full_volume_rohv.at(i), full_volume_sqi.at(i), full_volume_snr.at(i), Texture.at(0), full_volume_cleanID.at(i), full_volume_diffprob.at(i), full_volume_vrad.at(i).undetect , radar_name, i, true);

		//full_volume_diffprob.at(i).gain=100.0;

	      }else{
		radarelab::algo::Cleaner::evaluateClassID(full_volume_z.at(i), full_volume_wrad.at(i), full_volume_vrad.at(i), full_volume_cleanID.at(i), full_volume_vrad.at(i).undetect, radar_name, i);
	      }
              
	      task="Cleaner base";
	      double new_value=full_volume_z.at(last).nodata;
	      cout<<"novalue"<<new_value<<endl;
	      if (argc == 4 && ! sscanf(argv[3], " --Use_undetect") ) new_value=full_volume_z.at(last).undetect;
	      for (unsigned ii = 0; ii < full_volume_z.at(i).beam_count; ++ii)
                for (unsigned ib = 0; ib < full_volume_z.at(i).beam_size; ++ib) {
		  
     	          if(full_volume_cleanID.at(i)(ii,ib) ) 
		    full_volume_z.at(i)(ii,ib)= new_value;
		  //cout<<"full_clean_ID(i)(ii,ib)= "<<full_volume_cleanID.at(i)(ii,ib)<<endl;
        	}
	    	      
	      }
	  

	}

      std::cout<<"Finito Cleaner, salvo risultati"<<std::endl;
	volume::ODIMStorer storer;
	storer.replace_quantity((Volume<double>*)(&full_volume_z));
	cout<<"replaced quantity"<<endl;
	storer.store_quantity_uchar((Volume<unsigned char>*)(&full_volume_cleanID));
	cout<<"stored_quantity_uchar"<<endl;
	storer.store_quantity_fp((Volume<double>*)(&full_volume_diffprob));
	cout<<"stored quantity"<<endl;
	storer.store(argv[2]);
	cout<<endl<<"Fine"<<endl;
}
