
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
	if (argc < 1)
	{
		fprintf(stderr, "Usage: %s <h5-volume-input> [--Use_undetect] \n", argv[0]);
		exit(1);
	}
	
	std::string pathname = argv[1];	

	printf("il nome del mio file è %s\n", pathname.c_str());

	volume::ODIMLoader loader_all;

	volume::Scans<double> full_volume_z;
	volume::Scans<double> full_volume_zdr;
	volume::Scans<double> full_volume_vrad;
	volume::Scans<double> full_volume_wrad;
	volume::Scans<double> full_volume_sqi;
	volume::Scans<unsigned char> full_volume_cleanID;

	volume::Scans<double> Z_SD2D;
	Z_SD2D.quantity="DBZH_SD2D";
	volume::Scans<double> Z_SDRay_9;
	Z_SDRay_9.quantity="DBZH_SDRay_9";
	volume::Scans<double> Z_SDRay_21;
	Z_SDRay_21.quantity="DBZH_SDRay_21";
	volume::Scans<double> Z_SDAz;
	Z_SDAz.quantity="DBZH_SDAz";
	volume::Scans<double> ZDR_SD2D;
	ZDR_SD2D.quantity="ZDR_SD2D";
	volume::Scans<double> Z_VD;
	Z_VD.quantity="Z_VD";

	std::string task;
	bool is_zdr=false;

	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_WRAD,&full_volume_wrad);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_SQI,&full_volume_sqi);



	loader_all.load(argv[1]);

    
        for (unsigned i = 0; i < full_volume_z.size(); ++i){
  	    volume::Scans<double> Input, Input_zdr, Texture;
	    Input.push_back(full_volume_z.at(i));
// Texture Z 2D 
	    radarelab::volume::textureSD( Input, Texture, 1000. , 3.,false);
	    Texture.at(0).nodata=65535.;
	    Texture.at(0).undetect=0.;
	    Z_SD2D.push_back(Texture.at(0));
// Texture Z 1D Ray 9 bins
	    radarelab::volume::textureSD( Input,Texture, Input.at(0).cell_size*9 , 360./Input.at(0).beam_count,true);
	    Texture.at(0).nodata=65535.;
	    Texture.at(0).undetect=0.;
	    Z_SDRay_9.push_back(Texture.at(0));
// Texture Z 1D Ray 21 bins
	    radarelab::volume::textureSD( Input,Texture, Input.at(0).cell_size*21 , 360./Input.at(0).beam_count,true);
	    Texture.at(0).nodata=65535.;
	    Texture.at(0).undetect=0.;
	    Z_SDRay_21.push_back(Texture.at(0));
// Texture Z 1D Az 5 rays
	    radarelab::volume::textureSD( Input,Texture, Input.at(0).cell_size , 5*360./Input.at(0).beam_count,true);
	    Texture.at(0).nodata=65535.;
	    Texture.at(0).undetect=0.;
	    Z_SDAz.push_back(Texture.at(0));
// Texture ZDR 2D
	    if( ! full_volume_zdr.empty()){
	    Input_zdr.push_back(full_volume_zdr.at(i));
	    radarelab::volume::textureSD( Input_zdr, Texture, 1000. , 3.,false);
	    Texture.at(0).nodata=65535.;
	    Texture.at(0).undetect=0.;
	    ZDR_SD2D.push_back(Texture.at(0));

	    is_zdr=true;
	    }
	    // Z_VD
	    if(i<full_volume_z.size()-1){
	      //radarelab::PolarScan<double> Z_VDi(full_volume_z.at(i));

	      //PolarScan<double> Z_VDi(full_volume_z.at(i)); //.beam_count, full_volume_z.at(i).beam_size, full_volume_z.at(i).undetect);
	      //cout<<"z undetect, gain="<<full_volume_z.at(i).undetect<<", "<<full_volume_z.at(i).gain<<endl;
	      //cout<<"beam_count, beam_size, nodata, undetect,gain for Z_VDi= "<<Z_VDi.beam_count<<", "<<Z_VDi.beam_size<<", "<<Z_VDi.nodata<<", "<<Z_VDi.undetect<<", "<<Z_VDi.gain<<endl;
	      //radarelab::volume::textureVD(full_volume_z.at(i), full_volume_z.at(i+1), Z_VDi, true);
	      volume::Scans<double> Input2; 
	      Input2.push_back(full_volume_z.at(i+1));
	      radarelab::volume::textureVD(Input, Input2, Texture, true);
	      Texture.at(0).nodata=65535.;
	      Texture.at(0).undetect=0.;
	      Z_VD.push_back(Texture.at(0));
	      cout<<"Z_VD.size="<<Z_VD.size()<<endl;
	    }
	}
      std::cout<<"Finito Cleaner, salvo risultati"<<std::endl;
	volume::ODIMStorer storer;
	storer.store_quantity_fp((Volume<double>*)(&Z_SD2D));
	storer.store_quantity_fp((Volume<double>*)(&Z_SDRay_9));
	storer.store_quantity_fp((Volume<double>*)(&Z_SDRay_21));
	storer.store_quantity_fp((Volume<double>*)(&Z_SDAz));
	if( is_zdr ) storer.store_quantity_fp((Volume<double>*)(&ZDR_SD2D));
	storer.store_quantity_fp((Volume<double>*)(&Z_VD));
	storer.store(argv[1]);
	cout<<endl<<"Fine"<<endl;

	//cout<<"vol00 = "<<Z_VD.at(5)(30,40)<<" da raw2="<<full_volume_z.at(5)(30,40)<<" - raw1"<<full_volume_z.at(5)(30,40)<<endl;
}


//
//    radarelab::volume::Scans<double>   Z_S,  SD2D,SD_Ray,SD_Az;
//    Z_S.push_back(scan_z);
//    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);

//    radarelab::volume::textureSD( Z_S,SD_Ray, scan_z.cell_size*9 , 360./scan_z.beam_count,true);
//    radarelab::volume::textureSD( Z_S,SD_Az, scan_z.cell_size , 5*360./scan_z.beam_count,true);
//
