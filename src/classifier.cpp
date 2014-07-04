/*
 * =====================================================================================
 *
 *       Filename:  classifier.cpp
 *
 *    Description:  implementation of classifier.h methods
 *
 *        Version:  1.0
 *        Created:  23/06/2014 15:30
 *       Revision:  
 *       Compiler:  gcc
 *
 *         Author:  Davide Ori 
 *   Organization:  Dept.Physics University of Bologna
 *
 * =====================================================================================
 */


#include "classifier.h"

using namespace cumbac;
using namespace volume;
using namespace std;

PROB::PROB(double z,double zdr,double rhohv, double phidp, double sdz, double sdphidp)
{
	printf("cicci\n");
}

HCA_Park::HCA_Park(double Z, double ZDR, double RHOHV, double LKDP, double SDZ, double SDPHIDP)
	: z(Z),zdr(ZDR),rhohv(RHOHV),lkdp(LKDP),sdz(SDZ),sdphidp(SDPHIDP)
{
	PROB Pij(z,zdr,rhohv,lkdp,sdz,sdphidp);
//	CONF Qi;
}

classifier::classifier(const string& file, const Site& site):pathname(file)
{
	printf("il nome del mio file è %s\n", pathname.c_str());

	// TODO: I think that only 1 loader is needed and could be reused to load all variables
	// Just reset load_info for every variable to ensure coherent loading (same beams)
	volume::ODIMLoader loader_z(site, false, false, 1024);
	volume::ODIMLoader loader_zdr(site, false, false, 1024);
	volume::ODIMLoader loader_rhohv(site, false, false, 1024);
	volume::ODIMLoader loader_phidp(site, false, false, 1024);
	volume::ODIMLoader loader_vrad(site, false, false, 1024);
	loader_z.coherent_loader=true;
	loader_zdr.coherent_loader=true;
	loader_rhohv.coherent_loader=true;
	loader_phidp.coherent_loader=true;
	loader_vrad.coherent_loader=true;

	bool file_ok=true;

	volume::LoadInfo load_info_z;
	volume::LoadInfo load_info_zdr;
	volume::LoadInfo load_info_rhohv;
	volume::LoadInfo load_info_phidp;
	volume::LoadInfo load_info_vrad;

	loader_z.load_info = &load_info_z;
	loader_z.vol_z = &vol_z;
	file_ok = file_ok && loader_z.load(pathname,"DBZH");

	loader_zdr.load_info = &load_info_zdr;
	loader_zdr.vol_z = &vol_zdr;
	file_ok = file_ok && loader_zdr.load(pathname,"ZDR");

	loader_rhohv.load_info = &load_info_rhohv;
	loader_rhohv.vol_z = &vol_rhohv;
	file_ok = file_ok && loader_rhohv.load(pathname,"RHOHV");
	
	loader_phidp.load_info = &load_info_phidp;
	loader_phidp.vol_z = &vol_phidp;
	file_ok = file_ok && loader_phidp.load(pathname,"PHIDP");

	loader_vrad.load_info = &load_info_vrad;
	loader_vrad.vol_z = &vol_vrad;
	file_ok = file_ok && loader_vrad.load(pathname,"VRAD");

	if(file_ok) printf("Everything is fine\n");
	else printf("Something went wrong during loading %s\n",pathname.c_str());
}

void classifier::compute_lkdp()
{
	// TODO: la seguente è la traduzione del metodo per il calcolo di lkdp di Park et al. (2009)
	// capire se si vuole fare diverso

	vol_lkdp_2km.quantity.name="LKDP";
	vol_lkdp_2km.quantity.units="dB°/km";
	vol_lkdp_2km.quantity.nodata=-9999.;
	vol_lkdp_2km.quantity.undetect=-9999.;
	vol_lkdp_6km.quantity.name="LKDP";
	vol_lkdp_6km.quantity.units="dB°/km";
	vol_lkdp_6km.quantity.nodata=-9999.;
	vol_lkdp_6km.quantity.undetect=-9999.;
	
	printf("calcolo kdp 2km\n");
	vol_lkdp_2km.moving_average_slope(vol_phidp_2km,2000.);
	vol_lkdp_2km*=1000.;
	printf("calcolo kdp 6 km\n");
	vol_lkdp_6km.moving_average_slope(vol_phidp_6km,6000.);
	vol_lkdp_6km*=1000.;

	double lkdp=0;
	for(unsigned el=0; el<vol_phidp.size();el++)
	{
		PolarScan<double>& lkdp2 = vol_lkdp_2km.scan(el);
		PolarScan<double>& lkdp6 = vol_lkdp_6km.scan(el);
		for(unsigned az=0; az<lkdp2.beam_count;az++)
		{
			for(unsigned rg=0; rg<lkdp2.beam_size;rg++)
			{
				lkdp=lkdp2.get(az,rg);
				lkdp2.set(az,rg,lkdp>0.001?10.*std::log10(lkdp):-30.);
				lkdp=lkdp6.get(az,rg);
				lkdp6.set(az,rg,lkdp>0.001?10.*std::log10(lkdp):-30.);
			}
		}
	}
}


void classifier::correct_phidp()
{
	// It is assumed that vol_phidp exist and has been initialized
	// moving window average over a range length of 2km and 6km as prescribed by Park et al. (2009)
	printf("filtro phidp 2 km\n");
	vol_phidp_2km.filter(vol_phidp,2000.);
	printf("filtro phidp 6 km\n");
	vol_phidp_6km.filter(vol_phidp,6000.);
}


void classifier::compute_derived_volumes()
{
	correct_phidp();

	// filtro i volumi
	printf("filtro Z 1 km\n");
	vol_z_1km.filter(vol_z,1000.);
	printf("filtro Zdr 2 km\n");
	vol_zdr_2km.filter(vol_zdr,2000.);
	printf("filtro rhohv 2 km\n");
	vol_rhohv_2km.filter(vol_rhohv,2000.);

	// calcolo le texture
	vol_sdz.textureSD(vol_z,1000.);
	vol_sdphidp.textureSD(vol_phidp,2000.);

/*
 * TODO: Now we should correct Z and Zdr for attenuation by adding estimated bias
 * dZ=0.04*vol_phidp_6km;
 * dZdr=0.004*vol_phidp_6km;
 */

	// calcolo lkdp
	compute_lkdp();
}

void classifier::HCA_Park_2009()
{
	printf("inizio HCA\n");
	for(unsigned el=0;el<vol_z.size();el++)
	{
		//HCA_Park hca(1.,2.,3.,4.,5.,6.);
	}
}
