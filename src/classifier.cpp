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

//#include "FIR_filter.h"
#include "classifier.h"
#include <radarlib/radar.hpp>

using namespace elaboradar;
using namespace volume;
using namespace std;
namespace odim = OdimH5v21;

PROB::PROB(double z,double zdr,double rhohv, double lkdp, double sdz, double sdphidp)
{
	this->resize(10,6);
	this->row(0)=prob_class(GC_AP,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(1)=prob_class(   BS,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(2)=prob_class(   DS,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(3)=prob_class(   WS,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(4)=prob_class(   CR,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(5)=prob_class(   GR,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(6)=prob_class(   BD,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(7)=prob_class(   RA,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(8)=prob_class(   HR,z,zdr,rhohv,lkdp,sdz,sdphidp);
	this->row(9)=prob_class(   RH,z,zdr,rhohv,lkdp,sdz,sdphidp);
}

Matrix2D<double> PROB::prob_class(EchoClass classe,double z, double zdr, double rhohv, double lkdp, double sdz, double sdphidp)
{
	Matrix2D<double> probability(1,6);
	probability=Matrix2D::Constant(1,6,0.);
	double f1=f_1(z);
	double f2=f_2(z);
	double f3=f_3(z);
	double g1=g_1(z);
	double g2=g_2(z);
	switch(classe)
	{
		case GC_AP:
			probability<<trap(15.,20.,70.,80.,z),trap(-4.,-2.,1.,2.,zdr),trap(0.5,0.6,0.9,0.95,rhohv),
			trap(-30.,-25.,10.,20.,lkdp),trap(2.,4.,10.,15.,sdz),trap(30.,40.,50.,60.,sdphidp);
			return probability;
		case BS:
			probability<<trap(5.,10.,20.,30.,z),trap(0.,2.,10.,12.,zdr),trap(0.3,0.5,0.8,0.83,rhohv),
			trap(-30.,-25.,10.,10.,lkdp),trap(1.,2.,4.,7.,sdz),trap(8.,10.,40.,60.,sdphidp);
			return probability;
		case DS:
			probability<<trap(5.,10.,35.,40.,z),trap(-0.3,0.,0.3,0.6,zdr),trap(0.95,0.98,1.,1.01,rhohv),
			trap(-30.,-25.,10.,20.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case WS:
			probability<<trap(25.,30.,40.,50.,z),trap(0.5,1.,2.,3.0,zdr),trap(0.88,0.92,0.95,0.985,rhohv),
			trap(-30.,-25.,10.,20.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case CR:
			probability<<trap(0.,5.,20.,25.,z),trap(0.1,0.4,3.,3.3,zdr),trap(0.95,0.98,1.,1.01,rhohv),
			trap(-5.,0.,10.,20.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case GR:
			probability<<trap(25.,35.,50.,55.,z),trap(-0.3,0,f1,f1+0.3,zdr),trap(0.9,0.97,1.,1.01,rhohv),
			trap(-30.,-25.,10.,20.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case BD:
			probability<<trap(20.,25.,45.,50.,z),trap(f2-0.3,f2,f3,f3+1.,zdr),trap(0.92,0.95,1.,1.01,rhohv),
			trap(g1-1.,g1,g2,g2+1.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case RA:
			probability<<trap(5.,10.,45.,50.,z),trap(f1-0.3,f1,f2,f2+0.5,zdr),trap(0.95,0.97,1.,1.01,rhohv),
			trap(g1-1.,g1,g2,g2+1.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case HR:
			probability<<trap(40.,45.,55.,60.,z),trap(f1-0.3,f1,f2,f2+0.5,zdr),trap(0.92,0.95,1.,1.01,rhohv),
			trap(g1-1.,g1,g2,g2+1.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		case RH:
			probability<<trap(45.,50.,75.,80.,z),trap(-0.3,0.,f1,f1+0.5,zdr),trap(0.85,0.9,1.,1.01,rhohv),
			trap(-10.,-4.,g1,g1+1.,lkdp),trap(0.,0.5,3.,6.,sdz),trap(0.,1.,15.,30.,sdphidp);
			return probability;
		default:
			cout<<"ERROR!!!   unknown echo type "<<classe<<endl;
			return probability;	// without it produces compile warnings and runtime error if reached
	}
}

double PROB::trap(double x1, double x2, double x3, double x4, double val)
{
	if(val<=x3&&val>=x2) return 1.;
	else if(val<x2&&val>x1) return val/(x2-x1)-x1/(x2-x1);
	else if (val<x4&&val>x3) return val/(x3-x4)-x4/(x3-x4);
	else return 0.; // (val<=x1||val>=x4)
}

HCA_Park::HCA_Park(double Z, double ZDR, double RHOHV, double LKDP, double SDZ, double SDPHIDP)
	: z(Z),zdr(ZDR),rhohv(RHOHV),lkdp(LKDP),sdz(SDZ),sdphidp(SDPHIDP)
{
	PROB Pij(z,zdr,rhohv,lkdp,sdz,sdphidp);
	CONF Qi;	// TODO: confidence vector calculation not implemented,
			// currently it uses a vector of ones.
	Matrix2D<double> Wij(10,6);

//		Z	Zdr	rhohv	lkdp	SDZ	SDphidp
	Wij <<	0.2,	0.4,	1.0,	0.0,	0.6,	0.8,	// GC_AP
		0.4,	0.6,	1.0,	0.0,	0.8,	0.8,	// BS
		1.0,	0.8,	0.6,	0.0,	0.2,	0.2,	// DS
		0.6,	0.8,	1.0,	0.0,	0.2,	0.2,	// WS
		1.0,	0.6,	0.4,	0.5,	0.2,	0.2,	// CR
		0.8,	1.0,	0.4,	0.0,	0.2,	0.2,	// GR
		0.8,	1.0,	0.6,	0.0,	0.2,	0.2,	// BD
		1.0,	0.8,	0.6,	0.0,	0.2,	0.2,	// RA
		1.0,	0.8,	0.6,	1.0,	0.2,	0.2,	// HR
		1.0,	0.8,	0.6,	1.0,	0.2,	0.2;	// RH
	Ai.resize(10);
	Ai=((Wij.array()*Pij.array()).matrix()*Qi).array()/(Wij*Qi).array();
}

classifier::classifier(const string& file, const Site& site):pathname(file)
{
	printf("il nome del mio file è %s\n", pathname.c_str());

	volume::ODIMLoader loader_all(site, false, 1024);

	volume::Scans<double> full_volume_z;
	volume::Scans<double> full_volume_zdr;
	volume::Scans<double> full_volume_rhohv;
	volume::Scans<double> full_volume_phidp;
	volume::Scans<double> full_volume_vrad;

	loader_all.request_quantity(odim::PRODUCT_QUANTITY_DBZH,&full_volume_z);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_ZDR,&full_volume_zdr);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_RHOHV,&full_volume_rhohv);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_PHIDP,&full_volume_phidp);
	loader_all.request_quantity(odim::PRODUCT_QUANTITY_VRAD,&full_volume_vrad);

	loader_all.load(pathname);
	printf("Non so se è andato tutto bene, ma almeno sono arrivato in fondo\n");

	volume::volume_resample<double>(full_volume_z, loader_all.azimuth_maps, vol_z, volume::merger_closest<double>);
	volume::volume_resample<double>(full_volume_zdr, loader_all.azimuth_maps, vol_zdr, volume::merger_closest<double>);
	volume::volume_resample<double>(full_volume_rhohv, loader_all.azimuth_maps, vol_rhohv, volume::merger_closest<double>);
	volume::volume_resample<double>(full_volume_phidp, loader_all.azimuth_maps, vol_phidp, volume::merger_closest<double>);
	volume::volume_resample<double>(full_volume_vrad, loader_all.azimuth_maps, vol_vrad, volume::merger_closest<double>);
}

void classifier::compute_lkdp()
{
	// TODO: la seguente è la traduzione del metodo per il calcolo di lkdp di Park et al. (2009)
	// capire se si vuole fare diverso

/*	// TODO: da reinserire negli opportuni polarscan
	vol_lkdp_2km.quantity.name="LKDP";
	vol_lkdp_2km.quantity.units="dB°/km";
	vol_lkdp_2km.quantity.nodata=-9999.;
	vol_lkdp_2km.quantity.undetect=-9999.;
	vol_lkdp_6km.quantity.name="LKDP";
	vol_lkdp_6km.quantity.units="dB°/km";
	vol_lkdp_6km.quantity.nodata=-9999.;
	vol_lkdp_6km.quantity.undetect=-9999.;
*/	
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
/*
	for(unsigned el=0; el<vol_phidp.size();el++)
	{
		cout<<"el= "<<el<<endl;
		double* re_in = new double[vol_phidp.scan(el).beam_size];
		FIR_filter fil(vol_phidp.scan(el).beam_size,re_in);
		for(unsigned az=0;az<vol_phidp.scan(el).beam_count;az++)
		{
			vol_phidp.scan(el).set(az,0,vol_phidp.scan(el).row(az).mean());
			for(unsigned rg=1;rg<vol_phidp.scan(el).beam_size;rg++)
			{
				if(vol_phidp.scan(el).get(az,rg)<-179.) vol_phidp.scan(el).set(az,rg,vol_phidp.scan(el).get(az,rg-1));
			}
			fil.feed(vol_phidp.scan(el).row_ptr(az));
			fil.perform();

			/*if(el==5&&az==65)
			{
				//for(unsigned rg=0;rg<vol_phidp.scan(el).beam_size;rg++) vol_phidp.scan(el).set(az,rg,rg==0?1:0);
				//fil.feed(vol_phidp.scan(el).row_ptr(az));
				//fil.perform();
				//fil.dump();
				cout<<endl<<endl;
				for(unsigned i=0;i<400;i++)
				{
					cout<<fixed<<vol_phidp.scan(el).get(az,i)<<"\t"<<re_in[i]/vol_phidp.scan(el).beam_size<<endl;
				}
				cout<<endl<<endl;
			}///
		}
		cout<<"fine el= "<<el<<endl;
		delete re_in;
	}
*/
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
	vector< vector< HCA_Park > > SCAN;
	vector< HCA_Park > BEAM;
	printf("inizio HCA\n");
	vol_Ai.resize(vol_z.size());
	double Z,Zdr,rhohv,lkdp,sdz,sdphidp;
	for(unsigned el=0;el<vol_z.size();el++)
	{
		cout<<"\tHCA el "<<el<<endl;
		SCAN.resize(vol_z.scan(el).beam_count);
		for(unsigned az=0;az<vol_z.scan(el).beam_count;az++)
		{
			BEAM.resize(vol_z.scan(el).beam_size);
			for(unsigned rg=0;rg<vol_z.scan(el).beam_size;rg++)
			{
				Z=vol_z_1km.scan(el).get(az,rg);
				Zdr=vol_zdr_2km.scan(el).get(az,rg);
				rhohv=vol_rhohv_2km.scan(el).get(az,rg);

				lkdp=Z>40?vol_lkdp_2km.scan(el).get(az,rg):vol_lkdp_6km.scan(el).get(az,rg);

				sdz=vol_sdz.scan(el).get(az,rg);
				sdphidp=vol_sdphidp.scan(el).get(az,rg);

				HCA_Park hca(Z,Zdr,rhohv,lkdp,sdz,sdphidp);
				BEAM[rg]=hca;
			}
			SCAN[az]=BEAM;
		}
		vol_Ai[el]=SCAN;
	}
	// Dopo aver calcolato i valori di aggregazione cerco il melting layer
	MeltingLayer ML(vol_z,vol_zdr,vol_rhohv,vol_Ai);
	cout<<"uscito da ML"<<endl;
	//TODO:check aggregation values
	//TODO:check hard thresholds
	unsigned elev=2;
	unsigned azim=75;
	cout<<"GC\tBS\tDS\tWS\tCR\tGR\tBD\tRA\tHR\tRH"<<endl;
	for(unsigned rg=0;rg<vol_Ai[elev][azim].size();rg++)
	{
		cout<<fixed<<vol_Ai[elev][azim][rg].Ai[GC_AP]<<"\t"<<vol_Ai[elev][azim][rg].Ai[BS]<<"\t"<<
		vol_Ai[elev][azim][rg].Ai[DS]<<"\t"<<vol_Ai[elev][azim][rg].Ai[WS]<<"\t"<<
		vol_Ai[elev][azim][rg].Ai[CR]<<"\t"<<vol_Ai[elev][azim][rg].Ai[GR]<<"\t"<<
		vol_Ai[elev][azim][rg].Ai[BD]<<"\t"<<vol_Ai[elev][azim][rg].Ai[RA]<<"\t"<<
		vol_Ai[elev][azim][rg].Ai[HR]<<"\t"<<vol_Ai[elev][azim][rg].Ai[RH]<<endl;
	}

	
}
