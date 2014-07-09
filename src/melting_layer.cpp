/*
 * =====================================================================================
 *
 *       Filename:  melting_layer.cpp
 *
 *    Description:  implementation of melting_layer.h methods
 *
 *        Version:  1.0
 *        Created:  9/07/2014 15:30
 *       Revision:  
 *       Compiler:  gcc
 *
 *         Author:  Davide Ori 
 *   Organization:  Dept.Physics University of Bologna
 *
 * =====================================================================================
 */


#include "melting_layer.h"

using namespace cumbac;
using namespace volume;
using namespace std;

double diff_height(PolarScan<double>& scan, unsigned rg_start, unsigned rg_end)
{
	//TODO
	return 1.;
}

void increment(MLpoints& matrix,PolarScan<double>& scan, unsigned az_idx, unsigned rg_idx)
{
	//TODO
	unsigned m_rg_idx=0;
	matrix(az_idx,m_rg_idx)++;
}

MeltingLayer::MeltingLayer(Volume<double>& vol_z,Volume<double>& vol_zdr,Volume<double>& vol_rhohv)
{
	vol_z_0_5km.filter(vol_z,500.);
	vol_zdr_1km.filter(vol_zdr,1000.);
	vol_rhohv_1km.filter(vol_rhohv,1000.);

	//TODO: correzione attenuazione phidp
	//TODO: altro preprocessing Ryzhkov 2005b

	MLpoints melting_points(1.0,10.,vol_z.beam_count,100);
	unsigned curr_rg=0;
	bool confirmed=false;
	
	for(unsigned el=0;el<vol_rhohv_1km.size();el++)
	{
		PolarScan<double>& rho=vol_rhohv.scan(el);
		if(rho.elevation>4.&&rho.elevation<10.)
		{
			PolarScan<double>& z=vol_z_0_5km.scan(el);
			PolarScan<double>& zdr=vol_zdr_1km.scan(el);
			for(unsigned rg=0;rg<rho.beam_size;rg++)	//TODO: check for climatological boundaries in ML height
				for(unsigned az=0;az<rho.beam_count;az++)
				{
					if(rho(az,rg)>=0.9&&rho(az,rg)<=0.97) // TODO: check also GC_AP & BS
					{
						curr_rg=rg;
						while(curr_rg<z.beam_size && diff_height(z,rg,curr_rg)<500 && !confirmed)
						{
							if(z(az,rg)>30 && z(az,rg)<47 && zdr(az,rg)>0.8 && zdr(az,rg)<2.5)
							{
								confirmed=true;
							}
							curr_rg++;
						}
						if(confirmed)
						{
							increment(melting_points,z,az,rg);
						}
					}
					confirmed=false;
				}
		}
	}
}
