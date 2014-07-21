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


//#include "melting_layer.h"
#include "classifier.h"

#define kea 8494666.666667	// c'è qualcosa in geo_par.h

using namespace cumbac;
using namespace volume;
using namespace std;

double height(PolarScan<double>& scan, unsigned rg)
{
	double range=(double)rg*scan.cell_size;
	double h=sqrt(range*range+kea*kea+2.*kea*range*sin(scan.elevation*M_PI/180.))-kea;	//meters
	return h/1000.;	//km
}

double diff_height(PolarScan<double>& scan, unsigned rg_start, unsigned rg_end)
{
	return fabs(height(scan, rg_end)-height(scan, rg_start));
}

void increment(MLpoints& matrix,PolarScan<double>& scan, unsigned az_idx, unsigned rg_idx)
{
	unsigned m_h_idx=matrix.h_idx(height(scan,rg_idx));
	unsigned m_az_idx=matrix.deg2idx((double)az_idx*360./scan.beam_count);
	matrix(m_h_idx,m_az_idx)++;
}

void MLpoints::box_top_bottom(double box_width_deg, double bot_th, double top_th, std::vector<double>& ML_b, std::vector<double>& ML_t)
{
	if(bot_th<0.||bot_th>1.) cout<<"ERROR bot_th must be 0<%<1 "<<endl;
	if(top_th<0.||top_th>1.) cout<<"ERROR top_th must be 0<%<1 "<<endl;
	if(top_th<bot_th) cout<<"ERROR top_th must be > than bot_th"<<endl;
	unsigned width=1+2*std::floor(0.5*box_width_deg*this->cols()/360.);
	unsigned half=0.5*(width-1);
	unsigned box_count=0;
	double bottom_lim;
	double top_lim;
	for(unsigned az=0;az<this->cols();az++)
	{
		ML_b[az]= -99.;
		ML_t[az]= -99.;
		box_count=0;
		unsigned round_bm;
		for(int bm=az-half;bm<az+half+1;bm++)
		{
			if(bm<0) round_bm=this->cols()-bm;
			else if(bm>=this->cols()) round_bm=bm-this->cols();
			else round_bm=bm;
			
			for(unsigned h=0;h<this->rows();h++)box_count+=(*this)(h,round_bm);
		}
		bottom_lim=bot_th*box_count;
		top_lim=top_th*box_count;
		if(box_count!=0)
		{
			box_count=0;
			for(unsigned h=0;h<this->rows();h++)
			{
				for(unsigned bm=az-half;bm<az+half+1;bm++)
				{
					if(bm<0) round_bm=this->cols()-bm;
					else if(bm>=this->cols()) round_bm=bm-this->cols();
					else round_bm=bm;

					box_count+=(*this)(h,round_bm);
				}
				if(ML_b[az]<0 && box_count>bottom_lim)ML_b[az]=this->Hmin+h*(this->Hmax-this->Hmin)/this->rows();
				if(ML_t[az]<0 && box_count>top_lim)   ML_t[az]=this->Hmin+h*(this->Hmax-this->Hmin)/this->rows();
			}
		}
	}
}

MeltingLayer::MeltingLayer(Volume<double>& vol_z,Volume<double>& vol_zdr,Volume<double>& vol_rhohv, vector< vector< vector< HCA_Park> > >& HCA)
{
	cout<<"\tInizio melting Layer"<<endl;
	vol_z_0_5km.filter(vol_z,500.);
	vol_zdr_1km.filter(vol_zdr,1000.);
	vol_rhohv_1km.filter(vol_rhohv,1000.);

	//TODO: correzione attenuazione con phidp
	//TODO: altro preprocessing Ryzhkov 2005b ??? sull'articolo non c'è nulla

	MLpoints melting_points(0.,10.,vol_z.beam_count,200);
	ML_top.resize(vol_z.beam_count);
	ML_bot.resize(vol_z.beam_count);
	unsigned curr_rg=0;
	bool confirmed=false;
	
	for(unsigned el=0;el<vol_rhohv_1km.size();el++)
	{
		cout<<"\t\t ML el "<<el<<endl;
		PolarScan<double>& rho=vol_rhohv.scan(el);
		if(rho.elevation>4.&&rho.elevation<10.)
		{
			PolarScan<double>& z=vol_z_0_5km.scan(el);
			PolarScan<double>& zdr=vol_zdr_1km.scan(el);
			for(unsigned rg=0;rg<rho.beam_size;rg++)	//TODO: check for climatological boundaries in ML height
				for(unsigned az=0;az<rho.beam_count;az++)
				{
					if(rho(az,rg)>=0.9 && rho(az,rg)<=0.97 && HCA[el][az][rg].meteo_echo())
					{
						curr_rg=rg;
						while(curr_rg<z.beam_size && diff_height(z,rg,curr_rg)<0.5 && !confirmed)
						{
							if(z(az,rg)>30 && z(az,rg)<47 && zdr(az,rg)>0.8 && zdr(az,rg)<2.5 && height(z,rg)>melting_points.Hmin)
							{
								confirmed=true;
							}
							curr_rg++;
						}
						if(confirmed)
						{
							increment(melting_points,z,az,rg);
							melting_points.count++;
							confirmed=false;
						}
					}
				}
		}
	}

	cout<<"I punti ML trovati sono "<<melting_points.count<<endl;
	melting_points.box_top_bottom(20.,0.2,0.8,ML_bot,ML_top);
	cout<<"Altezza ML"<<endl;
	for(unsigned i=0;i<ML_bot.size();i++)	cout<<ML_bot[i]<<"\t"<<ML_top[i]<<endl;
	//TODO: fill empty azimuths
}