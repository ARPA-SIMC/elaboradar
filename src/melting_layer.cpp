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

#include "classifier.h"
#include "algo/elabora_volume.h"

using namespace elaboradar;
using namespace volume;
using namespace std;

void increment(MLpoints& matrix,PolarScan<double>& scan, unsigned az_idx, unsigned rg_idx)
{
	unsigned m_h_idx=matrix.h_idx(scan.height(rg_idx));
	unsigned m_az_idx=matrix.deg2idx((double)az_idx*360./scan.beam_count);
	matrix(m_h_idx,m_az_idx)++;
	matrix.count++;
}

void MLpoints::box_top_bottom(double box_width_deg, double bot_th, double top_th, std::vector<double>& ML_b, std::vector<double>& ML_t)
{
	if(bot_th<0.||bot_th>1.) cout<<"ERROR bot_th must be 0<%<1 "<<endl;
	if(top_th<0.||top_th>1.) cout<<"ERROR top_th must be 0<%<1 "<<endl;
	if(top_th<bot_th) cout<<"ERROR top_th must be > than bot_th"<<endl;
	int width=1+2*std::floor(0.5*box_width_deg*this->cols()/360.);
	int half=0.5*(width-1);
	unsigned box_count=0;
	double bottom_lim;
	double top_lim;
	for(unsigned az=0;az<this->cols();az++)
	{
		ML_b[az]= -99.;
		ML_t[az]= -99.;
		box_count=0;
		unsigned round_bm;
		width=az+half+1;
		//cout<<az<<"\t";
		for(int bm=az-half;bm<width;bm++)
		{
			if(bm<0) round_bm=this->cols()+bm;
			else if(bm>=this->cols()) round_bm=bm-this->cols();
			else round_bm=bm;
			for(unsigned h=0;h<this->rows();h++)box_count+=(*this)(h,round_bm);
		}
		bottom_lim=bot_th*box_count;
		top_lim=top_th*box_count;
		//cout<<box_count<<endl;
		if(box_count>=100)
		{
			box_count=0;
			for(unsigned h=0;h<this->rows();h++)
			{
				for(int bm=az-half;bm<width;bm++)
				{
					if(bm<0) round_bm=this->cols()+bm;
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

void MeltingLayer::seek4mlfile(time_t now, MLpoints& mlp)
{
	unsigned max_time=16*60;
	ostringstream filename;
	fstream file;
	unsigned m_h_idx;
	unsigned az, app;
	double h;
	while(max_time)
	{
		max_time-=60;
		now-=60;
		filename<<now<<".ml";
		file.open(filename.str());
		if(file.is_open())
		{
			cout<<endl<<"Apro MLfile "<<filename.str()<<endl;
			while(file)
			{
				file>>app;
				if(file.eof()) break;
				file>>az>>h;
				m_h_idx=mlp.h_idx(h);
				mlp(m_h_idx,az)++;
				mlp.count++;
			}
			file.close();
		}
		filename.seekp(0);
	}
}

void MeltingLayer::fill_empty_azimuths()
{
	// X FACILE riempo con la media
	// DIFFICILE interpolo fra settori buoni
	// TODO se tutti i punti sono -99 la media è 0 di default!! verifica a posteriori la numerosità del campione della media
	Statistic<double> mean_bot;
	Statistic<double> mean_top;
	for(unsigned i=0;i<bot.size();i++)
	{
		if(bot[i]!=-99) mean_bot.feed(bot[i]);
		if(top[i]!=-99) mean_top.feed(top[i]);
	}
	mean_bot.compute_mean();
	mean_top.compute_mean();

	if(mean_bot.N==0)
	{
		mean_bot.mean=3.0;	//TODO: i valori di default dovrebbero essere noti da metodi scientifici (modello, T al suolo, ecc...)
		mean_top.mean=3.5;
	}
	for(unsigned i=0;i<bot.size();i++)
	{
		if(bot[i]==-99) bot[i]=mean_bot.mean;
		if(top[i]==-99) top[i]=mean_top.mean;
	}	
}

MeltingLayer::MeltingLayer(Volume<double>& vol_z,Volume<double>& vol_zdr,Volume<double>& vol_rhohv, 
								vector< vector< vector< HCA_Park> > >& HCA)
{
	cout<<"\tInizio melting Layer"<<endl;
//	filter(vol_z,vol_z_0_5km,1000.,0.,false);	// TODO: se tengo questo range di filtro, semplificare la struttura e riusare i vol_1km vol_2km già filtrati
//	filter(vol_zdr,vol_zdr_1km,2000.,0.,false);
//	filter(vol_rhohv,vol_rhohv_1km,2000.,0.,false);

	Volume<double> dummy;
	filter(vol_z,vol_z_0_5km,500.,3.,false);
	filter(vol_zdr,dummy,1000.,3.,false);
	filter(dummy,vol_zdr_1km,1000.,3.,false);
	filter(vol_rhohv,dummy,1000.,3.,false);
	filter(dummy,vol_rhohv_1km,1000.,3.,false);

	cout<<"filtrati"<<endl;
	//correzione attenuazione con phidp fatta a priori da chi ha invocato
	//altro preprocessing Ryzhkov 2005b ??? sull'articolo non c'è nulla

	double MAX_ML_H=4.5;	//TODO: check for climatological boundaries in ML height
	double MIN_ML_H=0.;

	MLpoints melting_points(MIN_ML_H,MAX_ML_H,vol_z.beam_count,100);
	top.resize(vol_z.beam_count);
	bot.resize(vol_z.beam_count);
	unsigned curr_rg=0;
	bool confirmed=false;

	ostringstream ML;
	ostringstream DATE;
	DATE<<vol_z.load_info->acq_date<<".ml";
	ofstream OUT;
	OUT.open(DATE.str());
	
	for(unsigned el=0;el<vol_rhohv_1km.size();el++)
	{
		cout<<"\t\t ML el ";
		PolarScan<double>& rho=vol_rhohv_1km.scan(el);
		if(rho.elevation>3.&&rho.elevation<10.)	//TODO abbasso la soglia minima di un'elevazione
		{
			cout<<el<<endl;
			PolarScan<double>& z=vol_z_0_5km.scan(el);
			PolarScan<double>& zdr=vol_zdr_1km.scan(el);
			for(unsigned az=0;az<rho.beam_count;az++)
			{
				Statistic<double> average;
				double min_zdr=0.8;
				for(unsigned rg=0;rg<rho.beam_size;rg++)
				{
					if(zdr(az,rg)>0.2)average.feed(zdr(az,rg));
				}
				if(average.N) min_zdr=average.compute_mean();
				else min_zdr=0.8;
				
				for(unsigned rg=0;rg<rho.beam_size;rg++)
				{
					//if(el==5)cout<<rg<<" "<<rho.beam_size<<"\t"<<az<<" "<<rho.beam_count<<endl;
					if(rho(az,rg)>=0.9 && rho(az,rg)<=0.95 && HCA[el][az][rg].meteo_echo() && z.height(rg)>MIN_ML_H && z.height(rg)<MAX_ML_H)	//TODO diminuisco la soglia minima di rho da 0.9 a 0.85 e la massima da 0.97 a 0.95
					{
						curr_rg=rg;
						while(curr_rg<z.beam_size && z.diff_height(rg,curr_rg)<0.5 && !confirmed)
						{
							//if(el==5&&az==165&&rg==448)cout<<curr_rg<<endl;
							//if(el==4)if(az==85)if(rg>200)if(rg<250)cout<<rg<<" "<<rho(az,rg)<<" "<<z(az,rg)<<" "<<zdr(az,rg)<<endl;
							if(z(az,curr_rg)>20. && z(az,curr_rg)<47 && zdr(az,curr_rg)>min_zdr &&  //TODO cambio la soglia minima da 0.8 a metodo wise
								zdr(az,curr_rg)<2.5 && z.height(rg)>melting_points.Hmin)	// TODO diminuisco soglia minima z da 30 a 20
							{
								confirmed=true;
							}
							curr_rg++;
						}
						if(confirmed)
						{
							increment(melting_points,z,az,rg);
							confirmed=false;
							ML<<el<<"\t"<<az<<"\t"<<z.height(rg)<<endl;
						}
					}
				}
			}
		}
	}
	OUT<<ML.str();
	OUT.close();
	cout<<endl<<"I punti ML trovati sono "<<melting_points.count<<endl;
	cout<<"Cerco altri file di ML"<<endl;
	seek4mlfile(vol_z.load_info->acq_date, melting_points);
	cout<<"Ora i punti ML sono "<<melting_points.count<<endl;
	melting_points.box_top_bottom(20.,0.2,0.8,bot,top);
	fill_empty_azimuths();
	
//	cout<<"Altezza ML"<<endl;
//	for(unsigned i=0;i<bot.size();i++)cout<<bot[i]<<"\t"<<top[i]<<endl;
}
