#include"volume.h"
#include"statistics.h"

using namespace std;

namespace elaboradar {

namespace {bool phidp=false;}

template<typename T>
inline bool good(const PolarScan<T>& scan, unsigned row, unsigned col)
{
	if(phidp)
		if(scan(row,col)!=scan.nodata&&scan(row,col)!=scan.undetect)
			return true;
		else return false;
	else
		if(scan(row,col)!=scan.nodata)
			return true;
		else return false;
}

template<typename T>
PolarScan<T> make_slope_scan(const PolarScan<T>& raw, unsigned win)
{
	unsigned half_win=0.5*(win-1);
	PolarScan<T> scan(raw);
	Statistic<T> fit;
	int pre,post;
	for(unsigned i=0;i<raw.rows();i++)
	{
		fit.clear();
		for(unsigned j=0;j<half_win+1;j++)
			if(good(raw,i,j))
				fit.feed(0*raw.cell_size,raw(i,j));
		scan.set(i,0,fit.compute_slope());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);

		for(unsigned j=1;j<raw.beam_size;j++)
		{
			pre=j-half_win-1;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,j))
					fit.slim(pre*raw.cell_size,raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,j))
					fit.feed(post*raw.cell_size,raw(i,post));
			scan.set(i,j,fit.compute_slope());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}

template<typename T>
PolarScan<T> make_rms_scan(const PolarScan<T>& raw, unsigned win)
{
	unsigned half_win=0.5*(win-1);
	PolarScan<T> scan(raw);
	Statistic<T> rms;
	int pre,post;
	for(unsigned i=0;i<raw.rows();i++)
	{
		rms.clear();
		for(unsigned j=0;j<half_win+1;j++)
			if(good(raw,i,j))
				rms.feed(raw(i,j));
		scan.set(i,0,rms.compute_dev_std());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);

		for(unsigned j=1;j<raw.beam_size;j++)
		{
			pre=j-half_win-1;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,j))
					rms.slim(raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,j))
					rms.feed(raw(i,post));
			scan.set(i,j,rms.compute_dev_std());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}

template<typename T>
PolarScan<T> make_filter_scan(const PolarScan<T>& raw, unsigned win)
{
	unsigned half_win=0.5*(win-1);
	PolarScan<T> scan(raw);
	Statistic<T> filter;
	int pre,post;
	for(unsigned i=0;i<raw.rows();i++)
	{
		filter.clear();
		for(unsigned j=0;j<half_win+1;j++)
			if(good(raw,i,j))
				filter.feed(raw(i,j));
		scan.set(i,0,filter.compute_mean());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);
		for(unsigned j=1;j<raw.beam_size;j++)
		{
			pre=j-half_win-1;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,j))
					filter.slim(raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,j))
					filter.feed(raw(i,post));
			scan.set(i,j,filter.compute_mean());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}

namespace volume {

template<typename T>
void moving_average_slope(const Volume<T>& raw, Volume<T>& vol, double slope_range)
{
	unsigned window_size;
	vol.clear();
	if(strcmp(raw.quantity.c_str(),"PHIDP")) phidp=false;
	else phidp=true;
	//this->quantity=raw.quantity.quantity_slope(); // TODO: è complesso ma si potrebbe
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*slope_range/raw.scan(i).cell_size);
		vol.push_back(make_slope_scan(raw.scan(i),window_size));
	}
}

template<typename T>
void textureSD(const Volume<T>& raw, Volume<T>& vol, double filter_range)
{
	unsigned window_size;
	vol.clear();
	if(strcmp(raw.quantity.c_str(),"PHIDP")) phidp=false;
	else phidp=true;
	//this->quantity=raw.quantity.quantity_rms(); // TODO: è complesso ma si potrebbe 
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		vol.push_back(make_rms_scan(raw.scan(i),window_size));
	}
}

template<typename T>
void filter(const Volume<T>& raw, Volume<T>& vol, double filter_range)
{
	unsigned window_size;
	vol.clear();
	if(strcmp(raw.quantity.c_str(),"PHIDP")) phidp=false;
	else phidp=true;
	//this->quantity=raw.quantity.quantity_slope(); // TODO: è complesso ma si potrebbe 
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		vol.push_back(make_filter_scan(raw.scan(i),window_size));
	}
}

}	// namespace volume
}	// namespace elaboradar
