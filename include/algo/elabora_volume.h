#include"volume.h"
#include"statistics.h"

using namespace std;

namespace elaboradar {

namespace {bool check_undetect=false;}

template<typename T>
inline bool good(const PolarScan<T>& scan, unsigned row, unsigned col)
{
	if(check_undetect)
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
			pre=j-half_win;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,pre))
					fit.slim(pre*raw.cell_size,raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,post))
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
			pre=j-half_win;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,pre))
					rms.slim(raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,post))
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
			pre=j-half_win;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,pre))
					filter.slim(raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,post))
					filter.feed(raw(i,post));
			scan.set(i,j,filter.compute_mean());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}

namespace volume {

template<typename T>
void moving_average_slope(const Volume<T>& raw, Volume<T>& vol, double slope_range, bool force_check_undetect=false)
{
	unsigned window_size;
	vol.clear();
	check_undetect=force_check_undetect;
	//this->quantity=raw.quantity.quantity_slope(); // TODO: è complesso ma si potrebbe
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*slope_range/raw.scan(i).cell_size);
		vol.push_back(make_slope_scan(raw.scan(i),window_size));
	}
}

template<typename T>
void textureSD(const Volume<T>& raw, Volume<T>& vol, double filter_range, bool force_check_undetect=false)
{
	unsigned window_size;
	vol.clear();
	check_undetect=force_check_undetect;
	//this->quantity=raw.quantity.quantity_rms(); // TODO: è complesso ma si potrebbe 
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		vol.push_back(make_rms_scan(raw.scan(i),window_size));
	}
}

template<typename T>
void filter(const Volume<T>& raw, Volume<T>& vol, double filter_range, bool force_check_undetect=false)
{
	unsigned window_size;
	vol.clear();
	check_undetect=force_check_undetect;
	//this->quantity=raw.quantity.quantity_slope(); // TODO: è complesso ma si potrebbe 
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		vol.push_back(make_filter_scan(raw.scan(i),window_size));
	}
}

template<typename T>
void lin2dB(Volume<T>& lin, Volume<T>& dB)
{
	//this->quantity=lin.quantity.lin2dB(); // TODO: not yet implemented
	dB.clear();
	for(unsigned i=0;i<lin.size();i++)
	{
		dB.push_back(PolarScan<T>(lin.scan(i).beam_count,lin.scan(i).beam_size,0.));
		dB[i].cell_size = lin[i].cell_size;
		dB[i].elevation = lin.scan(i).elevation;
		dB[i].block(0,0,lin.scan(i).beam_count,lin.scan(i).beam_size)=lin.scan(i).log10();		
		dB[i].array()*=10.;
	}
}

template<typename T>
void dB2lin(Volume<T>& dB, Volume<T>& lin)
{
	//this->quantity=DB.quantity.dB2lin(); // TODO: not yet implemented
	lin.clear();
	for(unsigned i=0;i<dB.size();i++)
	{
		lin.push_back(PolarScan<T>(dB.scan(i).beam_count,dB.scan(i).beam_size,0.));
		lin[i].cell_size = dB[i].cell_size;
		lin[i].elevation = dB.scan(i).elevation;
		lin[i].block(0,0,dB.scan(i).beam_count,dB.scan(i).beam_size) = dB.scan(i)*0.1;
		lin[i].block(0,0,dB.scan(i).beam_count,dB.scan(i).beam_size) = lin.scan(i).exp10();
	}
}

template<typename T>
void lin2dB(Volume<T>& lin)
{
	//this->quantity=lin.quantity.lin2dB(); // TODO: not yet implemented
	for(unsigned i=0;i<lin.size();i++)
		lin[i].block(0,0,lin.scan(i).beam_count,lin.scan(i).beam_size)=lin.scan(i).log10();
	lin*=10.;
}

template<typename T>
void dB2lin(Volume<T>& dB)
{
	//this->quantity=DB.quantity.dB2lin(); // TODO: not yet implemented
	dB*=0.1;
	for(unsigned i=0;i<dB.size();i++)
		dB[i].block(0,0,dB.scan(i).beam_count,dB.scan(i).beam_size) = dB.scan(i).exp10();
}

}	// namespace volume
}	// namespace elaboradar
