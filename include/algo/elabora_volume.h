#include"volume.h"
#include"statistics.h"

namespace elaboradar {

template<typename T>
PolarScan<T> make_slope_scan(PolarScan<T>& raw, unsigned win)
{
	unsigned half_win=0.5*(win-1);
	PolarScan<T> scan(raw);
	Statistic<T> fit;
	for(unsigned i=0;i<raw.rows();i++)
	{
		fit.clear();
		for(unsigned j=0;j<half_win+1;j++)
			if(raw(i,j)!=raw.undetect&&raw(i,j)!=raw.nodata)
				fit.feed(0*raw.cell_size,raw(i,j));
		scan.set(i,0,fit.compute_slope());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);

		for(unsigned j=1;j<raw.beam_size+1;j++)
		{
			if(j-half_win-1>=0)
				if(raw(i,j-half_win-1)!=raw.undetect&&raw(i,j-half_win-1)!=raw.nodata)
					fit.slim((j-half_win-1)*raw.cell_size,raw(i,j-half_win-1));
			if(j+half_win<=raw.beam_size)
				if(raw(i,j+half_win)!=raw.undetect&&raw(i,j+half_win)!=raw.nodata)
					fit.feed(j*raw.cell_size,raw(i,j+half_win));
			scan.set(i,j,fit.compute_slope());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}

template<typename T>
void make_rms_scan(PolarScan<T>& raw, unsigned win)
{
	unsigned half_win=0.5*(win-1);
	PolarScan<T> scan(raw);
	Statistic<T> rms;
	for(unsigned i=0;i<raw.rows();i++)
	{
		rms.clear();
		for(unsigned j=0;j<half_win+1;j++)
			if(raw(i,j)!=raw.undetect&&raw(i,j)!=raw.nodata)
				rms.feed(raw(i,j));
		scan.set(i,0,rms.compute_dev_std());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);

		for(unsigned j=1;j<raw.beam_size+1;j++)
		{
			if(j-half_win-1>=0)
				if(raw(i,j-half_win-1)!=raw.undetect&&raw(i,j-half_win-1)!=raw.nodata)
					rms.slim(raw(i,j-half_win-1));
			if(j+half_win<=raw.beam_size)
				if(raw(i,j+half_win)!=raw.undetect&&raw(i,j+half_win)!=raw.nodata)
					rms.feed(raw(i,j+half_win));
			scan.set(i,j,rms.compute_dev_std());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}


namespace volume {

template<typename T>
Volume<T> moving_average_slope(Volume<T>& raw, double slope_range) // least-squares
{
	unsigned window_size;
	Volume<T> vol;
	//this->quantity=raw.quantity.quantity_slope(); // TODO: è complesso ma si potrebbe
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*slope_range/raw.scan(i).cell_size);
		vol.push_back(make_slope_scan(raw.scan(i),window_size));
	}
	return vol;
}

template<typename T>
Volume<T> textureSD(Volume<T>& raw, double filter_range) // least-squares
{
	unsigned window_size;
	Volume<T> vol;
	//this->quantity=raw.quantity.quantity_slope(); // TODO: è complesso ma si potrebbe 

	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		vol.push_back(make_rms_scan(raw.scan(i),window_size));
	}
	return vol;
}

}	// namespace volume
}	// namespace elaboradar
