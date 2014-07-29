#include"volume.h"
#include"statistics.h"

namespace elaboradar {

template<typename T>
PolarScan<T> make_slope_scan_range(PolarScan<T>& raw, unsigned win)
{
	unsigned half_win=0.5*(win-1);
	PolarScan<T> scan(raw);
	//this->push_back(PolarScan<T>(raw.beam_count,raw.beam_size,0.));
	//this->back().elevation = raw.elevation;
	//this->back().cell_size = raw.cell_size;

	//	using namespace stat;
	LinearFit<T> fit;
	for(unsigned i=0;i<raw.rows();i++)
	{
		for(unsigned j=0;j<half_win;j++)
		{
			for(unsigned k=0;k<(half_win+j+1);k++)
			{
				if((raw(i,k)!=raw.undetect)&&(raw(i,k)!=raw.nodata))
				{
					fit.feed(k*raw.cell_size,raw.get(i,k));
				}
			}
			if(fit.N)
			{
				scan.set(i,j,fit.compute_slope());
			}
			else scan.set(i,j,raw.nodata);
			fit.clear();

			for(unsigned k=0;k<(win-j-1);k++)
			{
				if((raw(i,raw.beam_size-win+1+j+k)!=raw.undetect)&&(raw(i,raw.beam_size-win+1+j+k)!=raw.nodata))
				{
					fit.feed(k*raw.cell_size,raw.get(i,raw.beam_size-win+1+j+k));
				}
			}
			if(fit.N)
			{
				scan.set(i,raw.beam_size-half_win+j,fit.compute_slope());
			}
			else scan.set(i,raw.beam_size-half_win+j,raw.nodata);
			fit.clear();
		}
		for(unsigned j=half_win;j<(raw.beam_size-half_win);j++)
		{
			for(unsigned k=0;k<win;k++)
			{
				if((raw(i,j-half_win+k)!=raw.undetect)&&(raw(i,j-half_win+k)!=raw.nodata))
				{
					fit.feed(k*raw.cell_size,raw.get(i,j-half_win+k));
				}
			}
			if(fit.N)
			{
				scan.set(i,j,fit.compute_slope());
			}
			else scan.set(i,j,raw.nodata);
			fit.clear();
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
							// intervenire con un metodo che 
							// determina la quantity della slope 
							// in funzione della quantity di raw.
							// Per adesso si suppone che la quantity
							// del volume di slope sia settata a priori 
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*slope_range/raw.scan(i).cell_size);
		vol.push_back(make_slope_scan_range(raw.scan(i),window_size));
	}
	return vol;
}



}	// namespace volume
}	// namespace elaboradar
