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

//===========================PolarScan manipulators==========================//

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
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);	// changed to undetect

		for(unsigned j=1;j<raw.beam_size;j++)
		{
			pre=j-half_win-1;
			post=j+half_win;
			if(pre>=0)
				if(good(raw,i,pre))
					fit.slim(pre*raw.cell_size,raw(i,pre));
			if(post<raw.beam_size)
				if(good(raw,i,post))
					fit.feed(post*raw.cell_size,raw(i,post));
			scan.set(i,j,fit.compute_slope());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);	// changed to undetect
		}
	}
	return scan;
}

template<typename T>
PolarScan<T> make_rms_scan(const PolarScan<T>& raw, unsigned len, unsigned wid)
{
	unsigned half_len=0.5*(len-1);
	unsigned half_wid=0.5*(wid-1);
	PolarScan<T> scan(raw);
	Statistic<T> rms;
	int pre_l,post_l;
	int pre_w,post_w;
	for(unsigned i=0;i<raw.rows();i++)
	{
		rms.clear();
		for(unsigned j=0;j<half_len+1;j++)
		{
			if(good(raw,i,j)) rms.feed(raw(i,j));
			for(unsigned k=1;k<=half_wid;k++) // 0 è già inclusa, se half_wid=0 non dovrebbe partire il for
			{
				pre_w=i-k;
				post_w=i+k;
				if(pre_w<0) 
					pre_w=raw.beam_count+pre_w;
				if(post_w>=raw.beam_count)
					post_w=post_w-raw.beam_count;
				if(good(raw,pre_w,j))
					rms.feed(raw(pre_w,j));
				if(good(raw,post_w,j))
					rms.feed(raw(post_w,j));
			}
		}
		scan.set(i,0,rms.compute_dev_std());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);	// changed to undetect
		for(unsigned j=1;j<raw.beam_size;j++)
		{
			pre_l=j-half_len-1;
			post_l=j+half_len;
			if(pre_l>=0)
			{
				if(good(raw,i,pre_l))
					rms.slim(raw(i,pre_l));
				for(unsigned k=1;k<=half_wid;k++)
				{
					pre_w=i-k;
					post_w=i+k;
					if(pre_w<0) 
						pre_w=raw.beam_count+pre_w;
					if(post_w>=raw.beam_count)
						post_w=post_w-raw.beam_count;
					if(good(raw,pre_w,pre_l))
						rms.slim(raw(pre_w,pre_l));
					if(good(raw,post_w,pre_l))
						rms.slim(raw(post_w,pre_l));
				}
			}
			if(post_l<raw.beam_size)
			{
				if(good(raw,i,post_l))
					rms.feed(raw(i,post_l));
				for(unsigned k=1;k<=half_wid;k++)
				{
					pre_w=i-k;
					post_w=i+k;
					if(pre_w<0) 
						pre_w=raw.beam_count+pre_w;
					if(post_w>=raw.beam_count)
						post_w=post_w-raw.beam_count;
					if(good(raw,pre_w,post_l))
						rms.feed(raw(pre_w,post_l));
					if(good(raw,post_w,post_l))
						rms.feed(raw(post_w,post_l));
				}
			}
			scan.set(i,j,rms.compute_dev_std());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}

template<typename T>
PolarScan<T> make_filter_scan(const PolarScan<T>& raw, unsigned len, unsigned wid)
{
	unsigned half_len=0.5*(len-1);
	unsigned half_wid=0.5*(wid-1);
	PolarScan<T> scan(raw);
	Statistic<T> filter;
	int pre_l,post_l;
	int pre_w,post_w;
	for(unsigned i=0;i<raw.rows();i++)
	{
		filter.clear();
		for(unsigned j=0;j<half_len+1;j++)
		{
			if(good(raw,i,j)) filter.feed(raw(i,j));
			for(unsigned k=1;k<=half_wid;k++)
			{
				pre_w=i-k;
				post_w=i+k;
				if(pre_w<0) 
					pre_w=raw.beam_count+pre_w;
				if(post_w>=raw.beam_count)
					post_w=post_w-raw.beam_count;
				if(good(raw,pre_w,j))
					filter.feed(raw(pre_w,j));
				if(good(raw,post_w,j))
					filter.feed(raw(post_w,j));
			}
		}
		scan.set(i,0,filter.compute_mean());
		if(scan(i,0)!=scan(i,0)) scan.set(i,0,raw.nodata);	// change to undetect
		for(unsigned j=1;j<raw.beam_size;j++)
		{
			pre_l=j-half_len-1;
			post_l=j+half_len;
			if(pre_l>=0)
			{
				if(good(raw,i,pre_l))
					filter.slim(raw(i,pre_l));
				for(unsigned k=1;k<=half_wid;k++)
				{
					pre_w=i-k;
					post_w=i+k;
					if(pre_w<0) 
						pre_w=raw.beam_count+pre_w;
					if(post_w>=raw.beam_count)
						post_w=post_w-raw.beam_count;
					if(good(raw,pre_w,pre_l))
						filter.slim(raw(pre_w,pre_l));
					if(good(raw,post_w,pre_l))
						filter.slim(raw(post_w,pre_l));
				}
			}
			if(post_l<raw.beam_size)
			{
				if(good(raw,i,post_l))
					filter.feed(raw(i,post_l));
				for(unsigned k=1;k<=half_wid;k++)
				{
					pre_w=i-k;
					post_w=i+k;
					if(pre_w<0) 
						pre_w=raw.beam_count+pre_w;
					if(post_w>=raw.beam_count)
						post_w=post_w-raw.beam_count;
					if(good(raw,pre_w,post_l))
						filter.feed(raw(pre_w,post_l));
					if(good(raw,post_w,post_l))
						filter.feed(raw(post_w,post_l));
				}
			}
			scan.set(i,j,filter.compute_mean());
			if(scan(i,j)!=scan(i,j)) scan.set(i,j,raw.nodata);
		}
	}
	return scan;
}
template<typename T>
PolarScan<T> make_gradient_azimuth_scan(const PolarScan<T>& raw)
{
	PolarScan<T> scan(raw);

	for(unsigned rg=0;rg<raw.beam_size;rg++)
	{
		//cout<<"0"<<" "<<rg<<endl;
		if(good(raw,raw.beam_count-1,rg) && good(raw,0,rg))
			scan(0,rg) = raw(raw.beam_count-1,rg)-raw(0,rg);
		else scan(0,rg) = 0.;
	}
	for(unsigned az=1;az<raw.beam_count;az++)		
		for(unsigned rg=0;rg<raw.beam_size;rg++)
		{
			//cout<<az<<" "<<rg<<endl;
			if(good(raw,az-1,rg) && good(raw,az,rg) )
				scan(az,rg) = raw(az-1,rg)-raw(az,rg);
			else scan(az,rg) = 0.;
		}
	scan*=(0.5*raw.beam_count/M_PI);
	return scan;
}

template<typename T>
PolarScan<T> make_gradient_elevation_scan(const PolarScan<T>& low, const PolarScan<T>& up)
{
	PolarScan<T> scan(up);
	for(unsigned az=0;az<up.beam_count;az++)		
		for(unsigned rg=0;rg<up.beam_size;rg++)
		{
			if(good(up,az,rg) && good(low,az,rg))
				scan(az,rg) = up(az,rg)-low(az,rg);
			else scan(az,rg) = 0.;
		}
	scan*=abs(up.elevation-low.elevation)*M_PI/180.;
	return scan;
}

//===========================Volume manipulators=============================//

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
		window_size=1+2*std::floor(0.5*slope_range/raw[i].cell_size);
		vol.push_back(make_slope_scan(raw[i],window_size));
	}
}

template<typename T>
void textureSD(const Scans<T>& raw, Scans<T>& vol, double filter_range, bool force_check_undetect=false)
{
	textureSD(raw,vol,filter_range,0.,force_check_undetect);
}

template<typename T>
void textureSD(const Scans<T>& raw, Scans<T>& vol, double filter_range, double filter_azimuth=0. ,bool force_check_undetect=false)
{
	unsigned window_length;
	unsigned window_width;
	vol.clear();
	check_undetect=force_check_undetect;
	vol.quantity=raw.quantity;
	vol.units=raw.units;
	for(unsigned i=0;i<raw.size();i++)
	{
		window_length=1+2*std::floor(0.5*filter_range/raw[i].cell_size);
		window_width=1+2*std::floor(0.5*filter_azimuth/(360./raw[i].beam_count));
		vol.push_back(make_rms_scan(raw[i], window_length, window_width));
	}
}

template<typename T>
void filter(const Volume<T>& raw, Volume<T>& vol, double filter_range, bool force_check_undetect=false)
{
	filter(raw,vol,filter_range,0.,force_check_undetect);
}

template<typename T>
void filter(const Volume<T>& raw, Volume<T>& vol, double filter_range, double filter_azimuth=0., bool force_check_undetect=false)
{
	unsigned window_length;
	unsigned window_width;
	vol.clear();
	check_undetect=force_check_undetect;
	vol.quantity=raw.quantity;
	vol.units=raw.units;
	for(unsigned i=0;i<raw.size();i++)
	{
		window_length=1+2*std::floor(0.5*filter_range/raw[i].cell_size);
		window_width=1+2*std::floor(0.5*filter_azimuth/(360./raw[i].beam_count));
		vol.push_back(make_filter_scan(raw[i], window_length, window_width));
	}
}

template<typename T>
void gradient_azimuth(const Volume<T>& raw, Volume<T>& vol, bool force_check_undetect=false)
{
	vol.clear();
	check_undetect=force_check_undetect;
	vol.quantity=raw.quantity;
	vol.units=raw.units;
	for(unsigned el=0;el<raw.size();el++)
		vol.push_back(make_gradient_azimuth_scan(raw[el]));
}

template<typename T>
void gradient_elevation(const Volume<T>& raw, Volume<T>& vol, bool force_check_undetect=false)
{
	vol.clear();
	check_undetect=force_check_undetect;
	vol.quantity=raw.quantity;
	vol.units=raw.units;
	vol.push_back(make_gradient_elevation_scan(raw[1]),raw[0]);
	for(unsigned el=1;el<raw.size();el++)
		vol.push_back(make_gradient_elevation_scan(raw[el-1],raw[el]));
}

template<typename T>
void lin2dB(Volume<T>& lin, Volume<T>& dB)
{
	//this->quantity=lin.quantity.lin2dB(); // TODO: not yet implemented
	dB.clear();
	for(unsigned i=0;i<lin.size();i++)
	{
		dB.push_back(PolarScan<T>(lin[i].beam_count,lin[i].beam_size,0.));
		dB[i].cell_size = lin[i].cell_size;
		dB[i].elevation = lin[i].elevation;
		dB[i].block(0,0,lin[i].beam_count,lin[i].beam_size)=lin[i].log10();		
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
		lin.push_back(PolarScan<T>(dB[i].beam_count,dB[i].beam_size,0.));
		lin[i].cell_size = dB[i].cell_size;
		lin[i].elevation = dB[i].elevation;
		lin[i].block(0,0,dB[i].beam_count,dB[i].beam_size) = dB[i]*0.1;
		lin[i].block(0,0,dB[i].beam_count,dB[i].beam_size) = lin[i].exp10();
	}
}

template<typename T>
void lin2dB(Volume<T>& lin)
{
	//this->quantity=lin.quantity.lin2dB(); // TODO: not yet implemented
	for(unsigned i=0;i<lin.size();i++)
		lin[i].block(0,0,lin[i].beam_count,lin[i].beam_size)=lin[i].log10();
	lin*=10.;
}

template<typename T>
void dB2lin(Volume<T>& dB)
{
	//this->quantity=DB.quantity.dB2lin(); // TODO: not yet implemented
	dB*=0.1;
	for(unsigned i=0;i<dB.size();i++)
		dB[i].block(0,0,dB[i].beam_count,dB[i].beam_size) = dB.scan(i).exp10();
}

}	// namespace volume
}	// namespace elaboradar
