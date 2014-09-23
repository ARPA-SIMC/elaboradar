#include<cmath>

namespace elaboradar {

template<typename T>
class Statistic
{
public:
	T sum_x;
	T sum_y;
	T sum_xy;
	T sum_x2;

	T slope;
	T intercept;
	T variance;
	T dev_std;
	T mean;
	T M2;
	T delta;
	unsigned N;

	Statistic():sum_x(0),sum_y(0),sum_xy(0),sum_x2(0),N(0),mean(0),M2(0){}

	void feed(T x, T y)
	{
		sum_x+=x;
		sum_y+=y;
		sum_xy+=x*y;
		sum_x2+=x*x;
		N++;
	}
	
	void slim(T x, T y)
	{
		sum_x-=x;
		sum_y-=y;
		sum_xy-=x*y;
		sum_x2-=x*x;
		N--;
	}

	void feed(T x)
	{
		sum_x+=x;
		sum_x2+=x*x;
		N++;

		delta = x-mean;
		mean = mean+delta/(double)N;
		M2 = M2+delta*(x-mean);
	}

	void slim(T x)
	{
		sum_x-=x;
		sum_x2-=x*x;
		N--;

		if(N)
		{
			delta = x-mean;
			mean = mean-delta/(double)N;
			M2 = M2-delta*(x-mean);
		}
		else
		{
			mean=0.;
			M2=0;
		}
	}

	void clear()
	{
		sum_x=0;
		sum_y=0;
		sum_xy=0;
		sum_x2=0;
		N=0;
		mean=0.;
		M2=0.;
   	}

	T compute_slope(unsigned minimum=2)
	{
		if(N>=minimum)
			slope = (N*sum_xy-sum_x*sum_y)/(N*sum_x*sum_x-sum_x2);
		else slope = sqrt(-1); // NaN
		return slope;
	}

	T compute_intercept(unsigned minimum=2)
	{
		if(N>=minimum)
			intercept = (sum_y*sum_x2-sum_x*sum_xy)/(N*sum_x2-sum_x*sum_x);
		else intercept = sqrt(-1);
		return intercept;
	}

	T compute_variance(unsigned minimum=1)
	{
		if(M2<0.1e-11)M2=0.;
		if(N>=minimum) return M2/(double)N;
		else return sqrt(-1);
	}
	
	T compute_dev_std(unsigned minimum=1)
	{
		if(N>=minimum) return sqrt(compute_variance());
		else return sqrt(-1);
	}
	
	T compute_mean(unsigned minimum=1)
	{
		if(N>=minimum)
			return mean;
		else return sqrt(-1);
	}
};

}	// namespace elaboradar
