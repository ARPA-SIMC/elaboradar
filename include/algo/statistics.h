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
	unsigned N;

	Statistic():sum_x(0),sum_y(0),sum_xy(0),sum_x2(0),N(0){}

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
	}

	void slim(T x)
	{
		sum_x-=x;
		sum_x2-=x*x;
		N--;
	}

	void clear()
	{
		sum_x=0;
		sum_y=0;
		sum_xy=0;
		sum_x2=0;
		N=0;
   	}

	T slope;
	T intercept;
	T variance;
	T dev_std;
	T mean;

	T compute_slope(unsigned minimum=2)
	{
		if(N>=minimum)
			slope = (N*sum_xy-sum_x*sum_y)/(N*sum_x*sum_x-sum_x2);
		else slope = slope/(slope-slope); // orribile modo di far ritornare NaN
		return slope;
	}

	T compute_intercept(unsigned minimum=2)
	{
		if(N>=minimum)
			intercept = (sum_y*sum_x2-sum_x*sum_xy)/(N*sum_x2-sum_x*sum_x);
		else return intercept = slope/(slope-slope);
		return intercept;
	}

	T compute_variance(unsigned minimum=2)
	{
		if(N>=minimum)
			variance = (sum_x2-sum_x*sum_x/(double)N)/((double)N-1.);
		else variance = variance/(variance-variance);
		return variance;
	}
	
	T compute_dev_std(unsigned minimum=2)
	{
		//compute_variance(minimum);
		//dev_std=sqrt(variance);
		if(N>=minimum)
			dev_std = sqrt(N*sum_x2-sum_x*sum_x)/N;
		else dev_std = dev_std/(dev_std-dev_std);
		return dev_std;
	}
	
	T compute_mean(unsigned minimum=2)
	{
		if(N>=minimum)
			mean = sum_x/(T)N;
		else mean = mean/(mean-mean);
		return mean;
	}
   
	T get_slope() {return slope;}
	T get_intercept() {return intercept;}
};

/*
template<typename T>
class LinearFit : public Statistic<T>
{
public:
	
	T slope;
	T intercept;

	LinearFit(){}

	T compute_slope()
	{
		if(N)
		{
			slope = (N*sum_xy-sum_x*sum_y)/(N*sum_x2-sum_x*sum_x);
			return slope;
		}
		else return slope/(slope-slope); // orribile modo di far ritornare NaN
	}

	T compute_intercept()
	{
		if(N)
		{
			intercept = (sum_y*sum_x2-sum_x*sum_xy)/(N*sum_x2-sum_x*sum_x);
			return intercept;
		}
		else return slope/(slope-slope);
	}
   
	T get_slope() {return slope;}
	T get_intercept() {return intercept;}
};
*/
}	// namespace elaboradar
