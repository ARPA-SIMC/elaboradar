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
