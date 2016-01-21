/*! \file statistics.h
 *  \ingroup radarelab 
 *  \brief Class to manage statistical information about streamed data
 */
#include <cmath>

namespace radarelab {
/*
 * ==================================================================
 *       Class: Statistic
 * Description: Generic class to perform basic statistics
 * ==================================================================
 */
/*! \class Statistic
 *  \brief Generic Class to perform statistical analysis
 *  Statistic object could be used as accumulator of data
 *  and performs statistical computations such as
 *  mean, variance, standard deviation, linear fit.
 */
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

/*! \fn Default constructor
 *  \brief Sets all member variables to zero
 */
	Statistic():sum_x(0),sum_y(0),sum_xy(0),sum_x2(0),mean(0),M2(0),N(0){}

/*! \fn feed
 *  \brief Feed accumulators with bidimensional numbers
 *  @param[in] x first coordinate of the bidimensional number
 *  @param[in] y second coordinate of the bidimensional number
 */
	void feed(T x, T y)
	{
		sum_x+=x;
		sum_y+=y;
		sum_xy+=x*y;
		sum_x2+=x*x;
		N++;
	}

/*! \fn slim
 *  \brief Take off a specified bidimensional number from the statistic
 *  @param[in] x first coordinate of the bidimensional number
 *  @param[in] y second coordinate of the bidimensional number
 */	
	void slim(T x, T y)
	{
		sum_x-=x;
		sum_y-=y;
		sum_xy-=x*y;
		sum_x2-=x*x;
		N--;
	}

/*! \fn feed
 *  \brief Feed accumulators with scalars
 *  @param[in] x scalar value to be accumulated
 */
	void feed(T x)
	{
		sum_x+=x;
		sum_x2+=x*x;
		N++;

		delta = x-mean;
		mean = mean+delta/(double)N;
		M2 = M2+delta*(x-mean);
	}

/*! \fn slim
 *  \brief Take off a specified scalar from the statistic
 *  @param[in] x scalar value to be deleted from the accumulator
 */
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

/*! \fn clear
 *  \brief Reset the Statistic object to a null state
 */
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

/*! \fn compute_slope
 *  \brief Compute least square linear fit (no exception if only scalars where provided)
 *  @param[in] minimum Minimum amount of accumulated data to perform the computation
 *  \return value of the slope. If \a N is less the \a minimum returns NaN
 */
	T compute_slope(unsigned minimum=2)
	{
		if(N>=minimum)
			slope = (N*sum_xy-sum_x*sum_y)/(N*sum_x*sum_x-sum_x2);
		else slope = sqrt(-1); // NaN
		return slope;
	}

/*! \fn compute_intercept
 *  \brief Compute least square linear fit (no exception if only scalars where provided)
 *  @param[in] minimum Minimum amount of accumulated data to perform the computation
 *  \return value of the intercept. If \a N is less the \a minimum returns NaN
 */
	T compute_intercept(unsigned minimum=2)
	{
		if(N>=minimum)
			intercept = (sum_y*sum_x2-sum_x*sum_xy)/(N*sum_x2-sum_x*sum_x);
		else intercept = sqrt(-1);
		return intercept;
	}

/*! \fn compute_variance
 *  \brief Compute variance of the distribution of \a x values
 *  @param[in] minimum Minimum amount of accumulated data to perform the computation
 *  \return value of the variance. If \a N is less the \a minimum returns NaN
 */
	T compute_variance(unsigned minimum=1)
	{
		if(M2<0.1e-11)M2=0.;
		if(N>=minimum) return M2/(double)N;
		else return sqrt(-1);
	}

/*! \fn compute_dev_std
 *  \brief Compute standard deviation of the distribution of \a x values
 *  @param[in] minimum Minimum amount of accumulated data to perform the computation
 *  \return value of the standard deviation. If \a N is less the \a minimum returns NaN
 */	
	T compute_dev_std(unsigned minimum=1)
	{
		if(N>=minimum) return sqrt(compute_variance());
		else return sqrt(-1);
	}

/*! \fn compute_mean
 *  \brief Compute mean of the distribution of \a x values
 *  @param[in] minimum Minimum amount of accumulated data to perform the computation
 *  \return value of the mean. If \a N is less the \a minimum returns NaN
 */
	T compute_mean(unsigned minimum=1)
	{
		if(N>=minimum)
			return mean;
		else return sqrt(-1);
	}
};

}	// namespace radarelab
