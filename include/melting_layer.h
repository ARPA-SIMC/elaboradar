/*!
 * \file
 * \brief Class to compute melting layer limits
 */

#include"volume.h"
//#include"volume/odim.h"
//#include"volume/loader.h"

#include<string>
#include<iostream>

namespace cumbac{
namespace volume{

/*
 * ==================================================================
 *	      Class: MLpoints
 *	Description: store ML points in a azimuth height matrix
 * ==================================================================
 */
/*!
 * 	\brief MLpoints
 * 	Melting Layer Points matrix AzH
 */
class MLpoints : public Matrix2D<unsigned>
{
public:
	double Hmin,Hmax;

	MLpoints(double minHeight,double maxHeight,unsigned az_count,unsigned height_count) 
	     : 	Matrix2D<unsigned>(Matrix2D::Constant(az_count,height_count,0)),
		Hmin(minHeight),Hmax(maxHeight) {}
};


/*
 * ==================================================================
 *	      Class: MeltingLayer
 *	Description: compute melting layer boundaries from polarimetry
 * ==================================================================
 */
/*!
 * 	\brief MeltingLayer
 * 	Melting Layer Detection Algorithm MLDA
 * 	Giangrande et al. 2008
 */
class MeltingLayer
{
public:
	double ML_top;
	double ML_bot;

	Volume<double> vol_z_0_5km;
	Volume<double> vol_zdr_1km;
	Volume<double> vol_rhohv_1km;

	MeltingLayer(Volume<double>& vol_z,Volume<double>& vol_zdr,Volume<double>& vol_rhohv);

};

} // namespace volume
} // namespace cumbac
