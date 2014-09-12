/*!
 * \file
 * \brief Class to compute hydrometeor classification
 */

#include"volume.h"
#include"volume/odim.h"
#include"volume/loader.h"
#include"volume/resample.h"
#include"site.h"

//#include "algo/polar_cart.h"
//#include"melting_layer.h"
//#include"FIR_filter.h"

#include<string>
#include<iostream>
#include<fstream>

namespace elaboradar {
namespace volume {

/*
 * ==================================================================
 *        Enum: EchoClass
 * Description: classes of radar echoes
 * ==================================================================
 */
/*!
 * 	\brief List classes of radar echoes
 * 	Classes defined in Park et al. (2009)
 * 	Overload << operator to print class names
 *
 * 	TODO: manage how to print odim volume of classes
 * 	TODO: This declarations deserve their own file
 */

enum EchoClass {
	GC_AP,	// ground clutter or anomalous propagation
	BS,	// biological scatterers
	DS,	// dry aggregated snow
	WS,	// wet snow
	CR,	// crystals of various orientation
	GR,	// graupel
	BD,	// big drops
	RA,	// light and moderate rain
	HR,	// heavy rain
	RH	// mixture of rain and hail
};

inline std::ostream& operator<<(std::ostream& oss, EchoClass& ECl)
{
	switch(ECl)
	{
		case GC_AP:
			oss<<"GC_AP";
			break;
		case BS:
			oss<<"BS";
			break;
		case DS:
			oss<<"DS";
			break;
		case WS:
			oss<<"WS";
			break;
		case CR:
			oss<<"CR";
			break;
		case GR:
			oss<<"GR";
			break;
		case BD:
			oss<<"BD";
			break;
		case RA:
			oss<<"RA";
			break;
		case HR:
			oss<<"HR";
			break;
		case RH:
			oss<<"RH";
			break;
		default:
			oss<<"unknown echo type "<<ECl;
	}
	return oss;
}

/*
 * ==================================================================
 * 	 Class: PROB
 * Description: Matrix of probability (trapezoidal)
 * ==================================================================
 */
/*!
 * 	\brief Given radar variables compute matrix of probability
 */
class PROB : public Matrix2D<double>
{
public:
/*!
 * Constructor
 */
	PROB(double z, double zdr, double rhohv, double lkdp, double sdz, double sdphidp);
	
private:
/*!
 * f1 function
 */
	double f_1(double Z) {return -0.05+2.5e-3*Z+7.5e-4*Z*Z;}
/*!
 * f2 function
 */
	double f_2(double Z) {return 0.68-4.81e-2*Z+2.92e-3*Z*Z;}
/*!
 * f3 function
 */
	double f_3(double Z) {return 1.42+6.67e-2*Z+4.85e-4*Z*Z;}
/*!
 * g1 function
 */
	double g_1(double Z) {return -44.0+0.8*Z;}
/*!
 * g2 function
 */
	double g_2(double Z) {return -22.0+0.5*Z;}
/*!
 * trapezoidal probability function
 */
	double trap(double x1, double x2, double x3, double x4, double val);
/*!
 * vector of trapezoidal probability for a specific class of echo
 */
	Matrix2D<double> prob_class(EchoClass classe,double z, double zdr, double rhohv, double lkdp, double sdz, double sdphidp);
};

/*
 * ==================================================================
 * 	 Class: CONF
 * Description: Confidence vector
 * ==================================================================
 */
/*!
 * 	\brief compute confidence vector of radar variables
 */
class CONF : public Matrix2D<double>
{
public:
	CONF() 
	{
		this->resize(6,1);
		*this<<1.,1.,1.,1.,1.,1.;
	}
};

/*
 * ==================================================================
 * 	 Class: HCA_Park
 * Description: compute aggregation values A_i and return class of echo
 * ==================================================================
 */
/*!
 * 	\brief compute HCA
 */
class HCA_Park
{
public:
	/*================== MEMBERS ====================*/
/*!
 * Passed input variables
 */
	double z,zdr,rhohv,lkdp,sdz,sdphidp;
/*!
 * Vector of aggregation classes
 */
	Eigen::VectorXd Ai;

	/*================== METHODS ====================*/
/*!
 * Default constructor do nothing
 */
	HCA_Park(){}
/*!
 * Constructor
 */
	HCA_Park(double Z, double ZDR, double RHOHV, double LKDP, double SDZ, double SDPHIDP);
/*!
 * Search for non meteorological echoes
 */
	inline bool non_meteo_echo()
	{
		unsigned idx;
		Ai.maxCoeff(&idx);
		if(idx==0||idx==1) return true;
		else return false;
	}
/*!
 * Search for meteo echoes
 */
	inline bool meteo_echo()
	{
		unsigned idx;
		Ai.maxCoeff(&idx);
		if(idx==0||idx==1) return false;
		else return true;
	}

	EchoClass echo()
	{
		unsigned idx;
		Ai.maxCoeff(&idx);
		EchoClass hca=(EchoClass)idx;
		return hca;
	}
};


/*
 * ==================================================================
 * 	 Class: Classifier
 * Description: compute hydrometeor classification Park et al. (2009)
 * ==================================================================
 */
/*!
 * 	\brief compute hydrometeor classification
 * 	
 * 	This class compute volumes of hydrometeor classification using
 * 	the HCA from Park et al. (2009).
 * 	Input data are Z, Zdr, rhohv and phidp and are provided through
 * 	odim volume data file.
 */

class classifier
{
public:
	/*================== MEMBERS ====================*/
/*!
 * Polar volume of class of hydrometeors
 */
	Volume<EchoClass> vol_hca;
/*!
 * Volume of aggregation vectors (I can't make volumes of vectors, only scalars are allowed)
 */
	std::vector< std::vector< std::vector<HCA_Park> > > vol_Ai;
/*!
 * Pathname of the loaded file (if provided)
 */
	std::string pathname;
/*!
 * Polar Volume of copolar reflectivity DBZH
 */
	Volume<double> vol_z;
/*!
 * Polar Volume of differential reflectivity ZDR
 */
	Volume<double> vol_zdr;
/*!
 * Polar Volume of correlation coefficient rhohv
 */
	Volume<double> vol_rhohv;
/*!
 * Polar Volume of differential phase shift PHIDP
 */
	Volume<double> vol_phidp;
/*!
 * Polar Volume of radial velocity
 */
	Volume<double> vol_vrad;
/*!
 * Polar Volume of signal to noise ratio
 */
	Volume<double> vol_snr;
/*!
 * Filtered copolar reflectivity over 1 km window
 */
	Volume<double> vol_z_1km;
/*!
 * Filtered Zdr over 2km window
 */
	Volume<double> vol_zdr_2km;
/*!
 * Filtered rhohv over 2km window
 */
	Volume<double> vol_rhohv_2km;
/*!
 * Filtered differential phase over a 2km window
 */
	Volume<double> vol_phidp_2km;
/*!
 * Filtered differential phase over a 6km window
 */
	Volume<double> vol_phidp_6km;
/*!
 * Specific differential phase averaged over a 2km window
 */
	Volume<double> vol_lkdp_2km;
/*!
 * Specific differential phase averaged over a 6km window
 */
	Volume<double> vol_lkdp_6km;
/*!
 * Texture parameter SD of reflectivity volume
 */
	Volume<double> vol_sdz;
/*!
 * Texture parameter SD of differential phase volume
 */
	Volume<double> vol_sdphidp;


	/*================== METHODS ====================*/
/*!
 * \brief Constructor from odim file
 */
/*!
 * Initialize basic input variables (vol_z vol_zdr vol_rhohv and vol_phidp)
 * \param [in] file - pathname of the odim volume file to be inspected
 * \param [in] site - site object
 */ 
	classifier(const std::string& file, const Site& site);
/*!
 * \brief Initialize derived input data
 * could be part of the constructor
 */
	void compute_derived_volumes();
/*!
 * \brief correct phidp for sparse nodata and undetected
 * moving average filter along beam path
 */
	void correct_phidp();
/*!
 * \brief correct rhohv and zdr for noise (Schuur et al. 2003)
 */
	void correct_for_snr();
/*!
 * \brief correct Z and Zdr for path attenuation
 */
	void correct_for_attenuation();
/*!
 * \brief Initialize vol_lkdp
 * 10log10 of the moving average slope of phidp along beam path
 */
	void compute_lkdp();
/*!
 * \brief Compute Echo Classes
 * Park et al. (2009) HCA algorithm
 */
	void HCA_Park_2009();
};

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
	double Hmin,Hmax;	// km
	unsigned count;

	MLpoints(double minHeight,double maxHeight,unsigned az_count,unsigned height_count) 
	     : 	Matrix2D<unsigned>(Matrix2D::Constant(height_count,az_count,0)),
		Hmin(minHeight),Hmax(maxHeight),count(0) {}

	double azimuth_deg(unsigned az_idx){return (double)az_idx*360./(double)this->cols();}
	double azimuth_rad(unsigned az_idx){return (double)az_idx*2.*M_PI/(double)this->cols();}
	unsigned rad2idx(double rad){return (unsigned)(rad*(double)this->cols()/(2.*M_PI));}
	unsigned deg2idx(double deg){return (unsigned)(deg*(double)this->cols()/360.);}

	double height(unsigned h_idx){return Hmin+(double)h_idx*(Hmax-Hmin)/(double)this->rows();}
	unsigned h_idx(double height){return (unsigned)((height-Hmin)*(double)this->rows()/(Hmax-Hmin));}

	void box_top_bottom(double box_width_deg, double bot_th, double top_th, std::vector<double>& ML_b, std::vector<double>& ML_t);
	void save2file();
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
	std::vector<double> top;
	std::vector<double> bot;

	Volume<double> vol_z_0_5km;
	Volume<double> vol_zdr_1km;
	Volume<double> vol_rhohv_1km;

	MeltingLayer(Volume<double>& vol_z,Volume<double>& vol_zdr,Volume<double>& vol_rhohv, 
			std::vector< std::vector< std::vector< HCA_Park> > >& HCA);

	void fill_empty_azimuths();
	
};




} // namespace volume
} // namespace elaboradar
