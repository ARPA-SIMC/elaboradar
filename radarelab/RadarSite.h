/**
 *  @file
 *  @ingroup radarelab
 *  @brief Radar Site description
 */
#ifndef RADARELAB_RADAR_SITE_
#define RADARELAB_RADAR_SITE_

#include <limits>
#include <iostream>
#include <string>

using namespace std;

 /*
 * \brief Class to store radar site coordinate
 */
class RadarSite{
  public:
	/*! lat_r */
 	/*! Latitude of the radar place */
//	/*! \sa checkLatitude*/
	double lat_r;
	/*! lon_r */
	/*! Longitude of radar place */
//	/*! \sa checkLongitude*/
	double lon_r;
	/*! height_r */
	/*! Height above the msl of radar position (antenna tower not included)*/
	double height_r = 0. ;
	/*! antennaTowerHeight */
	/*! Height of the electric antenna focus related above the surface*/
	double antennaTowerHeight = 0.;
	/*! Odim source attributea */
        std::string source ;

  public:
	/*! 
         *\brief Constructor
	 */
	/*!
	 * Inizialize radar site coordinates at \n
	 * Lat                 : 0. N - 
	 * Lon                 : 0. E - 
	 * Height              : 0. m - 
	 * Antennatower height : 0. m \n
	 * source 	       : ""\n 
	 */ 
	RadarSite() {
	  setRadarCoord(0.,0.,0.,0.);
	  source = "";
	}
	/*!
         *\brief Constructor with radar coordinates passed
	 *\param latr    - Radar site Latitude
	 *\param lonr    - Radar site Longitude
	 *\param heightr - Height above the msl of radar site (antenna tower not included)
	 *\param aTH     - Height of the electric antenna focus related above the surface
	 *\param source  - Source site description (Odim format) 
	 */
	RadarSite(double latr,double lonr,double heightr, double aTH, std::string source) 
	   : lat_r(latr), lon_r(lonr),height_r(heightr),antennaTowerHeight(aTH), source(source)
	{
	}

/** 
 *  Constructor
 *  Copy from another RadarSite
 *  @param [in] v - RadarSite object to be copied
 */  
    RadarSite(const RadarSite& v)
    {
        this->lat_r              = v.lat_r              ;
        this->lon_r              = v.lon_r              ;
        this->height_r           = v.height_r           ;
        this->antennaTowerHeight = v.antennaTowerHeight ;
	this->source             = v.source             ;
    }

	void setRadarCoord(float latr,float lonr,float heightr, float aTH) {
        this->lat_r              = latr              ;
        this->lon_r              = lonr              ;
        this->height_r           = heightr           ;
        this->antennaTowerHeight = aTH ;
	
	}
	double getTotalHeight () const { return (this->height_r + this->antennaTowerHeight) ; }
private:
};

#endif


