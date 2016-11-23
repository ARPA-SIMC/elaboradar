#include "cartproducts.h"
#include "assets.h"
#include <radarlib/radar.hpp>
#include <proj_api.h>

using namespace radarelab;
using namespace OdimH5v21;
namespace elaboradar {

CartProducts::CartProducts(const Volume<double>& volume, unsigned image_side, unsigned sample_square_size)
    : mapping(volume[0].beam_size),
      fullres(volume[0].beam_size),
      scaled(mapping, image_side, sample_square_size),
      z_out(image_side), z_fr(image_side*sample_square_size), top_1x1(image_side), qual_Z_1x1(image_side),
      quota_1x1(image_side), dato_corr_1x1(image_side),
      elev_fin_1x1(image_side), beam_blocking_1x1(image_side),
      neve_1x1(image_side), corr_1x1(image_side), conv_1x1(image_side)
{
    logging_category = log4c_category_get("radar.cart");

    quota_1x1.fill(128);

    LOG_INFO("Creazione Matrice Cartesiana");
    fullres.map_max_sample(volume[0], mapping);
    //assets.write_gdal_image(fullres.map_azimuth, "DIR_DEBUG", "map_azimuth", "PNG");
    //assets.write_gdal_image(fullres.map_range, "DIR_DEBUG", "map_range", "PNG");

    LOG_INFO("Creazione Matrice Cartesiana ridimensionata");
    scaled.map_max_sample(volume[0], fullres);

    FullsizeRes = volume.at(0).cell_size;
    ScaledRes   = FullsizeRes*sample_square_size;

}

void CartProducts::write_out(Assets& assets)
{
    assets.write_image(z_out, "OUTPUT_Z_LOWRIS_DIR", ".ZLR", "file output 1X1");
    assets.write_image(top_1x1, "DIR_QUALITY", ".top20_ZLR", "file top20");
    assets.write_image(qual_Z_1x1, "OUTPUT_Z_LOWRIS_DIR", ".qual_ZLR", "file qualita' Z");
    assets.write_image(quota_1x1, "DIR_QUALITY", ".quota_ZLR", "file qel1uota");
    assets.write_image(dato_corr_1x1, "DIR_QUALITY", ".anap_ZLR", "file anap");
    assets.write_image(elev_fin_1x1, "DIR_QUALITY", ".elev_ZLR", "file elev");
    assets.write_image(beam_blocking_1x1, "DIR_QUALITY", ".bloc_ZLR", "file bloc");
    assets.write_image(neve_1x1, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_image(corr_1x1, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_image(conv_1x1, "DIR_QUALITY", ".conv_ZLR", "punti convettivi");
}

void CartProducts::write_out(Assets& assets, unsigned image_side)
{
    assets.write_subimage(z_out, image_side, "OUTPUT_Z_LOWRIS_DIR", ".ZLR", "file output 1X1");
    assets.write_subimage(top_1x1, image_side, "DIR_QUALITY", ".top20_ZLR", "file top20");
    assets.write_subimage(qual_Z_1x1, image_side, "OUTPUT_Z_LOWRIS_DIR", ".qual_ZLR", "file qualita' Z");
    assets.write_subimage(quota_1x1, image_side, "DIR_QUALITY", ".quota_ZLR", "file qel1uota");
    assets.write_subimage(dato_corr_1x1, image_side, "DIR_QUALITY", ".anap_ZLR", "file anap");
    assets.write_subimage(elev_fin_1x1, image_side, "DIR_QUALITY", ".elev_ZLR", "file elev");
    assets.write_subimage(beam_blocking_1x1, image_side, "DIR_QUALITY", ".bloc_ZLR", "file bloc");
    assets.write_subimage(neve_1x1, image_side, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_subimage(corr_1x1, image_side, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_subimage(conv_1x1, image_side, "DIR_QUALITY", ".conv_ZLR", "punti convettivi");
}

void CartProducts::write_out(Assets& assets, unsigned image_side,std::string algos)
{
    assets.write_subimage(z_out, image_side, algos, "OUTPUT_Z_LOWRIS_DIR", ".ZLR", "file output 1X1");
    assets.write_subimage(top_1x1, image_side, algos, "DIR_QUALITY", ".top20_ZLR", "file top20");
    assets.write_subimage(qual_Z_1x1, image_side, algos, "OUTPUT_Z_LOWRIS_DIR", ".qual_ZLR", "file qualita' Z");
    assets.write_subimage(quota_1x1, image_side, algos, "DIR_QUALITY", ".quota_ZLR", "file qel1uota");
    assets.write_subimage(dato_corr_1x1, image_side, algos, "DIR_QUALITY", ".anap_ZLR", "file anap");
    assets.write_subimage(elev_fin_1x1, image_side, algos, "DIR_QUALITY", ".elev_ZLR", "file elev");
    assets.write_subimage(beam_blocking_1x1, image_side, algos, "DIR_QUALITY", ".bloc_ZLR", "file bloc");
    assets.write_subimage(neve_1x1, image_side, algos, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_subimage(corr_1x1, image_side, algos, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_subimage(conv_1x1, image_side, algos, "DIR_QUALITY", ".conv_ZLR", "punti convettivi");
}

   

    void CartProducts::write_odim(Assets& assets, unsigned image_side, std::string algos, OdimProdDefs &odimProd)
{
    const char* dir = getenv("OUTPUT_Z_LOWRIS_DIR");
    if (!dir)
    {
        LOG_INFO("$%s not set", "OUTPUT_Z_LOWRIS_DIR");
        throw runtime_error("required env var is not set");
    }
    SourceInfo OdimSource(assets.getRadarSite().source);
    
    string fname = string(dir) + "/" + OdimSource.OperaRadarNode + "_" + assets.fname_from_acq_time() + "_" + std::to_string(image_side) + "_"+algos+".h5";

    OdimH5v21::OdimFactory*	factory = NULL;
    OdimH5v21::ImageObject*	image	= NULL;
 
    factory = new OdimFactory();				/* creazione di una factory */	
    image   = factory->createImageObject(fname);		/* creazione di un oggetto generico */

    /* what (valori obbligatori e comuni a tutti gli oggetti) */
    image->setDateTime(assets.getAcqTime());
    image->setSource(OdimSource);

		/* where */
    std::string proj ="+proj=aeqd +lat_0=";
    proj = proj+std::to_string(assets.getRadarSite().lat_r)+"N +lon_0="+std::to_string(assets.getRadarSite().lon_r)+"E  +units=m +datum=WGS84";
    projPJ pj_aeqd, pj_latlong;
    std::string LatLon_def ("+proj=latlong +datum=WGS84");
    if (!(pj_aeqd = pj_init_plus(proj.c_str())) ) 
       exit(1);
    if (!(pj_latlong = pj_init_plus(LatLon_def.c_str())) )
       exit(1);
    double coord_min =  -(image_side * odimProd.prodRes *0.5) ;
    double coord_max =    image_side * odimProd.prodRes *0.5   ;
    double x[]={coord_min, coord_max, coord_min, coord_max};		// { LL , LR, UL, UR }
    double y[]={coord_min, coord_min, coord_max, coord_max};		// { LL , LR, UL, UR }
    if (pj_transform(pj_aeqd, pj_latlong,  4, 1, x, y, NULL ) != 0 ) exit(1000);
    image->setLL_Latitude (y[0]*RAD_TO_DEG);
    image->setLL_Longitude(x[0]*RAD_TO_DEG);
    image->setLR_Latitude (y[1]*RAD_TO_DEG);
    image->setLR_Longitude(x[1]*RAD_TO_DEG);
    image->setUL_Latitude (y[2]*RAD_TO_DEG);
    image->setUL_Longitude(x[2]*RAD_TO_DEG);
    image->setUR_Latitude (y[3]*RAD_TO_DEG);
    image->setUR_Longitude(x[3]*RAD_TO_DEG);

    image->setXSize(image_side);
    image->setYSize(image_side);
    image->setXScale(odimProd.prodRes);
    image->setYScale(odimProd.prodRes);
    image->setProjectionArguments(proj);

		/* how */
    image->setTaskOrProdGen(algos);
    image->setStartEpochs((unsigned int)assets.getAcqTime());
    image->setEndEpochs(  (unsigned int)assets.getAcqTime());
    image->setSystem(odimProd.System);
    image->setSoftware(odimProd.System);
    image->setSoftwareVer(odimProd.System);

/*   ------------------------
 *   QUI BISOGNA METTERE QUALCOSA CHE CREI IL PRODOTTO GIUSTO MA NON SO COME SI FA
 *
 *   AL MOEMNTO CREO SOLO QUELLO CHE MI SERVE
    switch (odimProd.ProductType) {
	extern RADAR_API const char* PRODUCT_SCAN;
	extern RADAR_API const char* PRODUCT_PPI; 
	extern RADAR_API const char* PRODUCT_CAPPI; 
	extern RADAR_API const char* PRODUCT_PCAPPI; 
	extern RADAR_API const char* PRODUCT_ETOP; 
	extern RADAR_API const char* PRODUCT_MAX; 
	extern RADAR_API const char* PRODUCT_RR; 
	extern RADAR_API const char* PRODUCT_VIL; 
	extern RADAR_API const char* PRODUCT_COMP; 
	extern RADAR_API const char* PRODUCT_VP; 
	extern RADAR_API const char* PRODUCT_RHI; 
	extern RADAR_API const char* PRODUCT_XSEC; 
	extern RADAR_API const char* PRODUCT_VSP; 
	extern RADAR_API const char* PRODUCT_HSP;
	extern RADAR_API const char* PRODUCT_RAY; 
	extern RADAR_API const char* PRODUCT_AZIM; 
	extern RADAR_API const char* PRODUCT_QUAL;
	extern RADAR_API const char* PRODUCT_LBM_ARPA;
	extern RADAR_API const char* PRODUCT_POH_ARPA;

 *
 */
    Product_LBM * dataset =  image->createProductLBM();
    dataset->setProduct("SURF");// Dato che ODIMh5 v2.2 codifica per i dati riferiti alla superfice il prodotto generico SURF sovrascrivo quindi il tipo di prodotto definito LBM 
    dataset->setStartDateTime(assets.getAcqTime());
    dataset->setEndDateTime  (assets.getAcqTime());
    Product_2D_Data* data = dataset->createQuantityData(odimProd.Quantity);
    data->setNodata(odimProd.Nodata);
    data->setUndetect(odimProd.Undetect);
    data->setOffset(odimProd.Offset);
    data->setGain(odimProd.QuantityDynamics/255.);

    unsigned xofs = (odimProd.prodField.cols() - image_side) / 2;
    unsigned yofs = (odimProd.prodField.rows() - image_side) / 2;
    OdimH5v21::DataMatrix <unsigned char> field (image_side,image_side,255);
    for (unsigned y = 0; y < image_side; ++y)
        for (unsigned x = 0; x < image_side; ++x)
            if (odimProd.prodField(y + yofs, x + xofs) == 0 ) field.elem(y,x) = 255 ;
            else if (odimProd.prodField(y + yofs, x + xofs) == 255) field.elem(y,x) = 254 ;
            else if (odimProd.prodField(y + yofs, x + xofs) == 1) field.elem(y,x) = 0 ;
            else field.elem(y,x) = odimProd.prodField(y + yofs, x + xofs);
    data->writeData(field);
    if (odimProd.SaveQuality) {
      OdimQuality * quality =data->createQualityData();
      quality->getWhat()->set(ATTRIBUTE_WHAT_OFFSET,		0.);
      quality->getWhat()->set(ATTRIBUTE_WHAT_GAIN,		0.01);
      quality->getHow() ->set(ATTRIBUTE_HOW_TASK,	"Anna Fornasiero");
      OdimH5v21::DataMatrix <unsigned char> Qfield (image_side,image_side,255);
      for (unsigned y = 0; y < image_side; ++y)
        for (unsigned x = 0; x < image_side; ++x)
            Qfield.elem(y,x) = odimProd.QualityField(y + yofs, x + xofs);
      quality->writeQuality(Qfield);
      delete quality;
    }
    delete data;
    delete dataset;
    delete image;
    delete factory;

}

OdimProdDefs::OdimProdDefs(radarelab::Image<unsigned char> & prodField, double prodRes) 
   :  prodField (prodField), System("ARPA-SIMC"), ProductType("SURF"),
    Quantity(PRODUCT_QUANTITY_DBZH), Nodata(255), Undetect (0), Offset(-20.),
    QuantityDynamics (80), SaveQuality(false), prodRes(prodRes)
{
	// Nothing to do
}

OdimProdDefs::OdimProdDefs(radarelab::Image<unsigned char> & prodField, radarelab::Image<unsigned char> & QualityField, double prodRes) 
   :  prodField (prodField), System("ARPA-SIMC"), ProductType("SURF"),
    Quantity(PRODUCT_QUANTITY_DBZH), Nodata(255), Undetect (0), Offset(-20.),
    QuantityDynamics (80), SaveQuality(true), QualityField(QualityField), prodRes(prodRes)

{
	// Nothing to do
}

    

}	// namespace elaboradar

