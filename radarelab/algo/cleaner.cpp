/*
 * =================================================================================
 *
 *       Filename:  volume_cleaner.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  18/02/2014 12:19:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =================================================================================
 */
#include "cleaner.h"
#include "elabora_volume.h"
#include "radarelab/image.h"
#include "radarelab/matrix.h"

namespace radarelab {
namespace algo {

using namespace std;
using namespace radarelab;

//--------------------------------------------------------------------------------
//  These methods use only VRAD and WRAD values to clean the beam
//
std::vector<bool> Cleaner::clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, int iray) const
{
    const unsigned beam_size = beam_z.rows();
    vector<bool> res(beam_size, false);
    bool in_a_segment = false;
    unsigned start, end;
    unsigned segment_length;
    bool before, after;
    unsigned counter = 0;

    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
//printf(" %4d %4d  %6.2f %6.2f %10.6f ",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin));
        if (!in_a_segment)
        {
            /* cerco la prima cella segmento da pulire*/
            if (beam_w(ibin) == W_threshold && beam_v(ibin) == bin_wind_magic_number)
            {
//printf(" 1  ----- START SEGMENT ------");
                in_a_segment = true;
                start = ibin;
                after = false;
                before = false;
            }
//	    else printf(" 0 ");
        } else {
            /* cerco la fine segmento da pulire*/
            if (beam_w(ibin) != W_threshold || beam_v(ibin) != bin_wind_magic_number || ibin == (beam_size - 1))
            {
                in_a_segment = false;
                end = ibin - 1;
                if (ibin == (beam_size - 1)) end = ibin;  // caso particolare per fine raggio
                /* Fine trovata ora procedo alla pulizia eventuale */
                segment_length = end - start;
                counter = counter + (unsigned)(segment_length);
		
		unsigned c_b=0;
		unsigned c_a=0;
                /* Cerco dati validi in Z prima del segmento */
                for (int ib = ibin - 12; ib < (signed)ibin; ++ib)
                    if (ib >= 0 && beam_z(ib) > Z_missing)
                  c_b++;                        
		if (c_b > 0.25*12) before = true;

                /* Cerco dati validi in Z dopo il segmento */
                for (unsigned ia = ibin + 1; ia <= ibin + 12; ++ia)
                    if (ia < beam_size && beam_z(ia) >= Z_missing)
			c_a++;                        
		if (c_a > 0.25*12) after = true;

//printf(" 0 ----- STOP SEGMENT ------ %4d  --  %4d    before %d %d  after %d %d ",segment_length,counter, c_b,before, c_a, after);
                if ((segment_length >= min_segment_length && !before && !after) ||
                        segment_length >= max_segment_length || counter > 100)
                {
                    /* qui pulisco */
                //            printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
                    for (unsigned ib = start; ib <= end; ++ib)
                        if( beam_z(ib) >  Z_missing)res[ib] = true;
                }
            }
//	    else printf(" 1 ");
        }
//printf("\n");
    }
    return res;
}

std::vector<unsigned char> Cleaner::eval_clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, int iray) const
{
    const unsigned beam_size = beam_z.rows();
    vector<unsigned char> res(beam_size, 0);
    bool in_a_segment = false;
    unsigned start = 0, end = 0;
    unsigned segment_length;
    bool before = false, after = false;
    unsigned counter = 0;

    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
//printf(" %4d %4d  %6.2f %6.2f %10.6f ",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin));
        if (!in_a_segment)
        {
            /* search the first radar bin's segment to be cleaned*/
            if (beam_w(ibin) == W_threshold && beam_v(ibin) == bin_wind_magic_number)
            {
//printf(" 1  ----- START SEGMENT ------");
                in_a_segment = true;
                start = ibin;
                after = false;
                before = false;
            }
        } else {
            /* search the last radar bin's segment to be cleaned*/
            if (beam_w(ibin) != W_threshold || beam_v(ibin) != bin_wind_magic_number || ibin == (beam_size - 1))
            {
                in_a_segment = false;
                end = ibin - 1;
                if (ibin == (beam_size - 1)) end = ibin;  // beam ended
                /* Fine trovata ora procedo alla pulizia eventuale */
                segment_length = end - start;
                counter = counter + (unsigned)(segment_length);
		
		unsigned c_b=0;
		unsigned c_a=0;
                /* Cerco dati validi in Z prima del segmento */
                for (int ib = ibin - 12; ib < (signed)ibin; ++ib)
                    if (ib >= 0 && beam_z(ib) > Z_missing)
                      c_b++;                        
		if (c_b > 0.25*12) before = true;

                /* Cerco dati validi in Z dopo il segmento */
                for (unsigned ia = ibin + 1; ia <= ibin + 12; ++ia)
                    if (ia < beam_size && beam_z(ia) >= Z_missing)
			c_a++;                        
		if (c_a > 0.25*12) after = true;

//printf(" 0 ----- STOP SEGMENT ------ %4d  --  %4d    before %d %d  after %d %d ",segment_length,counter, c_b,before, c_a, after);
                if ((segment_length >= min_segment_length && !before && !after) ||
                        segment_length >= max_segment_length || counter > 100)
                {
                    /* qui pulisco */
                //            printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
                    for (unsigned ib = start; ib <= end; ++ib)
                        if( beam_z(ib) >  Z_missing)res[ib] = 1;
                }
            }
//	    else printf(" 1 ");
        }
//printf("\n");
    }
    return res;
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//  This method use VRAD, WRAD, sdDBZH, sdZDR values to clean the beam
//
std::vector<bool> Cleaner::clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdzdr, PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& SD, int iray) const
{
    const unsigned beam_size = beam_z.rows();
    vector<bool> res(beam_size, false);
    bool in_a_segment = false;
    unsigned start = 0, end;
    unsigned segment_length;
    bool before, after;
    unsigned counter = 0;
    unsigned counter_trash = 0;
    unsigned counter_clutter =0;
    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
	bool is_clutter = false;
	bool is_trash = false;
	unsigned flag = 0 ;
//  In our systems (ARPA ER) interferences and other non meteo echo are characterised by the following steps
//
//  1)  Wind is not defined ad spectrumWidth is 0. with Z defined.
        if ( beam_w(ibin) == W_threshold && beam_v(ibin) == bin_wind_magic_number && beam_z (ibin) != Z_missing ) {
//  2)  Std<Dev of ZDR coulb be close to 0
           if( beam_sdzdr(ibin) <= 0.01 ){
		is_trash = true;
		flag=2;
           } else {
	     if (beam_z (ibin) >= 45. ){
//  2) inside thunderstorms (Z > 45) StdDev of Zdr and StdDev of Z are quite high) 
	       if ((ibin >100 && double(counter_trash)/double(ibin) >=0.5 &&  ( beam_sdzdr(ibin) >1 || beam_sd (ibin) > 5. )) ||
                   (beam_sdzdr(ibin) >4.0 && beam_sd (ibin) > 20.) ) {
	         is_trash = true;
	         flag=2;
               } else {
		 is_trash = false;
		 flag=0;
               }
           } else if ( (ibin >100 && double(counter_trash)/double(ibin) >=0.5 &&  ( beam_sdzdr(ibin) >1 || beam_sd (ibin) > 5. )) ||
                       (beam_sd (ibin) >2. && (beam_sdzdr(ibin) >2.0 || beam_sd (ibin) > 10. )) ) {
//  2) outside thunderstorms (Z > 45) StdDev of Zdr and StdDev of Z are lower  
	        is_trash = true;
	        flag=2;
             }
           }
         } else { 
//  3) Clutter is characterised by low value of VRAD and WRAD
	    if  ((beam_w(ibin) * fabs(beam_v(ibin)))  <= 0.25 && beam_z (ibin) != Z_missing ) {
	        is_clutter = true;
                flag = 1;	  
	    }
	 }
         if( is_clutter) counter_clutter ++;
         if( is_trash  ) counter_trash ++;

//printf(" %4d %4d  %6.2f %6.2f %10.6f %6.2f %6.2f ",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin), beam_sd(ibin),beam_sdzdr(ibin));
//printf("     -----    %2x %2x %2x %2x ",(unsigned char)((beam_z(ibin)-scan_z.offset)/scan_z.gain/256),
//(unsigned char)((beam_v(ibin)-scan_v.offset)/scan_v.gain/256),
//(unsigned char)((beam_w(ibin)-scan_w.offset)/scan_w.gain/256),
//(unsigned char)((beam_sd(ibin)-SD.offset)/SD.gain/256));
        if (!in_a_segment)
        {
            /* cerco la prima cella segmento da pulire*/
            if ( is_clutter || is_trash  )
            {
// printf(" %1d  ----- START SEGMENT ------",flag);
               in_a_segment = true;
                start = ibin;
                after = false;
                before = false;
            } 
//	    else printf(" %1d ",flag);
        } else {
            /* cerco la fine segmento da pulire*/
            if ( ! (is_clutter || is_trash ) || ibin == (beam_size - 1)) 
            {
                in_a_segment = false;
                end = ibin - 1;
                if (ibin == (beam_size - 1)) end = ibin;  // caso particolare per fine raggio
                /* Fine trovata ora procedo alla pulizia eventuale */
                segment_length = end - start+1;
                counter = counter + (unsigned)(segment_length);

/* 	il segmento è corto allora cerco nei dintorni dei dati validi, se li trovo non pulisco */
		if (segment_length <= 2*min_segment_length ){
              /* Cerco dati validi in Z prima del segmento */
	     	  int count=0;
                  for (int ib = ibin - 2*min_segment_length; ib < (signed)ibin; ++ib)
                    if (ib >= 0 && (beam_z(ib) > Z_missing && beam_w(ib) != W_threshold && ( beam_w(ib) > 0.5 || fabs(beam_v(ib)) > 0.5) ) )
                       count++;
                  if (double(count)/double(min(int(ibin),int(2*min_segment_length))) >=0.25) before = true;

              /* Cerco dati validi in Z dopo il segmento */
                  count = 0;
	          for (unsigned ia = ibin + 1; ia <= ibin + 2*min_segment_length; ++ia)
                   if (ia < beam_size && (beam_z(ia) > Z_missing && (beam_w(ia) != W_threshold && ( beam_w(ia) > 0.5 || fabs(beam_v(ia)) > 0.5))  ))
                        count ++;
                  if (double(count)/double(min(int(beam_size - ibin),int(2*min_segment_length))) >=0.25) after = true;
		}
// printf(" %1d ----- STOP SEGMENT ------ %4d  --  %4d    before %d   after %d ",flag, segment_length,counter, before,after);
                if ((segment_length >= min_segment_length && (!before || !after) ) ||  segment_length >= max_segment_length)
 //               if ((segment_length >= min_segment_length ) ||  segment_length >= max_segment_length)
                {
                    /* qui pulisco */
                    //         printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
                    for (unsigned ib = start; ib <= end; ++ib)
                        res[ib] = true;
                }
            } 
// 	    else printf(" %1d ",flag);

        }
//printf("   %4d %4d \n",counter_clutter,counter_trash);
    }
    return res;
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//  These methods use VRAD, WRAD, sdDBZH, values to clean the beam
//
std::vector<bool> Cleaner::clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& SD, int iray) const
{
    const unsigned beam_size = beam_z.rows();
    vector<bool> res(beam_size, false);
    bool in_a_segment = false;
    unsigned start, end;
    unsigned segment_length;
    bool before, after;
    unsigned counter = 0;

    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
//  printf(" %4d %4d  %6.2f %6.2f %10.6f %6.2f ",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin), beam_sd(ibin));
//printf("     -----    %2x %2x %2x %2x ",(unsigned char)((beam_z(ibin)-scan_z.offset)/scan_z.gain/256),
//(unsigned char)((beam_v(ibin)-scan_v.offset)/scan_v.gain/256),
//(unsigned char)((beam_w(ibin)-scan_w.offset)/scan_w.gain/256),
//(unsigned char)((beam_sd(ibin)-SD.offset)/SD.gain/256));
        if (!in_a_segment)
        {
            /* cerco la prima cella segmento da pulire*/
            if ( ((beam_w(ibin) == W_threshold && beam_v(ibin) == bin_wind_magic_number) ||(beam_w(ibin) * fabs(beam_v(ibin)) <= 0.25) )  && beam_z (ibin) != Z_missing  && beam_sd(ibin) > sd_threshold )
            {
// printf(" 1  ----- START SEGMENT ------");
                in_a_segment = true;
                start = ibin;
                after = false;
                before = false;
            } 
// 	    else printf(" 0 ");
        } else {
            /* cerco la fine segmento da pulire*/
            if ( ( ( beam_w(ibin) != W_threshold || beam_v(ibin) != bin_wind_magic_number) && (beam_w(ibin) * fabs(beam_v(ibin)) > 0.25) ) || ibin == (beam_size - 1) || beam_z(ibin) == Z_missing ||   beam_sd(ibin) <= sd_threshold ) 
            {
                in_a_segment = false;
                end = ibin - 1;
                if (ibin == (beam_size - 1)) end = ibin;  // caso particolare per fine raggio
                /* Fine trovata ora procedo alla pulizia eventuale */
                segment_length = end - start+1;
                counter = counter + (unsigned)(segment_length);

/* 	il segmento è corto allora cerco nei dintorni dei dati validi, se li trovo non pulisco */
		if (segment_length <= 2*min_segment_length ){
              /* Cerco dati validi in Z prima del segmento */
	     	  int count=0;
                  for (int ib = ibin - 2*min_segment_length; ib < (signed)ibin; ++ib)
                    if (ib >= 0 && (beam_z(ib) > Z_missing && beam_w(ib) != W_threshold && ( beam_w(ib) > 0.5 || fabs(beam_v(ib)) > 0.5) ) )
                       count++;
                  if (double(count)/double(min(int(ibin),int(2*min_segment_length))) >=0.25) before = true;

              /* Cerco dati validi in Z dopo il segmento */
                  count = 0;
	          for (unsigned ia = ibin + 1; ia <= ibin + 2*min_segment_length; ++ia)
                   if (ia < beam_size && (beam_z(ia) > Z_missing && (beam_w(ia) != W_threshold && ( beam_w(ia) > 0.5 || fabs(beam_v(ia)) > 0.5))  ))
                        count ++;
                  if (double(count)/double(min(int(beam_size - ibin),int(2*min_segment_length))) >=0.25) after = true;
		}
// printf(" 0 ----- STOP SEGMENT ------ %4d  --  %4d    before %d   after %d ",segment_length,counter, before,after);
                if ((segment_length >= min_segment_length && (!before || !after) ) ||  segment_length >= max_segment_length)
 //               if ((segment_length >= min_segment_length ) ||  segment_length >= max_segment_length)
                {
                    /* qui pulisco */
                    //         printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
                    for (unsigned ib = start; ib <= end; ++ib)
                        res[ib] = true;
                }
            } 
// 	    else printf(" 1 ");

        }
// printf("\n");
    }
    return res;
}


// Senza ZDR - Fuzzy logic
std::vector<unsigned char> Cleaner::eval_classID_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, int iray) const
{

    const unsigned beam_size = beam_z.rows();
    vector<unsigned char> res(beam_size, 0);
    vector<unsigned> counter (6,0) ;

	Matrix2D<double> Wij(6,6);
//		Z	V	Wrad	SD2d	SDRay	SDAz
	Wij <<	1.0,	1.0,	1.0,	0.5,	0.3,	0.3,	// METEO		4.1 
		0.8,	0.5,	0.5,	0.5,	0.5,	0.5,	// CLUTTER		3.3
		0.0,	0.2,	0.2,	0.6,	0.8,	0.5,	// INTERF. Strong	2.3
		0.0,	0.2,	0.2,	0.5,	0.6,	0.4,	// INTERF. Med.		2.2
		0.0,	0.2,	0.2,	0.6,	0.5,	0.5,	// INTERF. Weak		1.6
		0.6,	0.2,	0.2,	0.6,	0.5,	0.5;	// NOISE		2.6
// TOT =	2.4	2.3	2.3	3.4	3.2	2.5



    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
//printf(" %4d %4d  %6.2f %6.2f %10.6f  %10.6f %10.6f %10.6f",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin),beam_sd(ibin),beam_sdray(ibin),beam_sdaz(ibin));
	Matrix2D<double> Pij(6,6);
	Pij = Pij * 0.;
	ArrayXd  Class_WP(6);	
	Class_WP.setZero();
	if (beam_z(ibin)  == Z_missing) {
	  unsigned ID=0;
	  res[ibin]=ID;
	  counter[ID]++;
//printf("%2d %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %4d %4d %4d %4d %4d\n",res[ibin],Class_WP(0),Class_WP(1),Class_WP(2),Class_WP(3),Class_WP(4),Class_WP(5),counter[1],counter[2], counter[3], counter[4], counter[5]);
//printf("%2d %5.3f %4d %4d %4d %4d %4d\n",res[ibin],Class_WP(ID),counter[1],counter[2], counter[3], counter[4], counter[5]);
	  continue;
        }
	    
//eseguo un test unico per VRAD e WRAD e assegno tutte le prob assieme
	if (beam_v(ibin) != bin_wind_magic_number  ) {				//	VRAD
           Pij(0,1)=1.;	 						// METEO		 
	   double prob_v = trap (-0.5, -0.3, 0.3, 0.5, beam_v(ibin));	// CLUTTER		
	   Pij(1,1)=prob_v;	  					// INTERF. Strong	
           Pij(2,1)=prob_v;	  					// INTERF. Med.		
           Pij(3,1)=prob_v;	  					// INTERF. Weak		
           Pij(4,1)=prob_v;	  					// NOISE		
           Pij(5,1)=prob_v;	  					
	} else {
	   Pij(0,1)=0.3;	// METEO		
	   Pij(1,1)=1.;		// CLUTTER		  
           Pij(2,1)=1.;		// INTERF. Strong	  
           Pij(3,1)=1.;		// INTERF. Med.		  
           Pij(4,1)=1.;		// INTERF. Weak		  
           Pij(5,1)=1.;		// NOISE		  
	}

	if (beam_w(ibin) > W_threshold) {					//	WRAD
           Pij(0,2)=1.;	 						// METEO		 
	   double prob_w = trap (0.0001, 0., 0.2, 0.3, beam_w(ibin));	// CLUTTER		
	   Pij(1,2)=prob_w;	  					// INTERF. Strong	
           Pij(2,2)=prob_w;	  					// INTERF. Med.		
           Pij(3,2)=prob_w;	  					// INTERF. Weak		
           Pij(4,2)=prob_w;	  					// NOISE		
           Pij(5,2)=prob_w;	  					
	} else {
	   Pij(0,2)=0.3;	// METEO		
	   Pij(1,2)=1.;		// CLUTTER		  
           Pij(2,2)=1.;		// INTERF. Strong	  
           Pij(3,2)=1.;		// INTERF. Med.		  
           Pij(4,2)=1.;		// INTERF. Weak		  
           Pij(5,2)=1.;		// NOISE		  
	}

// METEO		
	Pij(0,0) = 1.;							//	Z
	Pij(0,3) = trap (0.5, 1., 10.,11., beam_sd(ibin));		//	SD_2D
	Pij(0,4) = trap (0.5, 1., 10.,11., beam_sdray(ibin));		//	SD_RAY
	Pij(0,5) = trap (0.5, 1., 10.,11., beam_sdaz(ibin));		//	SD_AZ

// CLUTTER		
	Pij(1,0) = trap (5., 15., 99., 99.9, beam_z(ibin));		//	Z
	Pij(1,3) = trap (4.5, 5., 99., 99.9, beam_sd(ibin));		//	SD_2D
	Pij(1,4) = trap (4., 4.5, 99., 99.9, beam_sdray(ibin));		//	SD_RAY
	Pij(1,5) = trap (4., 4.5, 99., 99.9, beam_sdaz(ibin));		//	SD_AZ
// INTERF. Strong	
	Pij(2,0) = 1.;							//	Z
	Pij(2,3) = trap (0.5, 1.5, 5., 7., beam_sd(ibin));		//	SD_2D
	Pij(2,4) = trap (0.0, 0.3, 2., 3., beam_sdray(ibin));		//	SD_RAY
	Pij(2,5) = trap (2.0, 4., 90.,90.9, beam_sdaz(ibin));		//	SD_AZ

// INTERF. Med.			
	Pij(3,0) = 1.;							//	Z
	Pij(3,3) = trap (4.0, 5., 90., 90.9, beam_sd(ibin));		//	SD_2D
	Pij(3,4) = trap (0., 0.5, 10., 15.0, beam_sdray(ibin));		//	SD_RAY
	Pij(3,5) = trap (3., 4., 90., 90.9, beam_sdaz(ibin));		//	SD_AZ

// INTERF. Weak		
	Pij(4,0) = 1.;							//	Z
	Pij(4,3) = trap (4., 7., 90., 90.9, beam_sd(ibin));		//	SD_2D
	Pij(4,4) = trap (0.5, 1.5, 3., 4., beam_sdray(ibin));		//	SD_RAY
	Pij(4,5) = trap (3., 4., 90., 90.9, beam_sdaz(ibin));		//	SD_AZ

// NOISE	
	double coeff = 1.;
	if (ibin >= 40 ) {
	   if ( ibin <= 160) coeff =  (160 - ibin)/ 120. ;
	     else coeff = 0.;
        }
	Pij(5,0) = trap(Z_missing-0.0001, Z_missing, 15., 20., beam_z(ibin));	//	Z
	Pij(5,3) = trap (5., 7., 99., 99.9, beam_sd(ibin));		//	SD_2D
	Pij(5,4) = trap (0., 0.001, 0.8 + coeff * 9.2 , 1. + coeff * 14.,beam_sdray(ibin));		//	SD_RAY
//	Pij(5,4) = trap (0., 0.001, 0.8, 1.,beam_sdray(ibin));		//	SD_RAY
//	Pij(5,4) = trap (0., 0.001, 10., 15.,beam_sdray(ibin));		//	SD_RAY
	Pij(5,5) = trap (1.5, 3., 99., 99.9,beam_sdaz(ibin));		//	SD_AZ

//---- fine calcolo probabilità
// Calcolo classe appartenenza

	Class_WP = ((Wij.array()*Pij.array()).matrix()*VectorXd::Ones(6)).array()/(Wij*VectorXd::Ones(6)).array();
	unsigned i,ID;
	Class_WP.maxCoeff(&i);
	ID=i;
	if (Class_WP(i) < 0.1 ) ID=6;
	res[ibin]=ID;
	counter[ID]++;
//printf("%2d %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %4d %4d %4d %4d %4d\n",res[ibin],Class_WP(0),Class_WP(1),Class_WP(2),Class_WP(3),Class_WP(4),Class_WP(5),counter[1],counter[2], counter[3], counter[4], counter[5]);
//	if (iray == 225 and ibin <= 50) {
//	std::cout <<"Valori per fuzzy - Z_missing "<<Z_missing<< "differenza " <<beam_z(ibin)-Z_missing << std::endl;
//	std::cout << Pij<<std::endl;
//        }
    }

    return res;
}

// Senza ZDR - Basato su test
std::vector<unsigned char> Cleaner::eval_clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, int iray) const
{

    const unsigned beam_size = beam_z.rows();
    vector<unsigned char> res(beam_size, 0);
    unsigned counter = 0;
    unsigned countIntWeak = 0;
    unsigned countIntStrong = 0;
    unsigned countIntMed = 0;
    unsigned countClutter = 0;
    unsigned countNoise = 0;

    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
printf(" %4d %4d  %6.2f %6.2f %10.6f  %10.6f %10.6f %10.6f",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin),beam_sd(ibin),beam_sdray(ibin),beam_sdaz(ibin));
	if ((beam_w(ibin) > W_threshold && beam_v(ibin) != bin_wind_magic_number && beam_sd(ibin) >= 1. && 
	    beam_sd(ibin) <= 10. ) || beam_z(ibin) == Z_missing) {
	// this should be a meteorological echo
	;
        }else {
	  res [ibin]=1;	// Not meteo but unclassified
	  counter++;
          if (beam_w(ibin) == W_threshold && beam_v(ibin) == bin_wind_magic_number)
          {
	    if (beam_sd(ibin) >= 5. && beam_sdray(ibin) >= 2 && beam_sdaz(ibin) > 4.){ 
		// this should be clutter
		res[ibin] = 2;
	        countClutter ++;
            }
	    if (beam_sd(ibin) >= 1. && beam_sd(ibin) <= 5. && beam_sdray(ibin) < 2. && beam_sdaz(ibin) > 4.){ 
		// this should be a strong Interference
		res[ibin] = 3;
	        countIntStrong ++;
            }
	    if (beam_sd(ibin) >= 3. && beam_sd(ibin) <= 7. && beam_sdray(ibin) < 2. && beam_sdaz(ibin) > 4.&& beam_sdaz(ibin) <7 ){ 
		// this should be a medium Interference
		res[ibin] = 4;
	        countIntMed ++;
            }
	    //if (beam_sd(ibin) >= 5. && beam_sd(ibin) <= 3. && beam_sdray(ibin) < 3. && beam_sdaz(ibin) < 3.){ 
	    if (beam_sd(ibin) >= 5. && beam_sdray(ibin) < 3. && beam_sdaz(ibin) < 3.){ 
		// this should be a weakInterference
		res[ibin] = 5;
	        countIntWeak ++;
            }
	    if (beam_sd(ibin) >= 10. && beam_sdray(ibin) < 2. && beam_sdaz(ibin) < 5. && beam_z(ibin) <10.){ 
		// this should be a noise
		res[ibin] = 6;
	        countNoise ++;
            }
        }
        }   // ELSE this should be not a meteo echo
printf("%2d %4d %4d %4d %4d %4d\n",res[ibin],countClutter,countIntStrong, countIntMed, countIntWeak, countNoise);
    }

    return res;

}
//----------------------------------------------------------------------------------



void Cleaner::evaluateCleanID(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v,PolarScan<unsigned char>& scan_cleanID, unsigned iel)
{
    //return evaluateCleanID(scan_z, scan_w, scan_v, scan_cleanID, scan_v.undetect,iel);
    return evaluateClassID(scan_z, scan_w, scan_v, scan_cleanID, scan_v.undetect,iel);
}

void Cleaner::evaluateCleanID(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v,PolarScan<unsigned char>& scan_cleanID, double bin_wind_magic_number, unsigned iel)
{

    if (scan_z.beam_count != scan_w.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_w beam_count");
    if (scan_z.beam_size != scan_w.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_w beam_size");

    if (scan_z.beam_count != scan_v.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_v beam_count");
    if (scan_z.beam_size != scan_v.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_v beam_size");

    Cleaner cleaner(scan_z.undetect, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);

    const unsigned beam_count = scan_z.beam_count;
    const unsigned beam_size = scan_z.beam_size;

  //  fprintf(stderr, "NEWCLEANER zmis %f, wthr %f, vmis %f, mn %f\n",
  //          cleaner.Z_missing, cleaner.W_threshold, cleaner.V_missing, cleaner.bin_wind_magic_number);

//    radarelab::volume::Scans<double>   Z_S,  SD2D;
//    Z_S.push_back(scan_z);
//    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);

    radarelab::volume::Scans<double>   Z_S,  SD2D,SD_Ray,SD_Az;
    Z_S.push_back(scan_z);
    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);

    radarelab::volume::textureSD( Z_S,SD_Ray, scan_z.cell_size*9 , 360./scan_z.beam_count,true);
    radarelab::volume::textureSD( Z_S,SD_Az, scan_z.cell_size , 5*360./scan_z.beam_count,true);

    for (unsigned i = 0; i < beam_count; ++i)
    {
        //vector<unsigned char> corrected = cleaner.eval_clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),i);
        vector<unsigned char> corrected = cleaner.eval_clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),SD2D[0].row(i), SD_Ray[0].row(i), SD_Az[0].row(i), i);
        //vector<bool> corrected = cleaner.clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),SD2D[0].row(i), SDZDR2D[0].row(i), scan_z, scan_w, scan_v, SD2D[0],i);
        for (unsigned ib = 0; ib < beam_size; ++ib)
	   scan_cleanID(i,ib)=corrected[ib];
    }
}

void Cleaner::evaluateClassID(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v,PolarScan<unsigned char>& scan_cleanID, double bin_wind_magic_number, unsigned iel)
{

    if (scan_z.beam_count != scan_w.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_w beam_count");
    if (scan_z.beam_size != scan_w.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_w beam_size");

    if (scan_z.beam_count != scan_v.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_v beam_count");
    if (scan_z.beam_size != scan_v.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_v beam_size");

    Cleaner cleaner(scan_z.undetect, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);

    const unsigned beam_count = scan_z.beam_count;
    const unsigned beam_size = scan_z.beam_size;

// compute texture volumes
    radarelab::volume::Scans<double>   Z_S,  SD2D,SD_Ray,SD_Az;
    Z_S.push_back(scan_z);
    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);
    radarelab::volume::textureSD( Z_S,SD_Ray, scan_z.cell_size*21 , 360./scan_z.beam_count,true);
    radarelab::volume::textureSD( Z_S,SD_Az, scan_z.cell_size , 5*360./scan_z.beam_count,true);


    for (unsigned i = 0; i < beam_count; ++i)
    {
        vector<unsigned char> corrected = cleaner.eval_classID_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),SD2D[0].row(i), SD_Ray[0].row(i), SD_Az[0].row(i), i);
        for (unsigned ib = 0; ib < beam_size; ++ib)
	   scan_cleanID(i,ib)=corrected[ib];
    }
}

void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, unsigned iel , bool set_undetect)
{
    return clean(scan_z, scan_w, scan_v, scan_v.undetect,iel);
}
void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v,  double bin_wind_magic_number, unsigned iel, bool set_undetect)
{

    if (scan_z.beam_count != scan_w.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_w beam_count");
    if (scan_z.beam_size != scan_w.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_w beam_size");

    if (scan_z.beam_count != scan_v.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_v beam_count");
    if (scan_z.beam_size != scan_v.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_v beam_size");

    Cleaner cleaner(scan_z.undetect, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);

    const unsigned beam_count = scan_z.beam_count;
    const unsigned beam_size = scan_z.beam_size;

  //  fprintf(stderr, "NEWCLEANER zmis %f, wthr %f, vmis %f, mn %f\n",
  //          cleaner.Z_missing, cleaner.W_threshold, cleaner.V_missing, cleaner.bin_wind_magic_number);

    radarelab::volume::Scans<double>   Z_S,  SD2D,SD_Ray,SD_Az;
    Z_S.push_back(scan_z);
    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);
    radarelab::volume::textureSD( Z_S,SD_Ray, scan_z.cell_size*9 , 360./scan_z.beam_count,false);
    radarelab::volume::textureSD( Z_S,SD_Az, scan_z.cell_size , 5*360./scan_z.beam_count,false);

//radarelab::gdal_init_once();
//printf("scrivo Z ");
//Matrix2D <double>img;
//img = (scan_z.array() - scan_z.offset )/ scan_z.gain /256 ;
//Matrix2D <unsigned char>img_tmp, z_clean;
//std::string ext;
//char pippo[200];
//sprintf(pippo, "_%02d.png",iel);
//ext=pippo;

//img_tmp=img.cast<unsigned char>();
//z_clean=img_tmp;
//radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_Z"+ext,  "PNG");

//printf("V ");
//img = (scan_v.array()-scan_v.offset)/scan_v.gain/256 ;
//img_tmp=img.cast<unsigned char>();
//radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_V"+ext,"PNG");
//printf("W ");
//img = (scan_w.array()-scan_w.offset)/scan_w.gain/256 ;
//img_tmp=img.cast<unsigned char>();
//radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_W"+ext,"PNG");
//printf("SD2d ");
//img = (SD2D[0].array()-SD2D[0].offset)/SD2D[0].gain/256 ;
//img_tmp=img.cast<unsigned char>();
//radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_SD2d"+ext,"PNG");
//printf("\n");

    double new_val=cleaner.Z_missing;
    if (set_undetect) new_val=scan_z.undetect;

    for (unsigned i = 0; i < beam_count; ++i)
    {
        // Compute which elements need to be cleaned
       // vector<bool> corrected = cleaner.clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),SD2D[0].row(i), scan_z, scan_w, scan_v, SD2D[0],i);
       vector<bool> corrected = cleaner.clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),i);

        for (unsigned ib = 0; ib < beam_size; ++ib)
            if (corrected[ib])
            {
                //scan_z(i, ib) = cleaner.Z_missing;
                scan_z(i, ib) = new_val;
    //            scan_w(i, ib) = cleaner.W_threshold;
    //            scan_v(i, ib) = cleaner.V_missing;
	    }
//	       img_tmp(i,ib)=255;
//	       z_clean(i,ib)=0;
//            } else img_tmp(i,ib)= 0 ;

    }
//radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_clean"+ext,"PNG");
//radarelab::write_image(z_clean,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_Zclean"+ext,"PNG");
}



void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& scan_zdr, unsigned iel, bool set_undetect )
{
    return clean(scan_z, scan_w, scan_v, scan_zdr, scan_v.undetect,iel);
}
void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& scan_zdr, double bin_wind_magic_number, unsigned iel, bool set_undetect )
{
    if (scan_z.beam_count != scan_w.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_w beam_count");
    if (scan_z.beam_size != scan_w.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_w beam_size");

    if (scan_z.beam_count != scan_v.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_v beam_count");
    if (scan_z.beam_size != scan_v.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_v beam_size");

    Cleaner cleaner(scan_z.undetect, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);

    const unsigned beam_count = scan_z.beam_count;
    const unsigned beam_size = scan_z.beam_size;

  //  fprintf(stderr, "NEWCLEANER zmis %f, wthr %f, vmis %f, mn %f\n",
  //          cleaner.Z_missing, cleaner.W_threshold, cleaner.V_missing, cleaner.bin_wind_magic_number);

    radarelab::volume::Scans<double>   Z_S,  SD2D;
    Z_S.push_back(scan_z);
    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);
    radarelab::volume::Scans<double>   ZDR_S,  SDZDR2D;
    ZDR_S.push_back(scan_zdr);
    radarelab::volume::textureSD( ZDR_S,SDZDR2D, 1000. , 3,false);

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// Mettere a true per fare grafica per debug o false per non fare grafica
//
// RICORDARSI DI TOGLIERE/METTERE COMMENTI DOPO CLEAN_BEAM
// -------------------------------------------------------------------
Matrix2D <unsigned char>img_tmp, z_clean;
Matrix2D <double>img;
std::string ext;
char pippo[200];
if (false){
  radarelab::gdal_init_once();
      
  printf("scrivo Z ");
  img = (scan_z.array() - scan_z.offset )/ scan_z.gain /256 ;
  sprintf(pippo, "_%02d.png",iel);
  ext=pippo;
  img_tmp=img.cast<unsigned char>();
  z_clean=img_tmp;
  radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_Z"+ext,  "PNG");

//  printf("V ");
//  img = (scan_v.array()-scan_v.offset)/scan_v.gain/256 ;
//  img_tmp=img.cast<unsigned char>();
//  radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_V"+ext,"PNG");
//  printf("W ");
//  img = (scan_w.array()-scan_w.offset)/scan_w.gain/256 ;
//  img_tmp=img.cast<unsigned char>();
//  radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_W"+ext,"PNG");
  printf("SD2d ");
  img = (SD2D[0].array()-SD2D[0].offset)/SD2D[0].gain/256 ;
  img_tmp=img.cast<unsigned char>();
  radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_SD2d"+ext,"PNG");
  printf("SDZDR2d ");
  img = (SDZDR2D[0].array()-SDZDR2D[0].offset)/SDZDR2D[0].gain/256 ;
  img_tmp=img.cast<unsigned char>();
  radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_SDZDR2d"+ext,"PNG");
  printf("\n");
}
    double new_val=cleaner.Z_missing;
    if (set_undetect) new_val=scan_z.undetect;

    for (unsigned i = 0; i < beam_count; ++i)
    {
        // Compute which elements need to be cleaned
        vector<bool> corrected = cleaner.clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i),SD2D[0].row(i), SDZDR2D[0].row(i), scan_z, scan_w, scan_v, SD2D[0],i);

        for (unsigned ib = 0; ib < beam_size; ++ib)
            if (corrected[ib])
            {
                //scan_z(i, ib) = cleaner.Z_missing;
                scan_z(i, ib) = new_val;
  //              scan_w(i, ib) = cleaner.W_threshold;
    //            scan_v(i, ib) = cleaner.V_missing;
            }
//	       img_tmp(i,ib)=255;
//	       z_clean(i,ib)=0;
//            } else img_tmp(i,ib)= 0 ;

    }
//radarelab::write_image(img_tmp,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_clean"+ext,"PNG");
//radarelab::write_image(z_clean,"/ponte/rad_svn/proc_operative/test_arch/rev_actual/radar/immagini/Cleaner/PPI_Zclean"+ext,"PNG");
}

void Cleaner::clean( radarelab::volume::Loader load_structure, double bin_wind_magic_number,unsigned iel, bool set_undetect)
{
  std::string Z_Quantity;

}


double Cleaner::trap(double x1, double x2, double x3, double x4, double val) const
{
	if(val<=x3&&val>=x2) return 1.;
	else if(val<x2&&val>x1) return val/(x2-x1)-x1/(x2-x1);
	else if (val<x4&&val>x3) return val/(x3-x4)-x4/(x3-x4);
	else return 0.; // (val<=x1||val>=x4)
}


}
}
