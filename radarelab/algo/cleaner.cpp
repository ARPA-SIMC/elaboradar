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
if(ibin <40 && false){
	printf(" %4d %4d  %6.2f %6.2f %10.6f %6.2f %6.2f ",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin), beam_sd(ibin),beam_sdzdr(ibin));
	printf("     -----    %2x %2x %2x %2x ",(unsigned char)((beam_z(ibin)-scan_z.offset)/scan_z.gain/256),
	(unsigned char)((beam_v(ibin)-scan_v.offset)/scan_v.gain/256),
	(unsigned char)((beam_w(ibin)-scan_w.offset)/scan_w.gain/256),
	(unsigned char)((beam_sd(ibin)-SD.offset)/SD.gain/256));
}       
        if (!in_a_segment)
        {
            /* cerco la prima cella segmento da pulire*/
            if ( is_clutter || is_trash  )
            {
// if(ibin <40)printf(" %1d  ----- START SEGMENT ------",flag);
               in_a_segment = true;
                start = ibin;
                after = false;
                before = false;
            } 
//            else  if(ibin <40)printf(" %1d ",flag);
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
//  if(ibin <40)printf(" %1d ----- STOP SEGMENT ------ %4d  --  %4d    before %d   after %d ",flag, segment_length,counter, before,after);
                if ((segment_length >= min_segment_length && (!before || !after) ) ||  segment_length >= max_segment_length)
 //               if ((segment_length >= min_segment_length ) ||  segment_length >= max_segment_length)
                {
                    /* qui pulisco */
//                              if(ibin <40)printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
                    for (unsigned ib = start; ib <= end; ++ib)
                        res[ib] = true;
                }
            } 
// 	    else  if(ibin <40)printf(" %1d ",flag);

        }
// if(ibin <40)printf("   %4d %4d \n",counter_clutter,counter_trash);
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

// CC: fuzzy logic
std::vector<unsigned char> Cleaner::eval_classID_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_zdr, const Eigen::VectorXd& beam_rohv, const Eigen::VectorXd& beam_sqi, const Eigen::VectorXd& beam_snr, const Eigen::VectorXd& beam_zvd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, const Eigen::VectorXd& beam_zdr_sd, int iray, const string radar, double v_ny, bool stamp, bool force_meteo) const
{

    const unsigned beam_size = beam_z.rows();
    vector<unsigned char> res(beam_size, 0);
    vector<double> diffs(beam_size, 0);
    unsigned prevID = 4;
    int Num_entries=0;
    int Num_echoes = 5;
    int Ntraps = 5; // 5 argomenti da passare a Trap : x1,x2,x3,x4,x5

    //leggo matrice dei pesi----------------------------------------------------
    string fin_w = "./matrix-"+radar+".txt";
    vector<string> w_vector;
    w_vector = read_matrix_from_txt(fin_w);
    Num_entries = w_vector.size()/Num_echoes;

    Matrix2D<double> Wij(Num_echoes,Num_entries);
    for(int i=0;i<Num_echoes;i++){ //itero colonna
      for(int j=0;j<Num_entries;j++){ //itero rriga
        Wij(i,j) = stod( w_vector[i*Num_entries+j]);
	//cout<<" W["<<i<<","<<j<<"]="<<Wij(i,j);
      }
      //cout<<" "<<endl;
    }

    //--------------------------------------------------------------------------
    //leggo matrice dei traps
    string fin_t = "./Trap-"+radar+".txt";
    vector<string> t_vector;
    t_vector = read_matrix_from_txt(fin_t);
    double Traps[Num_entries][Num_echoes][Ntraps];
    for(int i=0;i<Num_entries;i++){
      for(int j=0;j<Num_echoes;j++){ 
	for(int k=0;k<Ntraps;k++){
	  Traps[i][j][k] = stod( t_vector[i*(Num_echoes*Ntraps)+j*Num_echoes+k]);
	  //cout<<" Traps["<<i<<","<<j<<","<<k<<"]="<<Traps[i][j][k];
	}
	//cout<<endl;
      }
      //cout<<endl;
    }
    //------------------------------------------------------------------------
    
    vector<unsigned> counter (Num_entries,0) ; // non sono sicura di cosa delle dimensioni di questo counter
    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    //for (unsigned ibin = 87*4; ibin < 87*4+12; ++ibin)
    {
      //printf(" %4d %4d  %6.2f %6.2f %10.6f %10.6f  %10.6f %10.6f %10.6f",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin),beam_zdr(ibin),beam_sd(ibin),beam_sdray(ibin),beam_sdaz(ibin));
	Matrix2D<double> Pij(Num_echoes,Num_entries);
	Pij = Pij * 0.;
	//cout<<Pij<<endl;
	ArrayXd  Class_WP(Num_echoes);	// la dimensione di Class_WP deve essere Num_echoes, perchè alla fine ricavo un vettore con 6 valori, uno per ogni echo, passando per prodotto pesi prob 
	Class_WP.setZero();
	if (beam_z(ibin)  == Z_missing) {
	  unsigned ID=0;
	  res[ibin]=ID;
	  counter[ID]++;

	  continue;
        }

	//cout<<"V Nyquist: "<<v_ny<<endl;
	//cout<<sizeof(Traps)<<" "<<Wij.size()<<" "<<Pij.size()<<endl;
//eseguo un test unico per VRAD e WRAD e assegno tutte le prob assieme
	if (beam_v(ibin) != bin_wind_magic_number  ) {				//	VRAD
	  Pij(0,1)=trap(v_ny*Traps[0][0][0], v_ny*Traps[0][0][1], v_ny*Traps[0][0][2], v_ny*Traps[0][0][3], beam_v(ibin),v_ny*Traps[0][0][4]);	 						// METEO		 
	   double prob_v = trap (v_ny, v_ny, v_ny, v_ny, beam_v(ibin));	
	   //cout<<"prob_v computed: "<<prob_v<<endl;
	   //for(int e=1;e<Num_echoes;e++){ Pij(e,1) = prob_v };		
	   Pij(1,1)=trap(Traps[0][1][0],Traps[0][1][1],Traps[0][1][2],Traps[0][1][3], beam_v(ibin),v_ny*Traps[0][1][4]);	  					// CLUTTER	
           Pij(2,1)=prob_v;	  					// INTERF. Strong		
           Pij(3,1)=prob_v;	  					// INTERF. Med.		
           Pij(4,1)=prob_v;	  					// INTERF. Weak		
	} else {
	   Pij(0,1)=0.3;	// METEO		
	   Pij(1,1)=1.;		// CLUTTER		  
           Pij(2,1)=1.;		// INTERF. Strong	  
           Pij(3,1)=1.;		// INTERF. Med.		  
           Pij(4,1)=1.;		// INTERF. Weak		  
	}

	// WRAD
	Pij(0,2)= trap(Traps[1][0][0],Traps[1][0][1],Traps[1][0][2],Traps[1][0][3],beam_w(ibin));                   // METEO
	Pij(1,2)= trap(Traps[1][1][0],Traps[1][1][1],Traps[1][1][2],Traps[1][1][3],beam_w(ibin));                  // CLUTTER
	double prob_w = trap(Traps[1][2][0],Traps[1][2][1],Traps[1][2][2],Traps[1][2][3],beam_w(ibin));
	Pij(2,2) = prob_w;                                           // INTERF MULTIPLE
	Pij(3,2) = prob_w;                                           // INTERF. Med.
	Pij(4,2) = prob_w;                                           // INTERF. Weak
	
// METEO		
	Pij(0,0) = trap(Traps[2][0][0],Traps[2][0][1],Traps[2][0][2],Traps[2][0][3],beam_z(ibin), Traps[2][0][4]);				//	Z
	Pij(0,3) = trap (Traps[3][0][0],Traps[3][0][1],Traps[3][0][2],Traps[3][0][3], beam_sd(ibin));		//	SD_2D
	Pij(0,4) = trap (Traps[4][0][0],Traps[4][0][1],Traps[4][0][2],Traps[4][0][3], beam_sdray(ibin));		//	SD_RAY
	Pij(0,5) = trap (Traps[5][0][0],Traps[5][0][1],Traps[5][0][2],Traps[5][0][3], beam_sdaz(ibin));		//	SD_AZ	
        //if(Num_entries>6)
	Pij(0,6) = trap(Traps[6][0][0],Traps[6][0][1],Traps[6][0][2],Traps[6][0][3],beam_zdr(ibin));                 //      ZDR
	Pij(0,7) = trap(Traps[7][0][0],Traps[7][0][1],Traps[7][0][2],Traps[7][0][3],beam_rohv(ibin),Traps[7][0][4]);              //      ROHV
	Pij(0,8) = trap(Traps[8][0][0],Traps[8][0][1],Traps[8][0][2],Traps[8][0][3],beam_zdr_sd(ibin));              //      ZDR_SD2D
	Pij(0,9) = trap(Traps[9][0][0],Traps[9][0][1],Traps[9][0][2],Traps[9][0][3], beam_sqi(ibin));           // SQI
	Pij(0,10) = trap(Traps[10][0][0],Traps[10][0][1],Traps[10][0][2],Traps[10][0][3], beam_snr(ibin));           // SNR
	Pij(0,11) = trap(Traps[11][0][0],Traps[11][0][1],Traps[11][0][2],Traps[11][0][3],beam_zvd(ibin));          // DBZH_VD
	//cout<<"Num_entries control active"<<endl;

// CLUTTER		
	Pij(1,0) = trap (Traps[2][1][0],Traps[2][1][1],Traps[2][1][2],Traps[2][1][3], beam_z(ibin), Traps[2][1][4]);		//	Z
	Pij(1,3) = trap (Traps[3][1][0],Traps[3][1][1],Traps[3][1][2],Traps[3][1][3], beam_sd(ibin));		//	SD_2D
	Pij(1,4) = trap (Traps[4][1][0],Traps[4][1][1],Traps[4][1][2],Traps[4][1][3], beam_sdray(ibin));	//	SD_RAY
	Pij(1,5) = trap (Traps[5][1][0],Traps[5][1][1],Traps[5][1][2],Traps[5][1][3], beam_sdaz(ibin));		//	SD_AZ
	//if(Num_entries>6)
	Pij(1,6) = trap(Traps[6][1][0],Traps[6][1][1],Traps[6][1][2],Traps[6][1][3],beam_zdr(ibin),Traps[6][1][4] );            //      ZDR 
	Pij(1,7) = trap(Traps[7][1][0],Traps[7][1][1],Traps[7][1][2],Traps[7][1][3],beam_rohv(ibin),Traps[7][1][4]);                 //      ROHV
	Pij(1,8) = trap(Traps[8][1][0],Traps[8][1][1],Traps[8][1][2],Traps[8][1][3],beam_zdr_sd(ibin));              //      ZDR_SD2D
	Pij(1,9) = trap(Traps[9][1][0],Traps[9][1][1],Traps[9][1][2],Traps[9][1][3],beam_sqi(ibin));             // SQI
	Pij(1,10) = trap(Traps[10][1][0],Traps[10][1][1],Traps[10][1][2],Traps[10][1][3], beam_snr(ibin));           // SNR
	Pij(1,11) = trap(Traps[11][1][0],Traps[11][1][1],Traps[11][1][2],Traps[11][1][3],beam_zvd(ibin));//, trap(70.,70.,100.,100.,beam_zvd(ibin)) );          // DBZH_VD
	//if(radar=="GAT"&&beam_zvd(ibin)>=70.) Pij(1,11) = 1.0 ;
// INTERF. Strong	
	Pij(2,0) = trap(Traps[2][2][0],Traps[2][2][1],Traps[2][2][2],Traps[2][2][3],beam_z(ibin), Traps[2][2][4]);		//	Z
	Pij(2,3) = trap (Traps[3][2][0],Traps[3][2][1],Traps[3][2][2],Traps[3][2][3], beam_sd(ibin));		//	SD_2D
	Pij(2,4) = trap (Traps[4][2][0],Traps[4][2][1],Traps[4][2][2],Traps[4][2][3], beam_sdray(ibin));		//	SD_RAY
	Pij(2,5) = trap (Traps[5][2][0],Traps[5][2][1],Traps[5][2][2],Traps[5][2][3], beam_sdaz(ibin));		//	SD_AZ
	//if(Num_entries>6)
	Pij(2,6) = trap(Traps[6][2][0],Traps[6][2][1],Traps[6][2][2],Traps[6][2][3],beam_zdr(ibin),Traps[6][2][4]);           //      ZDR
	Pij(2,7) = trap(Traps[7][2][0],Traps[7][2][1],Traps[7][2][2],Traps[7][2][3],beam_rohv(ibin));                 //      ROHV
	Pij(2,8) = trap(Traps[8][2][0],Traps[8][2][1],Traps[8][2][2],Traps[8][2][3],beam_zdr_sd(ibin));              //      ZDR_SD2D
	Pij(2,9) = trap(Traps[9][2][0],Traps[9][2][1],Traps[9][2][2],Traps[9][2][3],beam_sqi(ibin),Traps[9][2][4]);             // SQI
        Pij(2,10) = trap(Traps[10][2][0],Traps[10][2][1],Traps[10][2][2],Traps[10][2][3], beam_snr(ibin));           // SNR
	Pij(2,11) = trap(Traps[11][2][0],Traps[11][2][1],Traps[11][2][2],Traps[11][2][3],beam_zvd(ibin));          // DBZH_VD
// INTERF. Med.			
	Pij(3,0) = trap(Traps[2][3][0],Traps[2][3][1],Traps[2][3][2],Traps[2][3][3],beam_z(ibin), Traps[2][3][4]);		//	Z
	Pij(3,3) = trap (Traps[3][3][0],Traps[3][3][1],Traps[3][3][2],Traps[3][3][3], beam_sd(ibin));		//	SD_2D
	Pij(3,4) = trap (Traps[4][3][0],Traps[4][3][1],Traps[4][3][2],Traps[4][3][3], beam_sdray(ibin));		//	SD_RAY
	Pij(3,5) = trap (Traps[5][3][0],Traps[5][3][1],Traps[5][3][2],Traps[5][3][3], beam_sdaz(ibin));		//	SD_AZ
	//if(Num_entries>6)
	Pij(3,6) = trap(Traps[6][3][0],Traps[6][3][1],Traps[6][3][2],Traps[6][3][3],beam_zdr(ibin),Traps[6][3][4]);           //      ZDR
	Pij(3,7) = trap(Traps[7][3][0],Traps[7][3][1],Traps[7][3][2],Traps[7][3][3],beam_rohv(ibin));                 //      ROHV
	Pij(3,8) = trap(Traps[8][3][0],Traps[8][3][1],Traps[8][3][2],Traps[8][3][3],beam_zdr_sd(ibin));             //      ZDR_SD2D
        Pij(3,9) = trap(Traps[9][3][0],Traps[9][3][1],Traps[9][3][2],Traps[9][3][3],beam_sqi(ibin),Traps[9][3][4]);             // SQI
	Pij(3,10) = trap(Traps[10][3][0],Traps[10][3][1],Traps[10][3][2],Traps[10][3][3], beam_snr(ibin));           // SNR
	Pij(3,11) = trap(Traps[11][3][0],Traps[11][3][1],Traps[11][3][2],Traps[11][3][3],beam_zvd(ibin));          // DBZH_VD
// INTERF. Weak		
	Pij(4,0) = trap(Traps[2][4][0],Traps[2][4][1],Traps[2][4][2],Traps[2][4][3],beam_z(ibin), Traps[2][4][4]);		//	Z
	Pij(4,3) = trap (Traps[3][4][0],Traps[3][4][1],Traps[3][4][2],Traps[3][4][3], beam_sd(ibin));		//	SD_2D
	Pij(4,4) = trap (Traps[4][4][0],Traps[4][4][1],Traps[4][4][2],Traps[4][4][3], beam_sdray(ibin));		//	SD_RAY
	Pij(4,5) = trap (Traps[5][4][0],Traps[5][4][1],Traps[5][4][2],Traps[5][4][3], beam_sdaz(ibin));		//	SD_AZ
	//if(Num_entries>6)
	Pij(4,6) = trap(Traps[6][4][0],Traps[6][4][1],Traps[6][4][2],Traps[6][4][3],beam_zdr(ibin),Traps[6][4][4]);            //      ZDR
	Pij(4,7) = trap(Traps[7][4][0],Traps[7][4][1],Traps[7][4][2],Traps[7][4][3],beam_rohv(ibin),Traps[7][4][4]);                 //      ROHV
	Pij(4,8) = trap(Traps[8][4][0],Traps[8][4][1],Traps[8][4][2],Traps[8][4][3],beam_zdr_sd(ibin));              //      ZDR_SD2D
        Pij(4,9) = trap(Traps[9][4][0],Traps[9][4][1],Traps[9][4][2],Traps[9][4][3],beam_sqi(ibin),Traps[9][4][4]);             // SQI
	Pij(4,10) = trap(Traps[10][4][0],Traps[10][4][1],Traps[10][4][2],Traps[10][4][3], beam_snr(ibin));           // SNR
	Pij(4,11) = trap(Traps[11][4][0],Traps[11][4][1],Traps[11][4][2],Traps[11][4][3],beam_zvd(ibin));          // DBZH_VD

//---- fine calcolo probabilità
// Calcolo classe appartenenza
	Class_WP = ((Wij.array()*Pij.array()).matrix()*VectorXd::Ones(Num_entries)).array()/(Wij*VectorXd::Ones(Num_entries)).array();
	unsigned i,ID;
	Class_WP.maxCoeff(&i);
	ID=i;
	if (Class_WP(i) < 0.1 ) ID=5;
	res[ibin]=ID;
	double diff = Class_WP[i]-Class_WP[0];
	diff = diff*100;
	//printf("ID %d \n",ID);
	//counter[ID]++;
	if(stamp){
	  printf("bin %4d ",ibin);
	  //printf("ID  %d -- diff=%8.3f\n",ID, diff);
	  //printf("prevID %")
	  //lo aggiungi come array
	  for(unsigned c=0;c<5;c++) printf("%5.3f ",Class_WP(c));
	    printf(" ID %d \n",ID);
	    printf(" %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",beam_z(ibin), beam_v(ibin), beam_w(ibin), beam_sd(ibin), beam_sdray(ibin), beam_sdaz(ibin), beam_zdr(ibin),
	      beam_rohv(ibin), beam_zdr_sd(ibin), beam_sqi(ibin), beam_snr(ibin), beam_zvd(ibin));
	  for (unsigned c=0;c<5; c++){
	    for (unsigned d=0; d<12; d++) printf(" %8.3f", Pij(c,d));
	    printf("\n");
	  }
	}
	//if((force_meteo)&&(diff<=50.0)&&(prevID==0)&&(ID==4)) res[ibin] = 0;
	if(radar=="GAT"){
	  if((force_meteo)&&(diff<=10.0)&&(prevID==0)&&(ID!=1)) res[ibin] = 0;
	}
	else if(radar=="SPC"){
	  if((force_meteo)&&(diff<=10.)&&(prevID==0)) res[ibin] = 0;
	
	  if((force_meteo)&&(diff<=20.0)&&(prevID==0)&&(ID==4)) res[ibin] = 0;
	}else cout<<"radar non specificato"<<endl;
	
	diffs[ibin]=diff;
	counter[ID]++;
	prevID = res[ibin];
    }

    return res;//, diffs};
	
    }

// CC: fuzzy logic senza zdr
std::vector<unsigned char> Cleaner::eval_classID_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v, const Eigen::VectorXd& beam_sd, const Eigen::VectorXd& beam_sdray, const Eigen::VectorXd& beam_sdaz, int iray, const string radar, double v_ny) const
{

    const unsigned beam_size = beam_z.rows();
    vector<unsigned char> res(beam_size, 0);
    int Num_entries=0;
    int Num_echoes = 5;
    int Ntraps = 5; // 5 argomenti da passare a Trap : x1,x2,x3,x4,x5

    //leggo matrice dei pesi
    string fin_w = "./matrix-"+radar+"-nozdr.txt";
    vector<string> w_vector;
    w_vector = read_matrix_from_txt(fin_w);
    Num_entries = w_vector.size()/Num_echoes;

    Matrix2D<double> Wij(Num_echoes,Num_entries);
    for(int i=0;i<Num_echoes;i++){ //itero colonna
      for(int j=0;j<Num_entries;j++){ //itero rriga
        Wij(i,j) = stod( w_vector[i*Num_entries+j]);
      }
    }

    //leggo matrice dei traps
    string fin_t = "./Trap-SPC-nozdr.txt";
    vector<string> t_vector;
    t_vector = read_matrix_from_txt(fin_t);
    double Traps[Num_entries][Num_echoes][Ntraps];
    for(int i=0;i<Num_entries;i++){
      for(int j=0;j<Num_echoes;j++){ 
	for(int k=0;k<Ntraps;k++){
	  Traps[i][j][k] = stod( t_vector[i*(Num_echoes*Ntraps)+j*Num_echoes+k]);
	  //cout<<" Traps["<<i<<","<<j<<","<<k<<"]="<<Traps[i][j][k];
	}
	//cout<<endl;
      }
      //cout<<endl;
    }
    //------------------------------------------------------------------------
    
    vector<unsigned> counter (Num_entries,0) ; // non sono sicura di cosa delle dimensioni di questo counter
    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    //for (unsigned ibin = 87*4; ibin < 87*4+12; ++ibin)
    {
      //printf(" %4d %4d  %6.2f %6.2f %10.6f %10.6f  %10.6f %10.6f %10.6f",iray,ibin , beam_z(ibin),beam_v(ibin),beam_w(ibin),beam_zdr(ibin),beam_sd(ibin),beam_sdray(ibin),beam_sdaz(ibin));
	Matrix2D<double> Pij(Num_echoes,Num_entries);
	Pij = Pij * 0.;
	//cout<<Pij<<endl;
	ArrayXd  Class_WP(Num_echoes);	// la dimensione di Class_WP deve essere Num_echoes, perchè alla fine ricavo un vettore con 6 valori, uno per ogni echo, passando per prodotto pesi prob 
	Class_WP.setZero();
	if (beam_z(ibin)  == Z_missing) {
	  unsigned ID=0;
	  res[ibin]=ID;
	  counter[ID]++;

	  continue;
        }
        
//eseguo un test unico per VRAD e WRAD e assegno tutte le prob assieme
	if (beam_v(ibin) != bin_wind_magic_number  ) {				//	VRAD
	  Pij(0,1)=trap(v_ny*Traps[0][0][0], v_ny*Traps[0][0][1], v_ny*Traps[0][0][2], v_ny*Traps[0][0][3], beam_v(ibin));	 		// METEO
	  double prob_v = trap (v_ny, v_ny, v_ny, v_ny, beam_v(ibin),v_ny);	
	   //cout<<"prob_v computed: "<<prob_v<<endl;
	   //for(int e=1;e<Num_echoes;e++){ Pij(e,1) = prob_v };		
	  Pij(1,1)=trap(Traps[0][1][0],Traps[0][1][1],Traps[0][1][2],Traps[0][1][3], beam_v(ibin));	  		        // CLUTTER		 	
           Pij(2,1)=prob_v;	  					// INTERF. Strong		
           Pij(3,1)=prob_v;	  					// INTERF. Med.		
           Pij(4,1)=prob_v;	  					// INTERF. Weak		
	} else {
	   Pij(0,1)=0.3;	// METEO		
	   Pij(1,1)=1.;		// CLUTTER		  
           Pij(2,1)=1.;		// INTERF. Strong	  
           Pij(3,1)=1.;		// INTERF. Med.		  
           Pij(4,1)=1.;		// INTERF. Weak		  
	}

	// WRAD
	Pij(0,2)= trap(Traps[1][0][0],Traps[1][0][1],Traps[1][0][2],Traps[1][0][3],beam_w(ibin));                   // METEO
	Pij(1,2)= trap(Traps[1][1][0],Traps[1][1][1],Traps[1][1][2],Traps[1][1][3],beam_w(ibin));                  // CLUTTER
	double prob_w = trap(Traps[1][2][0],Traps[1][2][1],Traps[1][2][2],Traps[1][2][3],beam_w(ibin));
	Pij(2,2) = prob_w;                                           // INTERF MULTIPLE
	Pij(3,2) = prob_w;                                           // INTERF. Med.
	Pij(4,2) = prob_w;                                           // INTERF. Weak
	
// METEO		
	Pij(0,0) = trap(Traps[2][0][0],Traps[2][0][1],Traps[2][0][2],Traps[2][0][3],beam_z(ibin), Traps[2][0][4]);				//	Z
	Pij(0,3) = trap (Traps[3][0][0],Traps[3][0][1],Traps[3][0][2],Traps[3][0][3], beam_sd(ibin));		//	SD_2D
	Pij(0,4) = trap (Traps[4][0][0],Traps[4][0][1],Traps[4][0][2],Traps[4][0][3], beam_sdray(ibin));		//	SD_RAY
	Pij(0,5) = trap (Traps[5][0][0],Traps[5][0][1],Traps[5][0][2],Traps[5][0][3], beam_sdaz(ibin));		//	SD_AZ

// CLUTTER		
	Pij(1,0) = trap (Traps[2][1][0],Traps[2][1][1],Traps[2][1][2],Traps[2][1][3], beam_z(ibin), Traps[2][1][4]);		//	Z
	Pij(1,3) = trap (Traps[3][1][0],Traps[3][1][1],Traps[3][1][2],Traps[3][1][3], beam_sd(ibin));		//	SD_2D
	Pij(1,4) = trap (Traps[4][1][0],Traps[4][1][1],Traps[4][1][2],Traps[4][1][3], beam_sdray(ibin));	//	SD_RAY
	Pij(1,5) = trap (Traps[5][1][0],Traps[5][1][1],Traps[5][1][2],Traps[5][1][3], beam_sdaz(ibin));		//	SD_AZ
// INTERF. Strong	
	Pij(2,0) = trap(Traps[2][2][0],Traps[2][2][1],Traps[2][2][2],Traps[2][2][3],beam_z(ibin), Traps[2][2][4]);		//	Z
	Pij(2,3) = trap (Traps[3][2][0],Traps[3][2][1],Traps[3][2][2],Traps[3][2][3], beam_sd(ibin));		//	SD_2D
	Pij(2,4) = trap (Traps[4][2][0],Traps[4][2][1],Traps[4][2][2],Traps[4][2][3], beam_sdray(ibin));		//	SD_RAY
	Pij(2,5) = trap (Traps[5][2][0],Traps[5][2][1],Traps[5][2][2],Traps[5][2][3], beam_sdaz(ibin));		//	SD_AZ
// INTERF. Med.			
	Pij(3,0) = trap(Traps[2][3][0],Traps[2][3][1],Traps[2][3][2],Traps[2][3][3],beam_z(ibin), Traps[2][3][4]);		//	Z
	Pij(3,3) = trap (Traps[3][3][0],Traps[3][3][1],Traps[3][3][2],Traps[3][3][3], beam_sd(ibin));		//	SD_2D
	Pij(3,4) = trap (Traps[4][3][0],Traps[4][3][1],Traps[4][3][2],Traps[4][3][3], beam_sdray(ibin));		//	SD_RAY
	Pij(3,5) = trap (Traps[5][3][0],Traps[5][3][1],Traps[5][3][2],Traps[5][3][3], beam_sdaz(ibin));		//	SD_AZ
// INTERF. Weak		
	Pij(4,0) = trap(Traps[2][4][0],Traps[2][4][1],Traps[2][4][2],Traps[2][4][3],beam_z(ibin), Traps[2][4][4]);		//	Z
	Pij(4,3) = trap (Traps[3][4][0],Traps[3][4][1],Traps[3][4][2],Traps[3][4][3], beam_sd(ibin));		//	SD_2D
	Pij(4,4) = trap (Traps[4][4][0],Traps[4][4][1],Traps[4][4][2],Traps[4][4][3], beam_sdray(ibin));		//	SD_RAY
	Pij(4,5) = trap (Traps[5][4][0],Traps[5][4][1],Traps[5][4][2],Traps[5][4][3], beam_sdaz(ibin));		//	SD_AZ
    
//---- fine calcolo probabilità
// Calcolo classe appartenenza
        
	Class_WP = ((Wij.array()*Pij.array()).matrix()*VectorXd::Ones(Num_entries)).array()/(Wij*VectorXd::Ones(Num_entries)).array();
	unsigned i,ID;
	Class_WP.maxCoeff(&i);
	ID=i;
	if (Class_WP(i) < 0.1 ) ID=5;
	res[ibin]=ID;
	//printf("ID %d \n",ID);
	counter[ID]++;

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



void Cleaner::evaluateCleanID(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v,PolarScan<unsigned char>& scan_cleanID,  unsigned iel)
{
    return evaluateCleanID(scan_z, scan_w, scan_v, scan_cleanID, scan_v.undetect,iel);
  //return evaluateClassID(scan_z, scan_w, scan_v, scan_cleanID, scan_v.undetect, radar, iel);
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

void Cleaner::evaluateClassID(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& scan_zdr, PolarScan<double>& scan_rohv, PolarScan<double>& scan_sqi, PolarScan<double>& scan_snr, PolarScan<double>& scan_zvd, PolarScan<unsigned char>& scan_cleanID, PolarScan<double>& scan_DiffProb, double bin_wind_magic_number,const string radar, unsigned iel, bool force_meteo)
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

    cout<<"VRAD OFFSET="<<scan_v.offset<<" , GAIN="<<scan_v.gain<<endl;    

// compute texture volumes
    radarelab::volume::Scans<double>   Z_S, SD2D, SD_Ray, SD_Az, ZDR_S, ZDR_SD2D;
    Z_S.push_back(scan_z);
    ZDR_S.push_back(scan_zdr);
    radarelab::volume::textureSD( Z_S,SD2D, 1000. , 3,false);
    radarelab::volume::textureSD( Z_S,SD_Ray, scan_z.cell_size*21 , 360./scan_z.beam_count,true);
    radarelab::volume::textureSD( Z_S,SD_Az, scan_z.cell_size , 5*360./scan_z.beam_count,true);

    radarelab::volume::textureSD( ZDR_S,ZDR_SD2D, 1000. , 3,false);

    bool stamp=false;
    for (unsigned i = 0; i <beam_count ; ++i) 
    {
      // bool stamp=false;
      //if(i==50) stamp=true;//311
      //vector<unsigned char> corrected = cleaner.eval_classID_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i), SD2D[0].row(i), scan_zdr.row(i), scan_rohv.row(i), scan_sqi.row(i), scan_snr.row(i), scan_zvd.row(i),  SD_Ray[0].row(i), SD_Az[0].row(i), ZDR_SD2D[0].row(i), i, radar, scan_v.offset, stamp);
      //auto [corrected, diff_prob] = cleaner.eval_classID_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i), SD2D[0].row(i), scan_zdr.row(i), scan_rohv.row(i), scan_sqi.row(i), scan_snr.row(i), scan_zvd.row(i),  SD_Ray[0].row(i), SD_Az[0].row(i), ZDR_SD2D[0].row(i), i, radar, scan_v.offset, stamp, force_meteo);
      vector<unsigned char> corrected = cleaner.eval_classID_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i), SD2D[0].row(i), scan_zdr.row(i), scan_rohv.row(i), scan_sqi.row(i), scan_snr.row(i), scan_zvd.row(i),  SD_Ray[0].row(i), SD_Az[0].row(i), ZDR_SD2D[0].row(i), i, radar, scan_v.offset, stamp, force_meteo);
      for (unsigned ib = 0; ib < beam_size; ++ib){
	   scan_cleanID(i,ib)=corrected[ib];
	   //scan_DiffProb(i,ib)=diff_prob[ib];
      }
    }
    //cout<<"Scan_DiffProb(183,848)="<<scan_DiffProb(183,848)<<endl;
    //cout << typeid(scan_DiffProb(0,100)).name() << endl; 
}

// senza zdr
  void Cleaner::evaluateClassID(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<unsigned char>& scan_cleanID, double bin_wind_magic_number, const string radar, unsigned iel)
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


    for (unsigned i = 0; i <beam_count ; ++i) 
    {
      vector<unsigned char> corrected = cleaner.eval_classID_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i), SD2D[0].row(i), SD_Ray[0].row(i), SD_Az[0].row(i), i, radar, scan_v.offset);
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
    if (set_undetect) new_val=scan_z.nodata;

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
    return clean(scan_z, scan_w, scan_v, scan_zdr, scan_v.undetect,iel,set_undetect);
}
void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, PolarScan<double>& scan_zdr, double bin_wind_magic_number, unsigned iel, bool set_undetect )
{
	cout<<"Chiamato cleaner "<<set_undetect<<endl;
    if (scan_z.beam_count != scan_w.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_w beam_count");
    if (scan_z.beam_size != scan_w.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_w beam_size");

    if (scan_z.beam_count != scan_v.beam_count)
        throw std::runtime_error("scan_z beam_count no equal to scan_v beam_count");
    if (scan_z.beam_size != scan_v.beam_size)
        throw std::runtime_error("scan_z beam_size no equal to scan_v beam_size");
    double z_val=scan_z.nodata;
    if(set_undetect) z_val=scan_z.undetect;
//Cleaner cleaner(scan_z.undetect, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);
    Cleaner cleaner(z_val, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);

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
//printf("valore new_val ---> %f\n",new_val);
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

double Cleaner::trap(double x1, double x2, double x3, double x4, double val, double x5) const
{
        if((val<=x3&&val>=x2)) return 1.;
	else if(val<x2&&val>x1) return val/(x2-x1)-x1/(x2-x1);
	else if (val<x4&&val>x3) return val/(x3-x4)-x4/(x3-x4);
        else if(val<=x5) return 1.;
	else return 0.; // (val<=x1||val>=x4)

}

vector<string> Cleaner::read_matrix_from_txt(string fin) const
{

  /*
  legge file txt e mette il contenuto in un vettore
   */
  vector<string> myVector;
  ifstream f(fin, ifstream::in);
  string line;
  
  if(f.is_open()){
    while(getline(f,line)){
      stringstream stream (line);
      while( getline(stream, line, ' ')){
	//inserisco controllo commenti che li gestisce se attaccati tipo //commento
	if(line[0]=='\\'){ continue;}
	else{
	  myVector.push_back(line);
	}
      }
    }
  }
  return myVector;
}


}
}
