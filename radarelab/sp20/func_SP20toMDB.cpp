/*=================================================================
 omstart func_SP20toMDB

  Questo file contiene le funzioni che leggono l'header del DBP
  nel nuovo formato SP20.Inoltre contiene la funzione che converte 
  l'header.Queste funzioni sono chiamate dal programma cum_bac_SP20_[]
 
 omend 
  ================================================================== */
#include <string.h>
#include <math.h>
#include "func_SP20read.h"

/*=============================================================
   HEADER DI VOLUME:  acquisizione e decodifica dei dati
    relativi al DBP e conversione nel vecchio formato.
  =============================================================*/

/* ==================== */
char get_DBP_scan_mode(unsigned char hd[8][4])
/* ==================== */
{
  return (hd [0][1]>>1 & 0x03);
}


/* ==================== */
float get_DBP_rot_vel(unsigned char hd[8][4])
/* ==================== */
{
  //  printf("%x %x\n", hd[0][1], hd[2][2]);

  return ( (  ((hd[0][1]>>3)&0x1f) + ((hd[2][2]>>6)&0x03)*0.1  ));

}


/* ==================== */
char get_DBP_pulse_lenght(unsigned char hd[8][4])
/* ==================== */
{
  return (hd [0][3]>>6 & 0x03);
}


/* ==================== */
float get_DBP_el_min(unsigned char hd[8][4])
/* ==================== */
{
  //  printf("el_min: %x %x \n", hd[1][1], hd[1][2]);

  return (  ( hd[1][1] | ((hd[1][2] &0x0f) <<8) )  *(360./4096.) );
}


/* ==================== */
float get_DBP_el_max(unsigned char hd[8][4])
/* ==================== */
{
  //  printf("el_max: %x %x \n", hd[2][1], hd[2][2]);

  return (  ( hd[2][1] | ((hd[2][2] &0x0f) <<8) )   *(360./4096.) );
}


/* ==================== */
float get_DBP_az_min(unsigned char hd[8][4])
/* ==================== */
{
  //  printf("az_min: %x %x \n", hd[4][1], hd[4][2]);

  return ( ( hd[4][1] | ((hd[4][2] &0x0f) <<8) ) *(360./4096.) );
}


/* ==================== */
float get_DBP_az_max(unsigned char hd[8][4])
/* ==================== */
{
  //  printf("az_max: %x %x \n", hd[5][1], hd[5][2]);

  return (  (hd[5][1] | ( (hd[5][2] &0x0f) <<8)) *(360./4096.) );
}


/* ==================== */
short get_DBP_pulse_number(unsigned char hd[8][4])
/* ==================== */
{
  return ( (short)( ((hd[4][2]>>4) & 0x0f) | ((hd[4][3] &0x1f) <<4) ));
}


/* ============================================ */
void get_DBP_flag_quantities (unsigned char hd[8][4], char *flag_quantities)
/* ============================================ */
{
  //misura di Z
  flag_quantities[0] = hd[1][3] & 0x01;
  //misura di Zdr
  flag_quantities[1] = hd[1][3]>>1 & 0x01;
  //misura di V
  flag_quantities[2] = hd[1][3]>>2 & 0x01;
  //misura di sV
  flag_quantities[3] = hd[1][3]>>3 & 0x01;

  return;

}



/* ==================== */
char get_DBP_cell_size (unsigned char hd[8][4])
/* ==================== */
{
  return ( (hd[1][3]>>4) & 0x07);
}
	   


/* ==================== */
char get_DBP_Dual_PRF (unsigned char hd[8][4])
/* ==================== */
{
  return ( (hd[1][3]>>7) & 0x01);
} 


//da verificare
/* ========================= */
char get_DBP_tipo_filtro(unsigned char hd[8][4])
/* ========================= */
{
  return ( (char) ((hd[2][3]>>4)&0x0f) );
}

/* =========================== */
char get_DBP_corr_clutter (unsigned char hd[8][4])
/* =========================== */
{
  return ( (char) ((hd[2][3]>>2)&0x01) );
}

/* ============================= */
char get_DBP_filtro_clutter (unsigned char hd[8][4])
/* ============================= */
{
  return ( (char) ((hd[2][3]>>3)&0x01) );
}

/* ================================== */
char get_DBP_anticlutter_mappato (unsigned char hd[8][4])
/* ================================== */
{
  return ( (char) ((hd[2][3])>>1 &0x01) );
}

//da verificare
/* ==================== */
char get_DBP_Z_range (unsigned char hd[8][4])
/* ==================== */
{
  return ( (char) (hd[3][2]>>6)&0x03 );
}


/* ================================= */
char get_DBP_calibration_factor (unsigned char hd[8][4])
/* ================================= */
{
  return( (char) (hd[3][3]&0x0f) + (((hd[3][3]>>4)&0x0f)*10) );
}



/* ========================= */
char get_DBP_stalo_code (unsigned char hd[8][4])
/* ========================= */
{
  //  printf("stalo_code esadecimale: %x\n",hd[7][3]);

  return ( hd[7][3] & 0x07 );
}





//per richiamare le singole routine
/*===========================================*/
void decode_header_DBP_SP20 (HD_DBP_SP20_RAW *hd_raw, HD_DBP_SP20_DECOD *hd_decod)
/*===========================================*/
{
  int i;
  
 
  hd_decod->scan_mode = get_DBP_scan_mode(hd_raw->frame);
  hd_decod->rot_vel = get_DBP_rot_vel(hd_raw->frame);
  hd_decod->pulse_lenght = get_DBP_pulse_lenght(hd_raw->frame);
  hd_decod->el_min = get_DBP_el_min(hd_raw->frame);
  hd_decod->el_max = get_DBP_el_max(hd_raw->frame);
  hd_decod->az_min = get_DBP_az_min(hd_raw->frame);
  hd_decod->az_max = get_DBP_az_max(hd_raw->frame);

  get_DBP_flag_quantities (hd_raw->frame, hd_decod->flag_quantities);

  hd_decod->cell_size = get_DBP_cell_size(hd_raw->frame);
  hd_decod->Dual_PRF = get_DBP_Dual_PRF(hd_raw->frame);
  hd_decod->filtro_clutter = get_DBP_filtro_clutter(hd_raw->frame);
  hd_decod->tipo_filtro = get_DBP_tipo_filtro(hd_raw->frame);
  hd_decod->anticlutter_mappato = get_DBP_anticlutter_mappato(hd_raw->frame);
  hd_decod->corr_clutter = get_DBP_corr_clutter(hd_raw->frame);
  hd_decod->Z_range = get_DBP_Z_range(hd_raw->frame);
  hd_decod->calibration_factor = get_DBP_calibration_factor(hd_raw->frame);
  hd_decod->pulse_number = get_DBP_pulse_number(hd_raw->frame);
  hd_decod->stalo_code = get_DBP_stalo_code(hd_raw->frame);
	if (hd_raw->num_ele >= 41 )
	      hd_decod->num_ele=hd_raw->num_ele-41; 
	else
	      hd_decod->num_ele=hd_raw->num_ele; 

  hd_decod->corr_pot=hd_raw->corr_pot;
 
  for(i=0; i<hd_decod->num_ele; i++)             //ciclo sul numero elementi significativi
    hd_decod->ele[i]=hd_raw->ele[i]*0.0879; 

  return ; 
}




/*converte nelle variabili del vecchio formato*/
/*==============================================*/
void convert_format(HD_DBP_SP20_DECOD *hd_decod, BEAM_HD_SP20_INFO *beam_info, T_MDB_data_header *old_header)   //raggio= beam_info
/*==============================================*/
{
  int i, k;

  old_header->norm.maq.acq_date=beam_info->time;
  old_header->norm.maq.scans_type=hd_decod->scan_mode;
  old_header->norm.maq.rot_vel=hd_decod->rot_vel;
  old_header->norm.maq.imp_duration=hd_decod->pulse_lenght;

  old_header->norm.maq.grand=0;
  for(k=0; k<4; k++)
  {
    old_header->norm.maq.grand=old_header->norm.maq.grand +
      (int)( pow((double)2,(double)k) * (int)hd_decod->flag_quantities[k]);
    // printf("\nflag_quantities[k]=%d  grand=%d\n  ",hd_decod->flag_quantities[k], old_header->norm.maq.grand);
  }
  old_header->norm.maq.resolution=hd_decod->cell_size;
  old_header->norm.maq.vel_range=hd_decod->Dual_PRF;
  old_header->norm.maq.declutter_rsp=hd_decod->filtro_clutter;

  //  printf("%1d 1d\n", old_header->norm.maq.declutter_rsp,hd_decod->filtro_clutter);

  old_header->norm.maq.filter_value=hd_decod->tipo_filtro;
  old_header->norm.maq.type_declutter=hd_decod->anticlutter_mappato;
  old_header->norm.maq.corr_Z=hd_decod->corr_clutter;
  old_header->norm.maq.spare[1]=hd_decod->Z_range;
  if(hd_decod->Z_range ==0) 
    old_header->norm.maq.spare[0]=2;
  if(hd_decod->Z_range ==1)
    old_header->norm.maq.spare[0]=1;
  if(hd_decod->Z_range ==2) 
    old_header->norm.maq.spare[0]=3;

  old_header->norm.maq.num_imp=hd_decod->pulse_number;
  old_header->norm.maq.num_el=hd_decod->num_ele; 
  //ciclo sul num di elementi significativi
  for(i=0; i<old_header->norm.maq.num_el ;i++) 
    old_header->norm.maq.value[i]=hd_decod->ele[i]/FATT_MOLT_EL; 
  return ; 
}


/* ========================================= */    
void PrintHeaderDBP(HD_DBP_SP20_DECOD *hd_decod)
/* ========================================= */
{
  int i;

  printf("HEADER DBP:\n");
  printf("%2d ",hd_decod->scan_mode);
  printf("%2d ",hd_decod->flag_quantities[0]);
  printf("%1d ",hd_decod->flag_quantities[1]);
  printf("%1d ",hd_decod->flag_quantities[2]);
  printf("%1d ",hd_decod->flag_quantities[3]);
  printf("%2d ",hd_decod->Z_range);
  printf("%2d ",hd_decod->Dual_PRF); 
  printf("%2d ",hd_decod->pulse_lenght);
  printf("%2d ",hd_decod->pulse_number);
  printf("%2d\n ",hd_decod->cell_size);
  printf("%9.5f ",hd_decod->az_min);
  printf("%9.5f ",hd_decod->az_max);
  printf("%8.5f ",hd_decod->el_min);
  printf("%8.5f\n ",hd_decod->el_max);
  printf("%5.2f ",hd_decod->rot_vel);               
  printf("%2d ",hd_decod->filtro_clutter);
  printf("%2d ",hd_decod->tipo_filtro);
  printf("%2d ",hd_decod->anticlutter_mappato);
  printf("%2d ",hd_decod->corr_clutter);
  printf("%2d ",hd_decod->calibration_factor);
  printf("%3d ",hd_decod->stalo_code);
  printf("%3d ",hd_decod->num_ele); 
  printf("%3d\n ",hd_decod->corr_pot);
  // ciclo sul numero di elementi significativi
  for(i=0;i<hd_decod->num_ele;i++)
    printf("%9.5f ",hd_decod->ele[i]); 

  printf("\n");
  return;
}



/*===============================================*/
void PrintOldHeader (T_MDB_data_header *old_header)
/*===============================================*/
{
  int i;

  printf("OLD HEADER :\n");
  printf("%2d ", old_header->norm.maq.scans_type);
  printf("%2d ", old_header->norm.maq.grand );
  printf("%2d ", old_header->norm.maq.spare[0]);
  printf("%2d ", old_header->norm.maq.vel_range);
  printf("%2d ", old_header->norm.maq.imp_duration);
  printf("%2d ", old_header->norm.maq.num_imp);
  printf("%2d\n ", old_header->norm.maq.resolution);
  printf("%5.2f ", (double)old_header->norm.maq.rot_vel);
  printf("%2d ", old_header->norm.maq.declutter_rsp);
  printf("%2d ", old_header->norm.maq.filter_value);
  printf("%2d ", old_header->norm.maq.type_declutter);
  printf("%2d ", old_header->norm.maq.corr_Z);
  printf("%2d ", old_header->norm.maq.num_el); 

  printf("\n"); 

  for(i=0; i<old_header->norm.maq.num_el; i++) 
  {
    printf("%9.5f ", old_header->norm.maq.value[i] * FATT_MOLT_EL);
  }
  printf("\n");  
  
}    

/*=============================================================
                          FINE   HEADER DI VOLUME
  =============================================================*/



/*=============================================================
                          FORMATO BEAM
  =============================================================*/

/* ==================== */
short check_flag(unsigned char hd[])
/* ==================== */
{
  if (hd[0]|hd[1]|hd[2])
    return(-1);
  return(0);
}

/* ======================== */
short BCD_decoding(char byte)
/* ======================== */
{
  unsigned char units,tenths;
  units=byte&0x0f;
  tenths=(byte>>4)&0x0f;
  return ((short)(tenths*10+units));
}

/* ====================== */
time_t get_date(unsigned char hd[], struct tm *tm_date)
/* ====================== */
{
  time_t differenza = 10000;
  time_t UTC_time;
  //  printf(" %2.2x %2.2x %2.2x %2.2x %2.2x %2.2x \n", hd[24], hd[23], hd[22], hd[20], hd[19], hd[18]);

  tm_date->tm_year=BCD_decoding(hd[24]);           /* Year     */
  if (tm_date->tm_year < 50) 
    tm_date->tm_year=tm_date->tm_year+100;
  
  tm_date->tm_mon=BCD_decoding(hd[23])-1;          /* Month    */
  tm_date->tm_mday=BCD_decoding(hd[22]);            /* Day      */

  /*------
    NON estraggo il dato ma lo ricalcolo con la funzione mktime
    date->tm_wday=BCD_decoding(hd->week_day);
    --------*/
  tm_date->tm_hour=BCD_decoding(hd[20]);           /* Hour     */
  tm_date->tm_min=BCD_decoding(hd[19]);       /* Minute   */
  tm_date->tm_sec=BCD_decoding(hd[18]);         /* Second   */
  
  /*------------
    A questo punto riempio tutto la struttuta tm_date
    ------------*/
  tm_date->tm_isdst=0;
  tm_date->tm_gmtoff=0;
  
  /*------------
    chiamo la mktime
    ------------*/
  /*  printf("--- %s\n",asctime(tm_date));*/


  /*-----------------------------------------
    calcolo la differenza tra l'ora locala della macchina su cui gira l'applicativo
    e il fuso GMT perche' mktime() restituisce i secondi lavorando sul fuso orario locale
    -----------------------------------------------*/
  differenza = differenza-mktime(gmtime(&differenza));
  UTC_time= mktime(tm_date)-differenza;
  tm_date=gmtime(&UTC_time);
  
  return (UTC_time);
}




/* ====================== */
void get_asc_date(struct tm *tm_date, char *asc_date)
/* ====================== */
{
  const char *mese[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
  //char *weekday[7]={"Mon","Tue","Wed","Thu","Fri","Sat","Sun"};
  const char *weekday[7]={"Sun","Mon","Tue","Wed","Thu","Fri","Sat"};
  
  sprintf(asc_date,"%s %s %2d %2.2d:%2.2d:%2.2d %4d",weekday[tm_date->tm_wday],
	  mese[tm_date->tm_mon],tm_date->tm_mday,tm_date->tm_hour,
	  tm_date->tm_min,tm_date->tm_sec,tm_date->tm_year+1900);
 /* sprintf(asc_date,"%s %s %2d %2.2d:%2.2d:%2.2d %4d",weekday[tm_date->tm_wday-1],
	  mese[tm_date->tm_mon],tm_date->tm_mday,tm_date->tm_hour,
	  tm_date->tm_min,tm_date->tm_sec,tm_date->tm_year+1900);
*/
  return;
}


/* ============================ */
char check_data_validity(unsigned char hd[])
/* ============================ */
{
  if ((hd[12]>>1)&0x01)
    return('V');
  return('N');
}

/* ======================== */
short get_real_power(unsigned char hd[])
/* ======================== */
{
  return((short)(hd[15]*4));
}

/* ======================== */
short get_eval_power(unsigned char hd[])
/* ======================== */
{
  return((short)(hd[16]*4));
}

/* ====================== */
char get_scan_mode(unsigned char hd[])
/* ====================== */
{
  switch ((hd[8]>>2)&0x03)
    {
    case 0:
      return('V');
      
    case 1:
      return('V');
      
    case 2:
      return('A');
      
    case 3:
      return('E');
      
    default:
      return('?');
    }
}


/* ======================= */
float get_elevation(unsigned char hd[])
/* ======================= */
{
  return((float)((hd[5] |(hd[6] &0x1f)<<8 )*ACF));
//  return((float)((hd[5] |(hd[6] &0x1f)<<8 )*ACF)/8.);

}

/* ===================== */
float get_azimuth(unsigned char hd[])
/* ===================== */
{
  return((float)((hd[3] |(hd[4] &0x1f)<<8 )*ACF));
//  return((float)((hd[3] |(hd[4] &0x1f)<<8 )*ACF)/8.);
}

/* ================================== */
void get_quantities(unsigned char hd[], char quantities[5])
/* ================================== */
{
  strcpy(quantities,"");
  if (hd[9]&0x01)
    strcat(quantities,"Z");
  if ((hd[9]>>1)&0x01)
    strcat(quantities,"D");
  if ((hd[9]>>2)&0x01)
    strcat(quantities,"V");
  if ((hd[9]>>7)&0x01)
    strcat(quantities,"S");
  return;
}

/* ================================== */
void get_flag_quantities(unsigned char hd[], char flag_quantities[])
/* ================================== */
{
  /* Z */  
  flag_quantities[0]=hd[9]&0x01;

  /* ZDR */
  flag_quantities[1]=(hd[9]>>1)&0x01;
  
  /* V */
  flag_quantities[2]=(hd[9]>>2)&0x01;
  
  /* sigma V */
  flag_quantities[3]=(hd[9]>>7)&0x01;

  return;
}

/* ===================== */
int check_dualPRF(unsigned char hd[])
/* ===================== */
{
  return ((hd[9]>>3)&0x01);
}

/* ==================== */
char get_Z_range(unsigned char hd[])
/* ==================== */
{
  switch (hd[14]&0x03)
    {
    case 0:
      return('M');
      
    case 1:
      return('H');
      
    case 2:
      return('L');
      
    default:
      return('?');
    }
}

/* ====================== */
char get_pulse_len(unsigned char hd[])
/* ====================== */
{
  switch (hd[8]&0x03)
    {
    case 0:
      return('L');
      
    case 1:
      return('M');
      
    case 2:
      return('S');
      
    default:
      return('?');
    }
}

/* ========================== */
   short get_pulse_num(unsigned char hd[])
/* ========================== */
{
  return((short)(hd[13] |(hd[14]>>6 &0x03)<<8 ));
}

/* ====================== */
   char get_cell_size(unsigned char hd[])
/* ====================== */
{
  return((hd[9]>>4)&0x07);
}

/* ====================== */
   short get_cell_num(unsigned char hd[])
/* ====================== */
{
  return((short)(hd[10] |(hd[11]>>6 &0x03)<<8 ));
}

/* ====================== */
   char get_channel(unsigned char hd[])
/* ====================== */
{
  return((hd[11]>>3)&0x07);
}

/* =========================================== */
   void get_controller_identity(unsigned char hd[], char controller[4])
/* =========================================== */
{
  switch ((hd[8]>>5)&0x03)
    {
    case 0:
      strcpy(controller,"LOC");
      break;
    case 1:
      strcpy(controller,"RM1");
      break;
    case 2:
      strcpy(controller,"RM2");
      break;
    default:
      strcpy(controller,"ERR");
    }
  return;
}

/* =========================================== */
char get_filtro_clutter(unsigned char *hd)
/* =========================================== */
{
  return((hd[12]>>7)&0x01);
}

/* =========================================== */
char get_tipo_filtro(unsigned char *hd)
/* =========================================== */
{
  return((hd[12]>>2)&0x0f);
}

/* =========================================== */
char get_anticlutter_mappato(unsigned char *hd)
/* =========================================== */
{
  return(hd[12]&0x01);
}

/* =========================================== */
char get_corr_clutter(unsigned char *hd)
/* =========================================== */
{
  return((hd[12]>>6)&0x01);
}


/* ================================= */
BEAM_HD_SP20_INFO * decode_header_sp20(unsigned char* hd_char, BEAM_HD_SP20_INFO *beam_info)   //info=beam_info
/* ================================= */
{
  char quantities[5];
  char controller[4];

  beam_info->time=get_date(hd_char,&(beam_info->tm_date));
  /*-----------------------------
   Nota - PPA - 29/09/2009
non sappiamo perchè ma la get_date_from_name mi restituisce la struttura tm_date errata.
La ricalcoliamo dal tempo in secondi
----------------------------------------*/
  memcpy(&(beam_info->tm_date),gmtime(&beam_info->time),sizeof(beam_info->tm_date));
  get_asc_date(&(beam_info->tm_date),&(beam_info->date[0]));

  beam_info->valid_data=check_data_validity(hd_char);

  beam_info->real_power=get_real_power(hd_char);

  beam_info->eval_power=get_eval_power(hd_char);

  beam_info->scan_mode=get_scan_mode(hd_char);


  beam_info->azimuth=get_azimuth(hd_char);

  beam_info->elevation=get_elevation(hd_char);

  get_quantities(hd_char,quantities);	
  
  strcpy(beam_info->quantities,quantities);
  
  get_flag_quantities(hd_char, beam_info->flag_quantities);
  
  
  beam_info->Z_range=get_Z_range(hd_char);
  
  beam_info->pulse_lenght=get_pulse_len(hd_char);
  
  beam_info->pulse_number=get_pulse_num(hd_char);
  
  beam_info->cell_size=get_cell_size(hd_char);
  
  beam_info->cell_num=get_cell_num(hd_char);
  
  beam_info->channel=get_channel(hd_char);	
  
  
  get_controller_identity(hd_char,controller);
  
  strcpy(beam_info->controller,controller);
  
  if (check_dualPRF(hd_char))
    beam_info->PRF='D';
  else
    beam_info->PRF='S';
  

  beam_info->filtro_clutter=get_filtro_clutter(hd_char);
  beam_info->tipo_filtro=get_tipo_filtro(hd_char);
  beam_info->anticlutter_mappato=get_anticlutter_mappato(hd_char);
  
  beam_info->corr_clutter=get_corr_clutter(hd_char);
  

  return(beam_info);
}


/* ================================= */
BEAM_HD_SP20_INFO * decode_header_sp20_date_real_eo_name(unsigned char *hd_char, BEAM_HD_SP20_INFO *beam_info, char *nome_file)
/* ================================= */
{
  
  char quantities[5];
  char controller[4];


  beam_info->time=get_date(hd_char,&(beam_info->tm_date));
  if (beam_info->time == -1) 
    beam_info->time=get_date_from_name((T_MDB_data_header*)hd_char,&(beam_info->tm_date),nome_file);
  /*-----------------------------
   Nota - PPA - 29/09/2009
non sappiamo perchè ma la get_date_from_name mi restituisce la struttura tm_date errata.
La ricalcoliamo dal tempo in secondi
----------------------------------------*/
  memcpy(&(beam_info->tm_date),gmtime(&beam_info->time),sizeof(beam_info->tm_date));

  get_asc_date(&(beam_info->tm_date),&(beam_info->date[0]));
  beam_info->valid_data=check_data_validity(hd_char);
  beam_info->real_power=get_real_power(hd_char);
  beam_info->eval_power=get_eval_power(hd_char);
  beam_info->scan_mode=get_scan_mode(hd_char);
  
  beam_info->azimuth=get_azimuth(hd_char);
  beam_info->elevation=get_elevation(hd_char);
  
  get_quantities(hd_char,quantities);	
  strcpy(beam_info->quantities,quantities);
  get_flag_quantities(hd_char, beam_info->flag_quantities);
  
  beam_info->Z_range=get_Z_range(hd_char);
  beam_info->pulse_lenght=get_pulse_len(hd_char);
  beam_info->pulse_number=get_pulse_num(hd_char);
  beam_info->cell_size=get_cell_size(hd_char);
  beam_info->cell_num=get_cell_num(hd_char);
  beam_info->channel=get_channel(hd_char);	
  
  get_controller_identity(hd_char,controller);
  strcpy(beam_info->controller,controller);
  if (check_dualPRF(hd_char))
    beam_info->PRF='D';
  else
    beam_info->PRF='S';
  
  beam_info->filtro_clutter=get_filtro_clutter(hd_char);
  beam_info->tipo_filtro=get_tipo_filtro(hd_char);
  beam_info->anticlutter_mappato=get_anticlutter_mappato(hd_char);
  beam_info->corr_clutter=get_corr_clutter(hd_char);
  
  return(beam_info);
}

/* ================================= */
BEAM_HD_SP20_INFO * decode_header_sp20_date_from_name(unsigned char *hd_char, BEAM_HD_SP20_INFO *beam_info, char *nome_file)
/* ================================= */
{
  
  char quantities[5];
  char controller[4];

  /*
  printf(" prima get_date_from_name --- ") ;
  printf("%2d %2d %4d %2d %2d\n",
	 beam_info->tm_date.tm_mday,beam_info->tm_date.tm_mon,beam_info->tm_date.tm_year,
	 beam_info->tm_date.tm_hour,beam_info->tm_date.tm_min);
  */

  beam_info->time=get_date_from_name((T_MDB_data_header*)hd_char,&(beam_info->tm_date),nome_file);
  
/*   printf(" dopo get_date_from_name --- ") ;
  printf("%2d %2d %4d %2d %2d\n",
	 beam_info->tm_date.tm_mday,beam_info->tm_date.tm_mon,beam_info->tm_date.tm_year,
	 beam_info->tm_date.tm_hour,beam_info->tm_date.tm_min);
*/
  /*-----------------------------
   Nota - PPA - 29/09/2009
non sappiamo perchè ma la get_date_from_name mi restituisce la struttura tm_date errata.
La ricalcoliamo dal tempo in secondi
----------------------------------------*/
  memcpy(&(beam_info->tm_date),gmtime(&beam_info->time),sizeof(beam_info->tm_date));
 /*  printf(" dopo get_date_from_name --- ") ;
  printf("%2d %2d %4d %2d %2d  -  %2d\n",
	 beam_info->tm_date.tm_mday,beam_info->tm_date.tm_mon,beam_info->tm_date.tm_year,
	 beam_info->tm_date.tm_hour,beam_info->tm_date.tm_min, beam_info->tm_date.tm_wday);
*/
  get_asc_date(&(beam_info->tm_date),&(beam_info->date[0]));
//printf(" data in chiaro -- %s\n", beam_info->date);
  beam_info->valid_data=check_data_validity(hd_char);
  beam_info->real_power=get_real_power(hd_char);
  beam_info->eval_power=get_eval_power(hd_char);
  beam_info->scan_mode=get_scan_mode(hd_char);
  
  beam_info->azimuth=get_azimuth(hd_char);
  beam_info->elevation=get_elevation(hd_char);
  
  get_quantities(hd_char,quantities);	
  strcpy(beam_info->quantities,quantities);
  get_flag_quantities(hd_char, beam_info->flag_quantities);
  
  beam_info->Z_range=get_Z_range(hd_char);
  beam_info->pulse_lenght=get_pulse_len(hd_char);
  beam_info->pulse_number=get_pulse_num(hd_char);
  beam_info->cell_size=get_cell_size(hd_char);
  beam_info->cell_num=get_cell_num(hd_char);
  beam_info->channel=get_channel(hd_char);	
  
  get_controller_identity(hd_char,controller);
  strcpy(beam_info->controller,controller);
  if (check_dualPRF(hd_char))
    beam_info->PRF='D';
  else
    beam_info->PRF='S';
  
  beam_info->filtro_clutter=get_filtro_clutter(hd_char);
  beam_info->tipo_filtro=get_tipo_filtro(hd_char);
  beam_info->anticlutter_mappato=get_anticlutter_mappato(hd_char);
  beam_info->corr_clutter=get_corr_clutter(hd_char);
  
  return(beam_info);
}

/* ========================================= */    
void PrintHeader(BEAM_HD_SP20_INFO *beam_info)
/* ========================================= */
{
  printf("HEADER BEAM:\n");
  printf("%s ",beam_info->date);
  printf("%li ",beam_info->time);
  printf("%c ",beam_info->scan_mode);
  printf("%s ",beam_info->quantities);
  printf("%c ",beam_info->Z_range);
  printf("%c ",beam_info->PRF); 
  printf("%c ",beam_info->pulse_lenght);
  printf("%i ",beam_info->pulse_number);
  printf("%i ",beam_info->cell_size);
  printf("%i ",beam_info->cell_num);
  printf("%c ",beam_info->valid_data);
  printf("%i ",beam_info->real_power);
  printf("%i ",beam_info->eval_power);
  printf("%9.5f ",beam_info->azimuth);
  printf("%8.5f ",beam_info->elevation);
  printf("%1d ",beam_info->channel);
  printf("%2d ",beam_info->filtro_clutter);
  printf("%2d ",beam_info->tipo_filtro);
  printf("%2d ",beam_info->anticlutter_mappato);
  printf("%2d ",beam_info->corr_clutter);

  printf("\n");
  return;
}


/*=============================================================
                      FINE FORMATO BEAM
  =============================================================*/





