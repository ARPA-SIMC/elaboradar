/*----------------------------------------------------------------------------*/
/*	INCLUDE file						       	      */
/*----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "func_SP20read.h" /* file contenente i prototipi 
		              delle funzioni necessarie per leggere e convertire 
		              i dati dal formato nuovo SP20 al vecchio MDB  */

/* omstart ReadHeaderSP20toMDB
   ========================================
   |  FUNZIONE ReadHeaderSP20toMDB        |  
   | legge i dati nel nuovo formato,      | 
   | header DBP e raggio, e li rende      |
   | compatibili con il vecchio.          |
   ======================================== 
omend */

extern int elev_array[NEL];
extern int nbeam_elev[NEL];
unsigned char header_beam_spc_char[40];

time_t get_date_from_name(T_MDB_data_header* old_data_header, struct tm *tm_date_ext, const char* nome_file)
{
  const char *DBP={"DBP2"};
  struct tm tm_date;
  const char *pointer;
  char s_date[100];
  int year,month,day,hour,min,sec;

  time_t differenza = 10000;
  time_t UTC_time;

  //    printf("sono qui\n");
  pointer=strstr(nome_file,DBP);
  strncpy(s_date,pointer+5,12);
  s_date[12]='\0';
  sscanf(s_date,"%2d%2d%4d%2d%2d",
	 &tm_date.tm_mday,&tm_date.tm_mon,&tm_date.tm_year,
	 &tm_date.tm_hour,&tm_date.tm_min);
  //   printf("data %s\n",s_date);
  
  tm_date.tm_year=tm_date.tm_year-1900;             /* Year     */
  tm_date.tm_mon=tm_date.tm_mon-1;                  /* Month    */
  tm_date.tm_sec=0;                                 /* Second   */
  
  /*------------
    A questo punto riempio tutto la struttuta tm_date
    ------------*/
  tm_date.tm_isdst=0;
  tm_date.tm_gmtoff=0;
  
  /*------------
    chiamo la mktime
    ------------*/
  //  printf("--- %s\n",asctime(&tm_date));
  //  return(1);


  /*-------
    copio la struttura tm_date nella struttura di input/output
    --------*/
  /*  printf("%2d %2d %4d %2d %2d\n",
	 tm_date.tm_mday,tm_date.tm_mon,tm_date.tm_year,
	 tm_date.tm_hour,tm_date.tm_min);
  */
   //memcpy(tm_date_ext,&tm_date,sizeof(tm_date));

  /*
  printf("%2d %2d %4d %2d %2d\n",
	 tm_date_ext->tm_mday,tm_date_ext->tm_mon,tm_date_ext->tm_year,
	 tm_date_ext->tm_hour,tm_date_ext->tm_min);
 
  */
  /*-----------------------------------------
    calcolo la differenza tra l'ora locala della macchina su cui gira l'applicativo
    e il fuso GMT perche' mktime() restituisce i secondi lavorando sul fuso orario locale
    -----------------------------------------------*/
   differenza = differenza-mktime(gmtime(&differenza));
  UTC_time=(mktime(&tm_date)+differenza);
  tm_date= *(gmtime(&UTC_time));
  memcpy(tm_date_ext,&tm_date,sizeof(tm_date));
  /*
  printf ("ancora -");
 printf("%2d %2d %4d %2d %2d  --  %2d\n",
	 tm_date_ext->tm_mday,tm_date_ext->tm_mon,tm_date_ext->tm_year,
	 tm_date_ext->tm_hour,tm_date_ext->tm_min , tm_date_ext->tm_wday);
   return (UTC_time);
  */
   return (UTC_time);
   // return (mktime(tm_date_ext)+differenza);

}


/*===============================================*/
int ReadHeaderSP20toMDB(T_MDB_data_header *old_header, FILE *fp)
/*===============================================*/
{

  HD_DBP_SP20_RAW hd_raw;                                     
  HD_DBP_SP20_DECOD hd_decod;                                 
  BEAM_HD_SP20_INFO beam_info;                                 
  int offset;
  
  if (fread(&hd_raw, sizeof(hd_raw), 1, fp) != 1)
    return NO_OK;
  decode_header_DBP_SP20 (&hd_raw, &hd_decod);
  //printf("------------\n read header \n");
  // PrintHeaderDBP(&hd_decod);  

  while(1){
    if (fread(header_beam_spc_char, sizeof(header_beam_spc_char), 1, fp) != 1)
      return NO_OK;
    decode_header_sp20(header_beam_spc_char, &beam_info);
    //PrintHeader(&beam_info);
    
    offset=strlen(beam_info.quantities)*beam_info.cell_num;
    fseek(fp,offset,SEEK_CUR);

    if (beam_info.valid_data == 'V')
    {    
      convert_format(&hd_decod, &beam_info, old_header);      
      //      PrintOldHeader (old_header);
      //  printf ("%d  -- %d \n",offset,sizeof(header_beam_spc_char));
      fseek(fp, -(offset+sizeof(header_beam_spc_char)),SEEK_CUR);
      break;
    }

  }
  return OK;

}

/*===============================================*/
int ReadHeaderSP20(HD_DBP_SP20_DECOD *new_header_deco, BEAM_HD_SP20_INFO *beam_info, FILE *fp)
/*===============================================*/
{
  HD_DBP_SP20_RAW hd_raw;                                     
  int offset;

  if (fread(&hd_raw, sizeof(hd_raw), 1, fp) != 1)
    return NO_OK;
  decode_header_DBP_SP20 (&hd_raw, new_header_deco);
  //  PrintHeaderDBP(&hd_decod);  

  while(1){
    if (fread(header_beam_spc_char, sizeof(header_beam_spc_char), 1, fp) != 1)
      return NO_OK;
  decode_header_sp20(header_beam_spc_char, beam_info);
  //  PrintHeader(&beam_info);

    offset=strlen(beam_info->quantities)*beam_info->cell_num;
    fseek(fp,offset,SEEK_CUR);

    if (beam_info->valid_data == 'V')
    {    
      //  PrintOldHeader (old_header);
      //  printf ("%d  -- %d \n",offset,sizeof(header_beam_spc_char));
      fseek(fp, -(offset+sizeof(header_beam_spc_char)),SEEK_CUR);
      break;
    }
  }
  return OK;

}

/*===============================================*/
int read_dbp_SP20(
        char *nome_file,
        struct VOL_POL vol_pol_locale[][NUM_AZ_X_PPI],
        T_MDB_data_header   *old_data_header,
        int tipo_dati,
        int nbeam_elev[])
/*===============================================*/
{
  BEAM_HD_SP20_INFO beam_info;                      
  BEAM_DATA data;
  T_MDB_ap_beam_header  old_beam_header;
  FILE  *file1;
  unsigned char dati[MAX_DIM];
  struct tm data_nome;

  int el_num,az_num,i;
  int new_az_num;
  // printf(" devo leggere %s\n",nome_file);
  clear_dbp(vol_pol_locale,nbeam_elev);
  if(access(nome_file,R_OK) == -1) return NO_OK;
  file1=fopen(nome_file,"rb");
  if(file1 == NULL )
  {
    fclose(file1);
    printf("errore apertura file \n");
    return NO_OK;
  }
 i = ReadHeaderSP20toMDB(old_data_header, file1);

  if( i == NO_OK)
  {
    fclose(file1);
    printf("errore lettura dati\n");    
    return NO_OK;
  }

  /*--------
    ATTENZIONE PRENDO LA DATA DAL NOME DEL FILE
    -------*/
   old_data_header->norm.maq.acq_date=get_date_from_name(old_data_header,&data_nome,nome_file);
 
//
    int kk;
kk=0;

   while(1)
  {
kk++;
    if(read_ray_SP20(&old_beam_header,dati,file1,tipo_dati)==NO_OK) break;
    el_num = elevation_index_MDB(old_beam_header.teta);
 
    if(el_num < NEL && old_beam_header.alfa < 4096)
    {
      az_num = azimut_index_MDB(old_beam_header.alfa);
      fill_beam(&vol_pol_locale[el_num][az_num],az_num,el_num,old_beam_header,dati,nbeam_elev);
      if(az_num*0.9-old_beam_header.alfa*FATT_MOLT_AZ < 0.)
        {
	  new_az_num = (az_num +1) %400;
	  fill_beam(&vol_pol_locale[el_num][new_az_num], new_az_num,el_num,old_beam_header,dati,nbeam_elev);
        }
      else if(az_num*0.9-old_beam_header.alfa*FATT_MOLT_AZ > 0.)
        {
	  new_az_num = (az_num -1+400) %400;
	  fill_beam(&vol_pol_locale[el_num][new_az_num], new_az_num, el_num, old_beam_header,dati,nbeam_elev);
        }
     }
//      printf(" ---- %d  \n",nbeam_elev[el_num]);
  }                               //end while
  fclose(file1);
   
  return OK;

}	                                 //end funzione read_dbp_dbp(nome_file)


/*==================================================*/
void conv_ray(float *f_ray, int m, BEAM_HD_SP20_INFO *beam_info, BEAM_DATA *data)
/*==================================================*/
{
  int p;
  float min_zeta;
 
  switch(beam_info->Z_range)
  {
    case 'M' :
      min_zeta=MIN20_Z;
      break;
    case 'L':
      min_zeta=MIN10_Z;
      break;
    case 'H':
      min_zeta=MIN30_Z;
      break;
    case 0 :
      min_zeta=MIN20_Z;
      break;
    case 1:
      min_zeta=MIN10_Z;
      break;
    case 2:
      min_zeta=MIN30_Z;
      break;
  }
  for (p=0; p<beam_info->cell_num; p++) 
    switch (m)
      {
      case INDEX_Z :     // Riflettività Z
	  f_ray[p]=data->beam[p]*RANGE_Z/255. + min_zeta;
	  break;
      case INDEX_ZDR :     // Riflettività differenziale ZDR
	  f_ray[p]=data->beam[p]*RANGE_ZDR/255. + MIN_ZDR;
	  break;
      case INDEX_V :     // Velocità V
	  if (data->beam_w[p] == -128) data->beam_w[p] = -127;
	  if ( beam_info->PRF == 'S')
	    f_ray[p] = data->beam_w[p] * RANGE_V / 127.*.5;
	  else
	    f_ray[p] = data->beam_w[p] * RANGE_V2 / 127.*.5;
	  break;
      case INDEX_SV :    // Spread - Sigma V
	  f_ray[p]=data->beam[p]*RANGE_SIG_V/RANGE_BYTE +MIN_SIG_V;
	  break;
      }
}


/*===============================================*/
void convert_format_beam(BEAM_HD_SP20_INFO *beam_info, T_MDB_ap_beam_header *old_beam_hd, int tipo_dati)
/*===============================================*/
{
  old_beam_hd->max_bin=beam_info->cell_num;
  old_beam_hd->teta=beam_info->elevation/FATT_MOLT_EL;
  old_beam_hd->alfa=beam_info->azimuth/FATT_MOLT_EL;
  old_beam_hd->tipo_gran=tipo_dati;

  return ; 
}

int DefIndiceDati(int tipo_dati, BEAM_HD_SP20_INFO *beam_info)
{
  int i, pos=0;

  // verifico che la grandezza cercata sia presente nel file in uso  
  if (!beam_info->flag_quantities[tipo_dati]) return -1;
  for (i=0;i<tipo_dati;i++)
    if (beam_info->flag_quantities[i])pos++;
  return pos;
}

/*===============================================*/
int ReadBeamSP20toMDB(T_MDB_ap_beam_header *old_beam_hd, unsigned char *dati, FILE *fp, int tipo_dati)
/*===============================================*/
{
  BEAM_HD_SP20_INFO beam_info;                    
  int offset;
  int indice_dati;

  beam_info.valid_data=0;
  while(!beam_info.valid_data){
    if (fread(header_beam_spc_char, sizeof(header_beam_spc_char), 1, fp) != 1)
      return NO_OK;
    decode_header_sp20(header_beam_spc_char, &beam_info);
//  convert_format_beam(&beam_info, old_beam_hd,tipo_dati);
//  printf( "#### %d %d @@@@",old_beam_hd->teta,old_beam_hd->alfa); 
    indice_dati=DefIndiceDati(tipo_dati,&beam_info);
    if (indice_dati == -1) return -1;

    fseek(fp,indice_dati*beam_info.cell_num,SEEK_CUR);
    if(fread(dati, beam_info.cell_num, 1, fp)!= 1)
      return NO_OK;
    offset=(strlen(beam_info.quantities)-indice_dati-1)*beam_info.cell_num;
    fseek(fp,offset,SEEK_CUR);
  }

  convert_format_beam(&beam_info, old_beam_hd,tipo_dati);
  return OK;
}


/*===============================================*/
int ReadBeamSP20(BEAM_HD_SP20_INFO *beam_info, unsigned char *dati, FILE *fp, int tipo_dati)
/*===============================================*/
{
  int offset;
  int indice_dati;

  beam_info->valid_data=0;
  while(!beam_info->valid_data){
    if (fread(header_beam_spc_char, sizeof(header_beam_spc_char), 1, fp) != 1)
      return NO_OK;
    decode_header_sp20(header_beam_spc_char, beam_info);
    
    indice_dati=DefIndiceDati(tipo_dati,beam_info);
    if (indice_dati == -1) return -1;

    fseek(fp,indice_dati*beam_info->cell_num,SEEK_CUR);
    if(fread(dati, beam_info->cell_num, 1, fp)!= 1)
      return NO_OK;
    offset=(strlen(beam_info->quantities)-indice_dati-1)*beam_info->cell_num;
    fseek(fp,offset,SEEK_CUR);
  }

  return OK;
}


/*===============================================*/
int read_ray_SP20(T_MDB_ap_beam_header *old_beam_hd, unsigned char *dati, FILE *fp, int tipo_dati)
/*===============================================*/
{
  T_MDB_ap_beam_header  beam_header_ray[3];
  
  unsigned char dati_ray[3][MAX_DIM];
  long position;
  int i,status;


  if( (status=ReadBeamSP20toMDB(old_beam_hd,dati,fp,tipo_dati))== OK)
  {
    if (old_beam_hd->tipo_gran != tipo_dati)
    {
      fclose(fp);
      _Exit(1);
    }
    return OK;
  }
  if (status == -1) printf (" Il file dati non contiene la grandezza cercata \n");
  return NO_OK;
}                          // end funzione 


          /*-------------------------------------------
	    | identifico a quale ppi appartiene il beam |
	    -------------------------------------------*/
/*===============================================*/
int elevation_index_MDB(short el)
/*===============================================*/
{
  int i;
  for (i=0; i<NEL; i++)
    if(el >= (elev_array[i]-6) && el < (elev_array[i]+5)) return i;
  return NEL;
}



/*----------------------------
 | identifico azimut del beam |
  ----------------------------*/
/*===============================================*/
int azimut_index_MDB(short az)
/*===============================================*/
{
  int i;
  float azimut;
  azimut = az*FATT_MOLT_AZ / .9;
  if(azimut - (int)azimut  <= .5)
    i = (int)azimut % 400;
  else
    i = ((int)azimut +1) % 400;
  return i;
}



             /*---------------------------------------------
	     | FUNCTION: clear_dbp                         |
	     |                                             |
	     ---------------------------------------------
	     | Serve per ripulire il volume dati di lavoro |
	     | Pulisce anche il vettore del numero di      |
	     | raggi trovato                               |
	     ---------------------------------------------*/
/*===============================================*/
void clear_dbp(
        struct VOL_POL vol_pol_locale [][NUM_AZ_X_PPI],
        int nbeam_elev[])
/*===============================================*/
{
  memset(vol_pol_locale,0,sizeof(*vol_pol_locale));
  memset(nbeam_elev,0,sizeof(*nbeam_elev));
  return;
}


/*====================================================================*/
void fill_beam(struct VOL_POL *raggio, int az_num, int el_num, T_MDB_ap_beam_header  old_beam_header, unsigned char *dati, int nbeam_elev[])
/*====================================================================*/
{
	int i;

	if(raggio->flag == 0)
        {
	  for(i=0; i< old_beam_header.max_bin; i++)
	  {
	    if(dati[i])
	      raggio->ray[i] = dati[i];
	    else
	      raggio->ray[i] = 1;
	  }
  	  nbeam_elev[el_num]++;
        }
	else
	  for(i=0; i< old_beam_header.max_bin; i++)
	    if(raggio->ray[i]<dati[i])
	      raggio->ray[i]=dati[i];
	raggio->flag=1;
	raggio->b_header.alfa =(short)(az_num*.9/FATT_MOLT_AZ);
	raggio->b_header.teta = elev_array[el_num];
	raggio->alfa_true =old_beam_header.alfa;
	raggio->teta_true =old_beam_header.teta;
	raggio->b_header.tipo_gran = old_beam_header.tipo_gran;
	raggio->b_header.max_bin = old_beam_header.max_bin;
}

