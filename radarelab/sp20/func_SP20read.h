/**
 *  @file
 *  @ingroup radarelab 
*/
/*----------------------------------------------------------------------------*/
/*	INCLUDE file						       	      */
/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <unistd.h>
#include <radarelab/sp20/struct_SP20read.h>
#include <radarelab/sp20/range1.h>
#include <radarelab/sp20/radar.h>
#include <radarelab/sp20/meteo_db.h>
#include <radarelab/sp20/read_dbp.h>
#include <math.h>

#define NUM_AZ_X_PPI 400 // numero di raggio per ogni PPI 

struct VOL_POL {
  T_MDB_ap_beam_header b_header;
  unsigned char ray[MAX_DIM];
  char flag;
  short alfa_true, teta_true;
};


/*----------------------------------------------------------------------------*/
/*	FUNCTION PROTOTYPE per la lettura dell'header e del beam  nel nuovo   */
/*        formato SP20 e la conversione nel vecchio MDB                       */
/*      Utilizzato dalle funzioni ReadHeaderSP20toMDB e ReadBeamSP20toMDB     */
/*----------------------------------------------------------------------------*/

int read_dbp_SP20(char *nome_file, struct VOL_POL vol_pol_locale[][NUM_AZ_X_PPI], T_MDB_data_header *old_data_header, int tipo_dati, int nbeam_elev[]);
int SP20ReadHeader();
int SP20ReadBeam();
int read_ray_SP20(T_MDB_ap_beam_header *old_beam_hd, unsigned char *dati, FILE* fp, int tipo_dati);
int elevation_index_MDB(short el);
int azimut_index_MDB(short az);
void conv_ray(float *f_ray, int m, BEAM_HD_SP20_INFO *beam_info, BEAM_DATA *data);
void convert_format_beam();
void clear_dbp(struct VOL_POL vol_pol_locale [][NUM_AZ_X_PPI], int nbeam_elev[]);
void fill_beam(struct VOL_POL *raggio, int az_num, int el_num, T_MDB_ap_beam_header old_beam_header, unsigned char* dati, int nbeam_elev[]);
int DefIndiceDati();
int ReadBeamSP20toMDB();
int ReadHeaderSP20toMDB(T_MDB_data_header *old_header, FILE *fp);
time_t get_date_from_name(T_MDB_data_header *old_data_header, struct tm *tm_date_ext, const char* nome_file);
BEAM_HD_SP20_INFO * decode_header_sp20(unsigned char* hd_char, BEAM_HD_SP20_INFO *beam_info);
void decode_header_DBP_SP20 (HD_DBP_SP20_RAW *hd_raw, HD_DBP_SP20_DECOD *hd_decod);
void PrintHeaderDBP(HD_DBP_SP20_DECOD *hd_decod);
void PrintHeader(BEAM_HD_SP20_INFO *beam_info);
BEAM_HD_SP20_INFO * decode_header_sp20_date_from_name(unsigned char *hd_char,BEAM_HD_SP20_INFO *beam_info,char *nome_file);
void convert_format(HD_DBP_SP20_DECOD *hd_decod, BEAM_HD_SP20_INFO *beam_info, T_MDB_data_header *old_header);
