/**
 *  @file
 *  @ingroup radarelab 
*/
/*====================================================================*/
/* Contiene la definizione delle strutture utilizzate in SP20read.c;  */ 
/* hd_DBP_functions.c e beam_functions.c                              */
/*====================================================================*/



#include <time.h>

#define MAX_CELL_NUM 1024    /* ?? = MAX_DIM ? */





struct SP20_BEAM_HEADER
{
  char flag[3]; 		/* bytes 1-3 */
  char az_LSB;			/* byte 4 */
  char az_MSB;			/* byte 5 */
  char el_LSB; 			/* byte 6 */
  char el_MSB;			/* byte 7 */
  char raw_LSB; 		/* byte 8 */
  
  unsigned pulse:2;
  unsigned scan_mode:2;
  unsigned          :1;
  unsigned ctrl_id:2;
  unsigned raw_MSB:1;		/* byte 9 */
  
  unsigned Z:1;
  unsigned Zdr:1;
  unsigned V:1;
  unsigned dualPRF:1;
  unsigned rbin_size:3;
  unsigned sigmaV:1; 		/*byte 10 */
  
  char rbins_LSB;		/* byte 11 */
  
  unsigned RSP_f:1;
  unsigned RSP_SRV_f:1;
  unsigned RSP_RTX_f:1;
  unsigned freq:3;
  unsigned rbins_MSB:2;  	/* byte 12 */
  
  unsigned fltr_type:1;
  unsigned vld_data:1;
  unsigned fltr_num:4;
  unsigned clttr_corr:1;
  unsigned fltr_enable:1; 	/* byte 13 */
  
  char pls_num_LSB;      	/* byte 14 */
  
  unsigned Z_range:2;
  unsigned servo:1;
  unsigned :3;
  unsigned pls_num_MSB:2; 	/* byte 15 */
  
  char real_pwr;		/* byte 16 */
  char eval_pwr;		/* byte 17 */
  char cents;			/* byte 18 */
  char second;			/* byte 19 */
  char minute;			/* byte 20 */
  char hour;			/* byte 21 */
  char week_day; 		/* byte 22 */
  char day;			/* byte 23 */
  char month;			/* byte 24 */
  char year;			/* byte 25 */
  
  unsigned :8;			/* byte 26 */
  
  char angle_offset_MSB;	/* byte 27 */
  char angle_offset_LSB;	/* byte 28 */
  
  char spare[12];		/* bytes 29-40 */
};
typedef struct SP20_BEAM_HEADER SP20_BEAM_HD;
extern unsigned char header_beam_spc_char[40];

   
struct HEADER_BEAM_SP20_INFO
{
  struct tm tm_date;
  time_t time;
  char date[25];
  char valid_data;
  short real_power;
  short eval_power;
  char scan_mode;
  float azimuth;
  float elevation;
  char quantities[5];
  char flag_quantities[4];
  char Z_range;
  char pulse_lenght;
  short pulse_number;
  char cell_size;
  short cell_num;
  char channel;
  char controller[4];
  char PRF;
  char filtro_clutter;
  char tipo_filtro;
  char anticlutter_mappato;
  char corr_clutter;
};
typedef struct HEADER_BEAM_SP20_INFO BEAM_HD_SP20_INFO;	  


struct HEADER_DBP_SP20
{
  char f11;
  char f12;                    /* 1 frame  */
  char f13;
  char f14;
  
  char f21;
  char f22;
  char f23;
  unsigned radZ:1;
  unsigned radZDR:1;
  unsigned radV:1;
  unsigned radSV:1;              /* 4 byte del 2 frame */
  unsigned risoluzione:3;
  unsigned vel_ambigua:1;
  
  char f31;
  char f32;                  /* 3 frame  */
  char f33;
  char f34;
  
  char f41;
  char f42;
  unsigned blank:6;          /* 4 frame */
  unsigned dinamica_Z:2;
  char f44;
  
  unsigned char frame [4][4];      /* ultimi 4 frame  */
  
  
  
  char num_el; 
  char corr_pot;
  char spare[2];
  short ele[40];
};
typedef struct HEADER_DBP_SP20 HD_DBP_SP20;





/*anto*/

struct HEADER_DBP_SP20_RAW
{
  unsigned char frame[8][4];
  char num_ele;
  char corr_pot;
  char spare[2];
  short ele[40];
};
typedef struct HEADER_DBP_SP20_RAW HD_DBP_SP20_RAW;


struct HEADER_DBP_SP20_DECOD
{
  char scan_mode;
  float rot_vel;
  char pulse_lenght;
  float el_min,el_max;
  float az_min,az_max;
  char flag_quantities[4];
  char cell_size;
  char Dual_PRF;
  char filtro_clutter;
  char tipo_filtro;
  char anticlutter_mappato;
  char corr_clutter;
  char Z_range;
  char calibration_factor;
  short pulse_number;
  char stalo_code;
  char controller[4];
  char num_ele;
  char corr_pot;
  float ele[40];
};

typedef struct HEADER_DBP_SP20_DECOD HD_DBP_SP20_DECOD;

union TYPE_OF_DATA
{
  unsigned char beam[MAX_CELL_NUM];
  char beam_w[MAX_CELL_NUM];
};
typedef union TYPE_OF_DATA BEAM_DATA;






