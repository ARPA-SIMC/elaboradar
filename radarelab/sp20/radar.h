/**
 *  @file
 *  @ingroup radarelab 
*/
/*#define ACF 0.043945313 */           /* angle conversion factor for SP20 radar data */
/*#define ACF 0.0879  */

#define ACF ((double) (360.0)/(double)(8192.0)) 

#define FATT_MOLT_EL ((double) 360./(double)4096.)
#define FATT_MOLT_AZ ((double) 360./(double)4096.)

#define HEADER_SIZE 40                /* SP20 radar data header size */
#define BEAM_HEADER_SIZE 40           /* SP20 radar data header size */

//#define NEL 10                        /* N elevazioni */
#define NEL 15                        /* N elevazioni deal261104*/


#define MAX_BIN 512     /* verificare il numero corretto 1024 o 512 */
#define MAX_DIM 512



   struct SP20_HEADER
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
   typedef struct SP20_HEADER SP20_HD;
   
   struct HEADER_INFO
          {
          char date[9];
	  char time[12];
	  char valid_data;
	  short real_power;
	  short eval_power;
	  char scan_mode;
	  float azimuth;
	  float elevation;
	  char quantities[5];
	  char Z_range;
	  char pulse_lenght;
	  short pulse_number;
	  char cell_size;
	  short cell_num;
	  char channel;
	  char controller[4];
          char PRF;
          };
   typedef struct HEADER_INFO HD_INFO;	  

	struct HEADER_FILE
	{
	char frame[8][4];
	char num_el;
	char corr_pot;
	char spare[2];
	short ele[40];
	};
   typedef struct HEADER_FILE FILE_INFO;
