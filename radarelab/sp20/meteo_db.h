/**
 *  @file
 *  @ingroup radarelab 
*/
/*---------------------------------------------------------------------------*
 * NAME/VERSION : meteo_db.h                                                   * 
 *                                                                           * 
 * C.I. :                                                                    *
 *                                                                           * 
 * AUTHOR : Michele Dassisti                                                 * 
 *                                                                           * 
 *                                                                           * 
 * DESCRIPTION :                                                             * 
 *                                                                           * 
 *    Yet implements a new meteo_db architecture whose main features are:    * 
 *                                                                           * 
 *    1. index file (old e_file) plus a data file : dense index and main file* 
 *       architecture.                                                       * 
 *                                                                           * 
 * 2. The index file is organized in two queue double lists of record (idx_rec) 
 *    , the available and free lists. Each record contains a couple of       * 
 *     pointers to the next and previous record in the queue. The record       
 *     information part embeds the most important attributes that a meteo       
 *     map can have, among them the starting point offset value of the map   * 
 *     data in the associated data file.                                     * 
 *                                                                           * 
 * 3. The data file has in its top position the related idx rec plus, if 
 *    necessary, a certain amount of more attribute information called       * 
 *    expansion data. After that header there's the map meteo data.          * 
 *                                                                           * 
 * 4. Extensive bit field data coding has been used in order to minimize     * 
 *    disk space use.                                                        * 
 *                                                                           * 
 * 5. Data file name are created at run time as soon as needed elaborating   * 
 *    the suitable attribute information in the idx_rec. This for saving space 
 *    on disk.                                                               * 
 *                                                                           * 
 * 6. Access to an item by attribute value requires a sequential search      * 
 *    along the pointers chain.                                              * 
 *                                                                           * 
 * 7. Insertion requires either a move of a idx rec from free list to the    * 
 *    available list bottom position or a deletion of the oldest  idx rec    * 
 *    in available list at top position.                                     * 
 *                                                                           * 
 * 8. Deletion requires a move from the available list to the free one plus 
 *    a data file deletion.                                                  * 
 *                                                                           * 
 * 9. A series of idx rec status flag has been implemented for saving maps   * 
 *    from deletion and  for automatic backup.                                   * 
 * CODING DATE : 2/5/91                                                      * 
 *                                                                           * 
 * MASTER LIBRARY TRANSFERT DATE :   3/5/91                                  * 
 *                                                                           * 
 *                                                                           * 
 *                                                                           * 
 *---------------------  MODIFICATION HISTORY  ------------------------------* 
 *                                                                           * 
 * MODIFICATION DATE :                                                       * 
 *   		3-jun-1991 F. Di Pace
 *                                                                           * 
 * MODIFICATION AUTHOR :                                                     * 
 *                                                                           * 
 * MODIFICATION DESCRIPTION :                                                * 
 *                                                                           * 
 *                                                                           * 
 *                                                                           * 
 * RATIONALE :                                                               * 
 *                                                                           * 
 *                                                                           * 
 *---------------------------------------------------------------------------*/

#ifndef METEO_DB_DEFINITION
#define METEO_DB_DEFINITION


#define MDB_SIGN "MDB_2.0"

/* open MDB MODES */

#define  UPDATE		1L
#define  WRITE   	UPDATE
#define  READ    	2L      
#define  CREATE		3L 

#define  MDB_UPDATE	1L 
#define  MDB_WRITE   	UPDATE
#define  MDB_READ    	2L      
#define  MDB_CREATE	3L 
#define  MDB_DELETE	4L 
#define  MDB_KEEP	5L
#define  MDB_INSERT	6L
#define  MDB_RECOVERY	7L 

#define MDB_PAR_SYS_KEEP 1L  /* usato per i dbp di declutter*/
#define MDB_PAR_KEEP 	 2L  /*usato per rat,utr,clb ecc.*/
#define MDB_PAR_BUSY     3L  /*usato per il backup su nastro*/
#define MDB_PAR_BACKUP   4L 

/* MDB_ compression falg */
#define MDB_ENCODE     0L 
#define MDB_DECODE     1L 
  

/* MDB ERRROR CODES */

#define   E_UNABLE    -1
#define   E_LOCK      -2
#define E_DMODE         -3 /* return by del fil */
         
/* MDB_FUNCTIONS command codes */

#define NEWEST          1 
#define OLDEST          2 
#define ACC_by_OFFSET	3
#define ACC_by_NAME	4 


                        
/* stato record */
#define	DELETED	0  /* e_rec  is logically deleted adn data file phisically erased */                   
#define UNLINKED 1 /* e_rec hasn't associate data file -- start up condition */		   
#define LINKED 	2  /* e_rec has a link with a data file */

#define DBP_GRAND_Z 1
#define DBP_GRAND_D 2
#define DBP_GRAND_V 4
#define DBP_GRAND_S 8

	/* Tipi ASM */
#define ASM_AZ 1
#define ASM_VER 2
#define ASM_RAD 4
	/* piani rife HVMI */
#define HVMI_XZ 1
#define HVMI_YZ 2
#define HVMI_XY 4
	/* tipo estrazione SRT/ART */
#define MOD_INT_AUT 0
#define MOD_INT_MAN 1

	/* piano riferimento */

#define PLAN_XZ 0 
#define PLAN_YZ 1 
#define PLAN_XY 2 
#define PLAN_Z  3 

	/* tipo di riempimento */
#define FILL_RAS 0 
#define FILL_NUM 1   
	/* tipo di compression */
#define NO_COMPRESSION 			 0 
#define HUFFMAN_COMPONENTS_COMPRESSION	 1         


                                    /*---------------------------------*
                                     * PARAMETERS DEFINITION           *
                                     *---------------------------------*/
                                     /*---------------------------------*      
                                     * STUFFS TYPE DEFINITION          *
                                     *---------------------------------*/
                
                                    /*---------------------------------*
                                     * data_files types                *
                                     *---------------------------------*/
#define MDB_max_e_type 9
typedef enum {DBP_IDX,               /* 0 */
	      RT_PROD_IDX,           /* 1 */
	      SEC_PROD_IDX,          /* 2 */
              INT_PROD_IDX,          /* 3 */
	      DECL_MAP_IDX,          /* 4 */
              SENS_DATA_SIAP_IDX,    /* 5 */
	      SENS_DATA_CAE_IDX,     /* 6 */
	      RAT_IDX,               /* 7 */                   
              CDS_A_SET_IDX,         /* 8 */
              UTRM_IDX                /* 9 */
	     } T_MDB_e_set_name;    

#define RASTER_IDX INT_PROD_IDX
                        
/**** e_file item number   ***/

#define 	TOT_DBP                 10 /*dbp completi nel db     */
#define 	E_DBP_NUM     		TOT_DBP*4
                                    /* num file nel dbp 1 per grande*/
#define 	E_RTP_NUM               40    
#define       	E_SEP_NUM               40              
/* num. prod integrati = cmm+srt+art number */
#define 	E_INP_NUM               10   
#define 	E_DEC_NUM               56   
#define 	E_SEN_NUM               50   
#define 	A_SET_NUM               30  /* num of a_set files from cds*/
#define 	E_PRO_NUM               E_RTP_NUM + E_SEP_NUM                   
#define 	E_RAT_NUM               200            
#define 	E_UTRM_NUM              200            
#define 	E_ARCHIVE_NUM           50            

/* dimensione file dei prodotti componenti un dato prodotto integrato */




/* data space parms */

#define 	QUANT_NUM               4      /* z,d,s,v           */
#define 	ONE_QUANT_DBP_SIZE      8  * 1024* 1024 /* in bytes */
#define 	ONE_DBP_SPACE           QUANT_NUM * ONE_QUANT_DBP_SIZE
/* deve essere = ONE_DBP_SPACE / 4 = index file size + data file size +
                       acq header infi size (T_MDB_acq_header*/

#define 	DBP_DATA_SPACE          TOT_DBP * ONE_DBP_SPACE   /* in bytes */
#define 	RTP_DATA_SPACE          E_RTP_NUM *400 *1024  /* 2  sri da 512 * 512 */
#define 	SEP_DATA_SPACE          E_SEP_NUM *400*1024       /* in bytes */
#define 	INP_DATA_SPACE          E_INP_NUM *100*1024       /* in bytes */
#define 	DEC_DATA_SPACE          1                  /* in bytes */
#define 	SEN_DATA_SPACE          1                  /* in bytes */
#define 	CDS_DATA_SPACE          10 * 8 *30         /* in bytes */
#define 	PRO_DATA_SPACE          RTP_DATA_SPACE + SEP_DATA_SPACE   
#define 	RAT_DATA_SPACE          E_RAT_NUM * 100 * 1024  
#define 	UTRM_DATA_SPACE         E_UTRM_NUM * 10 * 1024 +\
					E_ARCHIVE_NUM * 600 *1024
/* num record acquisizione per stazione */
#define 	SIAP_SEN_REC_NUM        72            
#define 	CAE_SEN_REC_NUM         144           
                                        

#define         DATA_REC_SIZE           64*1024


/* acquisition el number                */
#define 	EL_MAX_NUM          20             
                                           
#define 	MAX_BEAM_SIZE	    1024


/* max comp in a multifile */
#define    MAX_COMP		    60
#define    MULTI_COMP_MAX_NUM       MAX_COMP /* 4*4 */ 
#define    MULTI_FILE_MAX_NUM       4 
#define    MAX_CLASS_TYPE 	    30
#define	   NULL_INDEX		    0

                                        /*------------------------------
                                         * post el counters array     *
                                         * sizes. used for sep and    *
                                         * inp  .                     *
                                         *----------------------------*/ 

#define 	SEP_NUM                 14    /* different sep+ intp num */
#define 	ACQ_NUM                 999   /* max aqc count val in mod 999*/

                  




                                        /*------------------------------
                                         * sel_sri.c start sri set mark val
                                         * separes differant sri-name-sets
                                         *----------------------------*/ 


#define         START_OF_SET_MARK       0xff



/*** max num of componemt slot in array struct for integr product ****/
#define         COMP_MAX_NUM            270 



/*** max num of a_set name in meteo_db , cds area file selection  ****/
#define         A_SET_MAX_NUM            10 




 
                                        /*------------------------------
                                         * prod name ->index value    *
                                         * association in array       *
                                         *----------------------------*/ 

/* meteo_db.c data acq timeout **************************/

#define         DBP_DATA_TIMEOUT          5 *60  /*sec */
#define         RTP_DATA_TIMEOUT         5 *60  /*sec */






/* meteo_db.c find function out file name****************/



#define         FIND_FILE                "/tmp/find.tmp"



/* meteo_db.c select_sri  function out file name****************/

#define         SRI_SET_FILE             "/proj/tmp/sel_sri.tmp"
#define         VOL_ACQ_SRI              1           /* selection critaria */
#define         MATCH_SECT_SRI           2           /* selection critaria */


/* meteo_db.c queues conf parms **************************/



#define 	HIGH_DATA_SIZE		1
#define 	LOW_DATA_SIZE		26



/*queues selector coding -- prov. oper. ***********/
#define   QUERY          1
#define   INS_DBP        2   /* high queue */
#define   INS_RTPRO      3   /* high queue */
#define   INS_INTPRO     4
#define   INS_DECL       5
#define   INS_SECPRO     6
#define   INS_SD         7
#define   START          8
#define   STOP           9
#define   CLOSE_FILE     10


/* meteo_db queue msg  parms ****************************/

#define    PROD_INFO_MSG_SIZE       1
#define    FILE_NAME_MSG_SIZE       1
#define    CLOSE_FILE_MSG_SIZE      1

/* part_dbp.c parms *************************/


#define    MAX_NUM_AZ_BEAMS         500 /* num max di bam su una scans cic az*/
                                        /* serveper allocare un'array (n_elev
                                           acq corr , n_max beam in az )   */


/* max num of prod comp in inp_e_rec  ********************/

#define    INP_COMP_MAX_NUM       270 



/* inp prod type in norm_inp.c        ********************/

#define    INP_SRT         PROD_ID_SRT          
#define    INP_ART         PROD_ID_ART


/* rain values from siap net meteo_station ***************/

#define    RAIN_VAL_NUM    6


           
typedef 
 unsigned char     T_MDB_rt_prod_name [8];   /* 1 campo rt prod header da 286 */

        
typedef            
 unsigned char     T_MDB_e_file_name[16];	

typedef 
 unsigned char     T_MDB_d_file_name[16];	

typedef 
 struct                 
 {
	unsigned char     c[12];	
 }T_MDB_s_file_name;	

typedef                                                 
 	struct
	{
	 	unsigned char     c[16];
			
	}T_MDB_file_name; 
                             
typedef 
 unsigned char     T_MDB_file_path[50];	

typedef 
 unsigned char     T_MDB_date[18];	


typedef 
 int              T_MDB_num_date;	
                                        
typedef 
 char              T_MDB_sos_mark;/* separes sri-name-sets relative to 
                                          dillfernt time interval in 
                                          sri_set_file(output of the 
                                          sel_sri_function */    

typedef
 char       T_MDB_cds_annotation[21];
typedef
 char       T_MDB_annotation[52];
                   
/* used by get_name (output) and delete_data_file as input************/
/* list of data files associated to one idx_rec */

typedef struct       
{
   unsigned int	n_comp;
   T_MDB_d_file_name  	name[MULTI_FILE_MAX_NUM];    
   T_MDB_d_file_name  	idx_name[MULTI_FILE_MAX_NUM];    
} T_MDB_name_des;                  

/* makes up a_set file from cds */ 
typedef struct 
{                                  
        int	x1;
	int	y1;
	int	x2;
	int	y2;

}T_MDB_a_set_rec;

typedef                                                     
 char   T_MDB_cds_file[A_SET_NUM * sizeof(T_MDB_a_set_rec)];             
                                                                        
                                    /*---------------------------------*
                                     * EXTERNAL DATA TYPES             *
                                     *---------------------------------*/
                               
                        /**********************************/   
			/* Descrizione file indice e dati */
                        /**********************************/   
              

/***********************************************/
/*			index file descriptors */
/***********************************************/
/* link :4 int */
typedef struct           
{
	T_MDB_num_date	cre_date;
	unsigned int	d_file_size;
                                              
	unsigned int	busy_flag:1;
	unsigned int	keep_flag:1;
	unsigned int	sys_keep_flag:1;
	unsigned int	backup_flag:1;
	unsigned int	exp_flag:1;
	unsigned int 	extr_type:2;       /* 0: real time, 1: no */

	unsigned int	del_status:2; 	/* 
					   0: deleted 
					   1: unlinked;
                                           2: linked;
                                        */
	unsigned int	prod_type:7;
				   	/* 0.. 128 */

	unsigned int	qualif:4;
					/* 
					   1: Z
					   2: Zdr
					   4: V
					   8 sV		per cappi 
				 	*/ 
					/* per DBP grand presenti
					   1: Z                        
					   2: Zdr    
					   4: V
			  		   8 sV		
					*/                    

					/*
					   0: automatica, 1: progra per SRT/ART
                                        */  

	unsigned int 	type_compression:2;      
 
					/* 0: non compresso
					   1: compressione 1
					   2: compressione 2
					   3: compressione 3    
					*/
						
 	unsigned int	counter:10; 

 	unsigned int	radar:7; 
                                            
} T_MDB_link_compr;     
                       
/* derivabili: multifile, acq dependent, integrato, utente, mappa, prodotto */
typedef 
	struct           
	{
		T_MDB_num_date	cre_date;
	  	unsigned int	d_file_size;
		unsigned int	rec_offset;

		unsigned short	counter;
 
		unsigned char	prod_type;        
		unsigned char	prod_class;	/*
						 0	map acq
						 1 	pro acq
						 2 	pro_ut maq
						 3 	map Int
						 4 	pro int der
						 5 	pro
						*/
		unsigned char	qualif;
		unsigned char	type_compression;      
		unsigned char	extr_type;
		unsigned char	radar; 		/* 0 SPC */
	} T_MDB_general;     
          

typedef	struct
	{
		unsigned int scans_type:2;	
				    	/*                     
					  0: oraria
					  1: antioraria
				  	  2: sett azimutale
					  3: sett verticale
			 		*/
		unsigned int rot_vel:5;	
					/* 0..30 g/s step 1*/
 									
		unsigned int grand:4;	
					/* 
					  1: Z   
					  2: Zdr = dual polarization
					  4: V
					  8: sV
					 */

		unsigned int resolution:3;
					/*
					  0:   62.5
					  1:  125
					  2:  250
					  3:  500
					  4: 1000
					  5: 2000 metri
					*/
		unsigned int vel_range:1;
					/*
					  0: 16.5= no prf2
					  1: 49  = si prf2
					 */ 
 
		unsigned int declutter_rsp:1;
 					/* 0: no 1: si 	*/

		unsigned int type_declutter:1;
					/* 0: uniforme 1: mappato */
		unsigned int filter_value:4;	
											/* codice 0..15 */
		unsigned int corr_Z:1;		
					/* 0: no 1: si 	*/
		unsigned int num_imp:10;
					/* 0..1023		*/ 
/* int */

		unsigned int declutter_sw:3;
					/* 0: no,
			 	           1..3: declutter 1..3
					*/
		unsigned int quota_cut_sw:5;
                                                               
		unsigned int el_ini:12;                       
				     	/* 0..4095		*/
		unsigned int el_fin:12;
					/* 0..4095		*/ 
/* int */
		unsigned int imp_duration:2;		
					/*
					  0: lunga= 3.  micros = 300  prf
					  1: media= 1.5 micros = 600  prf
					  2: corta= 0.5 micros = 1200 prf
					*/
		unsigned int coverage:3;   			
					/* 
					   0: 16 km
					   1: 32
					   2: 64
					   3: 128
					   4: 256
					   5: 512
					 */

		unsigned int	acq_count:10;
		unsigned int	acq_duration:10; /* sec *8 */
	    	unsigned int	prec_z:7;        
/* int 3 */
	    	unsigned int num_el:8;	/* 0..256	*/
		/*
			     	unsigned int az_step:4;	
	  			unsigned int el_step:4;  	
		*/
		unsigned int az_ini:12;	       	/* 0..4095  */ 
		unsigned int az_fin:12;			/* 0..4095 */ 

	     	unsigned int	acq_date;

	} T_MAQ_compr;        
/*
dati derivabili:
	copertura effettiva = min(f(durata_impulso,velocita),copertura richiesta)
	precisione=f(numero_impulsi,risoluzione,velocita) 
*/ 	                               

typedef	struct
	{
	     	int		acq_date;    

	   	unsigned short	num_el;
		unsigned short	num_az;

		unsigned short	coverage;			
	   	unsigned short	num_imp;

		unsigned short  az_ini;
		unsigned short  az_fin;

		unsigned short  el_ini;
		unsigned short	el_fin;
               
   		unsigned short  acq_duration;        /* Tempo di esecuzione rat */
		unsigned short	acq_count;

		unsigned short  value[30];  /* valori per scans. sett. 0->4095 */


		unsigned char	grand;	
		unsigned char	declutter_rsp;
		unsigned char	type_declutter;
		unsigned char	filter_value;	

	 	unsigned char	corr_Z;		
		unsigned char	declutter_sw;
		unsigned char	quota_cut_sw;
		unsigned char	imp_duration;		

		unsigned char	resolution;
		unsigned char	rot_vel;	
		unsigned char	az_step;
		unsigned char	el_step; 

		unsigned char	scans_type;	
		unsigned char	vel_range;      
		unsigned char   dual;
		unsigned char   spare[1];

		unsigned char   prec[4];      /* precisioni*/                  
	}T_MAQ;     
                                       

/* MAP : 4 int */
typedef	struct
	{
		unsigned int alt_azim:12;	/* 0..4095                  
					   		0.0879121 gradi
					   		62.5 metri 
	   					*/	
 
		unsigned int quadr:3;	       	/* 0:pieno,
					       	   1,2,3,4
					       	*/	

		unsigned int type_map:1;	/* 0: numerica
						   1: raster
						*/

		unsigned int gray_level:2;	/* 0: 256
						   1: 16
						*/
		unsigned int num_pix_x:14;

		unsigned short num_pix_y;
		unsigned short num_pix_z;
		               
		unsigned char	latitude[3];
		unsigned char	intitude[3];

		unsigned short	map_x_resolution;
		unsigned short	map_y_resolution;
	 	unsigned short	map_z_resolution;
               
	}T_MDB_map_info_compr;
                                
#define T_MDB_base_compr T_MDB_map_info_compr
typedef	struct           
	{
 		unsigned int	map_x_resolution;
		unsigned int	map_y_resolution;
		unsigned int	map_z_resolution;

 		unsigned short	quadr; 	/* TOTAL per immagini raster */
 		unsigned short	alt_azim;

		unsigned short	type_map;
		unsigned short  gray_level;

	 	unsigned short  num_pix_x;
		unsigned short  num_pix_y;

		unsigned short  num_pix_z;
		unsigned short  az_ini;

 		unsigned short	az_fin;
		unsigned short 	az_ini2;

		unsigned short 	az_fin2;

		unsigned char	latitude[3];
		unsigned char	longitude[3];
	

 	}T_MDB_map_info;
               
#define T_MDB_base T_MDB_map_info;
               

typedef
	struct 
	{
		unsigned int 	end_intgr;
	       	unsigned short 	start_az,
				end_az;
		unsigned short	int_intgr;
		unsigned short	minuti_integrati;                
	}T_MDB_map_int;             
	
typedef
	struct
	{
		unsigned int 	end_intgr;
		unsigned short	int_intgr;
		unsigned short	spare;
	} T_MDB_pro_int;    

/* Mappe standard dipendenti dall'acquisizione : cappi, sri, vmi etc  */
typedef
	struct{                
	 	T_MAQ_compr		maq;
		T_MDB_map_info_compr	map;       
	}T_MDB_map_acq_e_rec;                      

#define T_MDB_rtp_e_rec T_MDB_map_acq_e_rec
#define T_MDB_sep_e_rec T_MDB_map_acq_e_rec      

/* Prodotto standard dipendente dall'acquisizione: DBP  ,Declutter_sw  */
typedef
	struct
	{                                          
		T_MAQ_compr    		maq;
	       	T_MDB_s_file_name		s_file;	
	}T_MDB_pro_acq_e_rec;                           

#define T_MDB_dbp_e_rec T_MDB_pro_acq_e_rec                              
#define T_MDB_clu_e_rec T_MDB_pro_acq_e_rec

/* Prodotto utente dipendente da modalita' acquisizione: RAT, MAQ , */
typedef
	struct                                  
	{
		T_MAQ_compr		maq;
	       	T_MDB_file_name		file;	
	}T_MDB_pro_ut_maq_e_rec;                         

#define T_MDB_rat_e_rec T_MDB_pro_ut_maq_e_rec
#define T_MDB_maq_e_rec T_MDB_pro_ut_maq_e_rec

/* Mappe standard integ nel tempo : SRT */
typedef
	struct{ 
		T_MDB_map_int 		int_inf;
		T_MDB_map_info_compr	map;
	}T_MDB_map_int_e_rec;           

#define T_MDB_srt_e_rec T_MDB_map_int_e_rec

/* Prodotto integr tempo e dipendente da descrittore : ART */

typedef
	struct                                  
	{
		T_MDB_pro_int  		int_inf;
	       	T_MDB_file_name		file_dip;	
	}T_MDB_pro_int_dip_e_rec;            

#define T_MDB_art_e_rec T_MDB_pro_int_dip_e_rec 


/* Prodotto standard indipendente dall'acquisizione: SENSORI , areole */
/* Prodotto utente indipendente dall'acquisizione: SAT */
typedef
	struct
	{
	       	T_MDB_file_name		file;	
	}T_MDB_pro_ut_e_rec;         
                                         
typedef    
	struct
	{       
	  T_MDB_link_compr		link;
	  union       
	  {
	    	T_MDB_map_acq_e_rec   	map_acq; 
		T_MDB_pro_acq_e_rec     pro_acq;       
	   	T_MDB_pro_ut_maq_e_rec  pro_ut_maq;
		T_MDB_map_int_e_rec	map_int;
		T_MDB_pro_int_dip_e_rec	pro_int_dip;
	       	T_MDB_pro_ut_e_rec	pro;        
	  }prod;        
	}T_MDB_idx_info;
                    
typedef struct   
{
	unsigned int	  	next;
	unsigned int		prev;
} T_MDB_idx_rec_head;
                                        
/* record del file indice */
typedef struct
	{
		T_MDB_idx_rec_head   head;
	        T_MDB_idx_info    	info;
	} T_MDB_idx_rec;          

/************************/
/* Struttura index file */
/************************/

typedef struct
{
		unsigned int avail_start; 
		unsigned int avail_end; 
		unsigned int free_start; 
		unsigned int free_end; 

	        unsigned int			e_rec_num;
		unsigned int                    tot_d_space;
		unsigned int                    av_d_space; 
		unsigned int                    acq_count;  

}T_MDB_idx_file_header;     

typedef
	struct
	{                  
	       	T_MDB_name_des 		names;	
		char	spare[4-sizeof(T_MDB_name_des)%4];
                      
	 	T_MDB_general 		general;
		char	spare1[4-sizeof(T_MDB_general)%4];

	 	T_MAQ			maq;   
		char	spare2[4-sizeof(T_MAQ)%4];

		T_MDB_map_info		map;       
/*		char	spare3[4-sizeof(T_MDB_map_info)%4];*/

	       	T_MDB_file_name		file_dip;   
		char	spare4[4-sizeof(T_MDB_file_name)%4];

		T_MDB_map_int  		int_inf;
		char	spare5[4-sizeof(T_MDB_map_int)%4];

     	       	T_MDB_annotation	annotation;
		char	spare6[4-sizeof(T_MDB_annotation)%4];
	}T_MDB_norm_e_rec;               

/*DATA TYPES THAT MADE UP THE EXPANSION INFO AFTER THE DATA FILE HEADER  */ 
typedef 
	struct
	{
	  	T_MDB_d_file_name          comp[INP_COMP_MAX_NUM];
	  	char			   spare[INP_COMP_MAX_NUM];
                char			   spare1[2];
	}T_MDB_srt_exp_attr;    

typedef
	struct
	{
	  	unsigned int	inf[128];       
	}T_MDB_sri_exp_attr;    
	
typedef 
	struct
	{          
		T_MDB_d_file_name          comp[INP_COMP_MAX_NUM];	
	  	char			   spare[INP_COMP_MAX_NUM];
                char 			   spare1[2];
		T_MDB_cds_file 	        cds;/* is the cds selected areas where 
	                                           art comes from*/	
	}T_MDB_art_exp_attr;    

	
#define STATUS_BIT_OFF		0
#define STATUS_BIT_ON		1
				/* Stati dei bit per il campo che definisce   */
				/* lo stato della classe: STATUS_mask         */
typedef struct {
   	unsigned int to_create 	: 1;
   	unsigned int to_destroy 	: 1;
   	unsigned int to_save 		: 1;
   	unsigned int to_load 		: 1;
	unsigned int saved 		: 1;
   	unsigned int modified 		: 1;
   	unsigned int loaded 		: 1;
   	unsigned int cannot_encode 	: 1;
   	unsigned int to_encode 	: 1;
   	unsigned int to_not_decode 	: 1;
   	unsigned int to_append 	: 1;
   	unsigned int resize 		: 3;
	unsigned int spare 		:18;
} STATUS_mask;


#define UNIX_CLASS		1
#define VMS_CLASS		2
#define MSDOS_CLASS		3
				/* Costanti che definiscono l'ambiente dove e'*/
				/* stata creata la classe o il file	      */

typedef struct {		/* Tale struttura definisce qualsiasi classe  */
   				/* che risiede in memoria o in un file	      */
	unsigned int pointer;	/* Puntatore alla posizione della classe:     */
   				/* indirizzo di memoria o offset sul file     */
	int space;		/* Size in byte della memoria allocata per la */
				/* classe				      */
	int size;		/* Size in byte della classe	      	      */
	STATUS_mask status_mask;/* Stato della classe			      */
	unsigned char type_compression;
				/* Eventuale codice del tipo di compressio-   */
				/* ne applicato sulla classe		      */
	unsigned char environ;  /* Ambiente operativo dove e' stato creata la */
				/* classe:				      */
				/* UNIX_CLASS - VMS_CLASS - MSDOS_CLASS       */
        /* 2014-01-15 renamed from class to class_code to make it valid C++ */
	unsigned char class_code;	/* Codice della classe 		              */
	unsigned char index;    /* Indice della componente nel file           */
} CLASS_pointer;


typedef struct {		/*             CLASS_HEADER_CLASS             */
   	unsigned char signature[9];
   				/* Firma di riconoscimento dell'header        */
   				/* costante MDB_SIGN			      */
	unsigned char environ;  /* Ambiente operativo dove e' stato creato il */
				/* file:			              */
				/* UNIX_CLASS - VMS_CLASS - MSDOS_CLASS       */
	unsigned short max_comp;/* MAX_COMP			              */
	unsigned short max_class_type;
				/* MAX_CLASS_TYPE			      */
	unsigned short n_comp;  /* Maggiore indice delle componenti contenute */
				/* nel file    				      */
	CLASS_pointer compDesc[MAX_COMP+1];
} CLASS_header;



typedef                 
	struct
	{     
		CLASS_header 		class_header;
		T_MDB_norm_e_rec	norm;   
	}T_MDB_data_header;

#define T_MDB_acq_header T_MDB_data_header


/*-------------------------------------------------------------*/
/***** strutture dati che arrivano dal 286; linea dati rt ******/
/*-------------------------------------------------------------*/

/* header dei dati polari estratti dal 286 */

typedef struct
		{
			 short 		coverage;
			 char           grand;          
			 char   	vel_range;
			 char           resolution;               
			 char           precision[4];
			 char           date[17];                     
			 short          prf;	
			 char           prf2;	
			 char		imp_duration;
			 char		dual_pol;
			 char		rot_vel;       
			 short          	acq_duration;
			 short          num_az; 
			 short          num_el; 
			 short          start_az;   
			 short          end_az;          
			 short          start_el;
			 short          end_el;  
			 char		scans_type;
			 char		declutter_rsp;
                         char           corr_z;
                         char           filter_value;
			 char           num_prod; 
                         char           declutter_sw;
                         char           el_list[EL_MAX_NUM];
                         char           az_step;   
                         char           el_step;   
		}T_MDB_286_acq_header;    


/* header dei prodotti estratti dal 286 */      
typedef struct
{
		 char     prod_type[8];	/*			*/
		 char     resolution;   /*			*/
		 char     dimension;    /*			*/
		 short    altitude;     /*                      */
		 char     ref_plane;	/*     1:xz, 2: yz, 3:xy, 4:z */
		 char     extr_type;    /* 			*/
		 short	  coverage;     /* Km                   */
}T_MDB_286_pro_header;
               
		/*--------------------------------------*   
		* HEADER dei beam 286		     	*
		*---------------------------------------*/ 

typedef struct
		{
		       	short          teta;        /* 0..184 */
			short          alfa;        /* 0..399 */
			short          tipo_gran;      
			short          max_bin;      

		}T_MDB_ap_beam_header;

	/*-----------------------------------------------------------------*/ 
	/* Pacchetto informativo fornito dal produttore prodotti secondari */   
	/*-----------------------------------------------------------------*/

struct coord_descr          
 {
   char		spare;
   char		gradi;          
   char		primi;
   char		secondi;
};
typedef struct
{           
   unsigned char   prod_type;
   unsigned char   quadr;

   unsigned char   extr_type;
   unsigned char   type_compression;
   unsigned short   gray_level;
  
   unsigned short  alt_azim;     /* quota in metri o az 0:4095 */  
   unsigned short  num_pix_x;      
   unsigned short  num_pix_y;      
   unsigned short  num_pix_z;      

   float	map_x_resolution; 
   float	map_y_resolution;
   float	map_z_resolution;

   struct 	coord_descr latitudine;    
   struct 	coord_descr intitudine;    
}T_MDB_sep_info;        

typedef struct
{           

		unsigned int 	newest_data;
		unsigned short	int_intgr;
		unsigned short	minuti_integrati;
}T_MDB_inp_info;   

	/*---------------------------------------------------------------*/
	/* Pacchetto informativo fornito dal produttore prodotti ART/SRT */ 
	/*---------------------------------------------------------------*/
/* used by auto_srt /man_srt producer */                  

typedef struct  { 
		   	unsigned int		d_file_size; 
			unsigned int		intgr_start_time; 
			unsigned int		intgr_end_time; 
			unsigned int		intgr_int; 
			unsigned int		coverage;
			float			resolution;
                        char                     calc_mode; 
                        char                     prod_type; 
		}T_MDB_inp_header;
               

		/*---------------------------------*/
		/* dbp_index file record           */  
		/*---------------------------------*/

typedef struct
		{     
			int           offset;      
			int           beam_size;      

		}T_MDB_dbp_idx_info;
typedef 
 T_MDB_dbp_idx_info        *T_MDB_idx_array[4]; 

			/* array di puntatori ad array record idx_info */	
                        /* dim = num grandezze */
                                    /*---------------------------------*
                                     * for each couple (sep, dbp_acq_num)
                                     * rhe array gives the last post 
                                     * elaboration counter value; used        
                                     * to create sep data file name.
                                     *---------------------------------*/

			
       
/**** defines components file name for int prod man_srt and art ***/

typedef T_MDB_d_file_name   T_MDB_comp_selection[COMP_MAX_NUM];    


/* used by get_name (output) and delete_data_file as input************/
/* list of data files associated to one idx_rec */


                                                                            
                                


/* used to store areas selection from cds , to calculate art integr prod*/
typedef struct  
{
	T_MDB_link_compr              link;      

}T_MDB_cds_e_rec ;

		    



typedef unsigned char	T_MDB_data_rec [ DATA_REC_SIZE];




                                    /*---------------------------------*
                                     * FIND tipes defs                 *
                                     *---------------------------------*/


typedef   struct
{ 

		    unsigned char                prod_name[5];  
		    unsigned char                num_char_to_comp;
} T_MDB_prod_id;


                                    /*---------------------------------*
                                     * SEL_SRI type defs (for auto_srt)*
                                     *---------------------------------*/


/* sri selection criteria in sel_sri.c */
typedef   struct
{ 
	char    sri_acq_type;    
	short   integration_int;
	short   start_az;         
	short   end_az;          
        
} T_MDB_sri_selection;


typedef   struct    /* tells how to write next sel_rec in update_sri_file proc*/
{ 
	int    p_in;            
} T_MDB_sri_set_file_header;




/* info returned from search in db by sel_sri.c */
typedef   struct
{ 
	short   tot_sri_db_num;     
	short   ok_sri_num;              
	short   ok_sri_sat_num;       

} T_MDB_sri_search_info;



typedef   struct
{ 
		    T_MDB_d_file_name       sri_name;      
		    unsigned short               coverage;         
                    T_MDB_num_date          time;
                    float                        temp_coeff; /* for srt calc*/
} T_MDB_sri_sel_data;



/*** used by a-set_list routinre to reteunn a_set name list in meteo_db */
typedef   struct
{ 
		    T_MDB_d_file_name       name[A_SET_MAX_NUM];      
		    int                          name_num;         
} T_MDB_a_set_list  ;


/* makes up a meteo station file in SEN branch of meteo_db */
typedef   struct
{ 
                    T_MDB_num_date          acq_time;
		    short                        amr;                          
		    short                        rain_value[RAIN_VAL_NUM];
} T_MDB_siap_sen_rec;

/* defines data code to select data in siap rec buffer */
typedef  unsigned char  T_MDB_sen_data_selection [RAIN_VAL_NUM];

                       



/* used in sen data part of meteo_db */
typedef struct
{
		int				p_in;    
		int				p_out; 
		int                            p_end; 
		int                            p_start; 
	 	unsigned int                    e_rec_size;  	
		unsigned int			e_rec_num;
		unsigned int                    tot_d_space;
		unsigned int                    av_d_space; 
		unsigned int                    acq_count;  
}T_MDB_ms_file_header;     

/* parte informazioni del record del file indice 
/typedef struct 
/{
/  	T_MDB_e_link               link; 
/        union {
/               T_MDB_sep_e_rec  	sep;
/               T_MDB_rtp_e_rec     rtp;
/               T_MDB_dbp_e_rec     dbp;
/               T_MDB_art_e_rec	art;
/               T_MDB_srt_e_rec	srt;
/               T_MDB_rat_e_rec	rat;
/              }prod_type;
/                
/} T_MDB_idx_info;
*/
                                         
typedef struct
{
  int		        data_fd[MULTI_FILE_MAX_NUM];  
  int		      	idx_fd[MULTI_FILE_MAX_NUM];               
  int		      	num_fd;

  int			mode;
  T_MDB_data_header    	*head;
  T_MDB_e_set_name     	e_type;
  unsigned int        	comp_read_off[MULTI_COMP_MAX_NUM];       
  unsigned int        	comp_off[MULTI_COMP_MAX_NUM];

}T_MDB_channel;

typedef struct
{
	int	num_fd;
	int	data_fd[MULTI_FILE_MAX_NUM];	  
	int	idx_fd[MULTI_FILE_MAX_NUM];	   

}T_MDB_multi_fd;

/****************************************************************/
/* 	FIND                                                    */
/****************************************************************/

#define    MAX_QUERY   10L
                          
#define    SUCCESS     1L   
#define    FAIL        0L 

#define    NO_ORD      0L 
#define    ACQ_OLDEST  1L 
#define    ACQ_NEWEST  2L  
#define    CRE_OLDEST  4L  
#define    CRE_NEWEST  8L  


#define    NO_TEST     0L   
#define    TEST_EQ     1L 
#define    TEST_NEQ    2L 
#define    TEST_GT     3L     
#define    TEST_GE     4L 
#define    TEST_LT     5L    
#define    TEST_LE     6L 
#define    TEST_LT_GT  7L 
#define    TEST_LT_GE  8L                         
#define    TEST_LE_GT  9L      
#define    TEST_LE_GE  10L                         
#define    TEST_GT_LT  11L      
#define    TEST_GT_LE  12L                         
#define    TEST_GE_LT  13L      
#define    TEST_GE_LE  14L                         


#define PROD_NIL  '\0'   
#define QUADR_NIL '\0'   
#define FILE_NIL  '\0'
#define QUOTA_NIL -1
#define AZ_NIL -1
#define INT_NIL   -1
#define DATE_NIL  0

typedef
	struct
	{
		char	prod[12+1]; 	/*nome ascii del prodotto        */
					/* vuoto non significativo */

		char	quadr[8+1];     /* quadrante */
					/* vuoto non significativo */

	        char    file[14+1];  	/* sens, file ART, etc)*/
					/* vuoto non significativo */

		float	quota; 		/* -4giga non significativo */

		float	azimuth;	/* -4giga non significativo */

		float	int_azimuth[2];	/* -4giga non significativo */
                            
		int	int_intgr; 	/* -4giga non significativo */

		int	extr_type; 	/* -4giga non significativo */

		int	cre_date;	/* 0 non significativo */

		int	x_copertura	;

		int	z_copertura	;

		int	map_x_resolution		;

		int	map_z_resolution  		;

		int	type_map 	;

		int	gray_level 	;

		char	latitude[4];    
		char	intitude[4];                                                                             

		int	acq_date;   	/* 0 non significativo */

	}T_MDB_find_rec;                                                 

typedef
	struct
	{                                                    
		int	prod_match ;      /* EQ - NEQ         */
	      	char	prod[12+1] ;   	  /*nome ascii del prodotto        */
			/*se blank no query*/
                                 
		int	quadr_match ;     /* EQ - NEQ         */
		char	quadr[8+1] ;      /* quadrante "ALL","HL",HR,LL,LR */
			/*se blank no query*/

		int	file_match ;      /* EQ - NEQ         */
	        char    file[14+1] ; 	  /* wildchar (sens, file ART, etc)*/
    			/*se blank no query*/

		int	quota_match ;      /* = > < inta est   */

		float	quota[2]    ; 	   /* intervallo quote */ 	

		int	azimuth_match ;    /* = > < inta est   */

		float	azimuth[2]    ;    /* intervallo azimut */

		int	int_azimuth_match ;    /* = > < inta est   */ 

		float	int_azimuth[2]    ;    /* intervallo azimut */
                                        
		int	extr_type_match   ;    /* EQ -NEQ   */
		int	extr_type  ; 

   		int	integ_match   ;    /* = > < inta est   */
		int	int_intgr[2]  ;    	  	                       
                                      
 		int	creat_match   ;    /* = > < inta est   */
		int	cre_date[2] ;    /* intervallo data creazione */

		int	acq_match     ;    /* = > < inta est   */
		int	acq_date[2]   ;    /* intervallo data acquisizione */

		int	x_copertura_match	;
		int	x_copertura[2]	   	;

		int	z_copertura_match	;
		int	z_copertura[2]		;

		int	map_x_resolution_match;
		int	map_x_resolution;

		int	map_z_resolution_match	;
		int	map_z_resolution  		;

		int	type_map_match	;
		int	type_map 	;

		int	gray_level_match	;
		int	gray_level 	;

		int	coord_match     ;    
		char	latitude[4];    
		char	intitude[4];                                                                             

		int	ord_match     ;                                     
		int	ord_flag      ;                                     
		int	ord_time      ; 
                

	}T_MDB_find_match_rec;                                                 

#endif
