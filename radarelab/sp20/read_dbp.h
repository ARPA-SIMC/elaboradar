//#define MAX_BIN 460
/*per uniformità lo definisco come 
negli include di cum_bac pari a 512*/
#define MAX_BIN 512

/*---------------------
 | Function prototypes |
  ---------------------*/
extern unsigned short swap2();
extern unsigned int swap4_uns();
extern unsigned int swap4();
extern char * readchar();
extern unsigned short readushort();
extern unsigned int readlong();
extern unsigned int readlong_unswap();
//extern int ReadHeader();
//extern int ReadBeam();
//extern void Read_DBP();
extern int ReadStructureDBP();

extern int read_dbp_SP20_to_DBP(); 
/*----------------------------
 | definizione tipi variabili |
  ----------------------------*/
struct BEAM
{
     T_MDB_ap_beam_header head_beam;
     unsigned char beam[MAX_BIN];
     char flag;
};

struct PPI
{ 
   struct BEAM beam[400];
};

struct DBP
{
   T_MDB_data_header head;
   struct PPI ppi[20]; 
   int n_beam; 
   int flag;
};










































































