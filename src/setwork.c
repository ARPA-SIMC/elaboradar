#include <setwork.h>

#include <stdlib.h> 
#include <errno.h> 
#include <string.h> 



int setwork(char *sito)
{  
  char *lv,*lvad,*v0,*v0ad,*vh,*vhad,*logv;
  int ier;

  /*settaggio var amb tempo GTM*/  
  if (getenv("TZ")==NULL)    ier=putenv("TZ=GTM");
  /*settaggio var amb files e directories di lavoro nel caso non siano settate*/		
  if (getenv("DIR_OUT_PP_BLOC")==NULL)  ier=putenv("DIR_OUT_PP_BLOC=./");
  if (getenv("OUTPUT_Z_DIR")==NULL)  ier=putenv("OUTPUT_Z_DIR=./");
  if (getenv("OUTPUT_RAIN_DIR")==NULL)  ier=putenv("OUTPUT_RAIN_DIR=./");
  if (getenv("BACINI_DIR")==NULL)  ier=putenv("BACINI_DIR=./");
  if (getenv("OUTPUT_Z_LOWRIS_DIR")==NULL)  ier=putenv("OUTPUT_Z_LOWRIS_DIR=./");
  if (getenv("DIR_QUALITY")==NULL)  ier=putenv("DIR_QUALITY=./");
  if (getenv("LISTA_FILE")==NULL)  ier=putenv("LISTA_FILE=lista_dbp");
  if (getenv("LAST_FILE")==NULL)   ier=putenv("LAST_FILE=anap_last");
  if (getenv("LOG_FILE")==NULL)  ier=putenv("LOG_FILE=log_file");
  if (getenv("TEST_VPR")==NULL)  ier=putenv("TEST_VPR=test_vpr");
  if (getenv("ANAP_STAT_FILE")==NULL)  ier=putenv("ANAP_STAT_FILE=stat_anap");
  if (getenv("BLOC_STAT_FILE")==NULL)  ier=putenv("BLOC_STAT_FILE=stat_bloc");
  if (getenv("ELEV_STAT_FILE")==NULL)  ier=putenv("ELEV_STAT_FILE=stat_elev");
  if (getenv("NOME_SIDME")==NULL) ier=putenv("NOME_SIDME=sidme.list");
  if (getenv("MATRICE_BACINI")==NULL) ier=putenv("MATRICE_BACINI=grid_bac.dat");
  if (getenv("BACINI_AUS_FILE")==NULL) ier=putenv("BACINI_AUS_FILE=Bac_aus_file");
  if (getenv("BACINI_BOLOGNA")==NULL) ier=putenv("BACINI_BOLOGNA=bacini_xdr_tmp");
  if (getenv("BACINI_HISTORY_FILE")==NULL) ier=putenv("BACINI_HISTORY_FILE=Bacini.History");

  if (!strcmp(sito,"SPC")){
  lv="LAST_VPR=last_vpr_SPC";
  lvad="LAST_VPR_RES=last_vpr_GAT";
  v0="VPR0_FILE=vpr_SPC";
  v0ad="VPR0_FILE_RES=vpr_GAT";
  vh="VPR_HEATING=vpr_heat_SPC";
  vhad="VPR_HEATING_RES=vpr_heat_GAT";
  logv="LOG_VPR=log_vpr_SPC";
  }
  if (!strcmp(sito,"GAT")){
  lv="LAST_VPR=last_vpr_GAT";
  lvad="LAST_VPR_RES=last_vpr_SPC";
  v0="VPR0_FILE=vpr_GAT";
  v0ad="VPR0_FILE_RES=vpr_SPC";
  vh="VPR_HEATING=vpr_heat_GAT";
  vhad="VPR_HEATING_RES=vpr_heat_SPC";
  logv="LOG_VPR=log_vpr_GAT";
  }

  if (getenv("LAST_VPR")==NULL) ier=putenv(lv);
  if (getenv("LAST_VPR_RES")==NULL) ier=putenv(lvad);
  if (getenv("VPR0_FILE")==NULL) ier=putenv(v0);
  if (getenv("VPR0_FILE_RES ")==NULL) ier=putenv(v0ad);
  if (getenv("VPR_HEATING")==NULL) ier=putenv(vh);
  if (getenv("VPR_HEATING_RES")==NULL) ier=putenv(vhad);
  if (getenv("LOG_VPR")==NULL) ier=putenv(logv);

  return 0;
}
