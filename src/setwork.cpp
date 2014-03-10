#include <setwork.h>

#include <stdlib.h> 
#include <errno.h> 
#include <string.h> 
#include <stdio.h>
#include "logging.h"

int setwork(char *sito)
{  
  char *lv,*lvad,*v0,*v0ad,*vh,*vhad,*logv;
  int ier;

  /*settaggio var amb tempo GTM*/  
  if (getenv("TZ")==NULL)    ier=putenv((char *)"TZ=GTM");
  /*settaggio var amb files e directories di lavoro nel caso non siano settate*/		
  if (getenv("DIR_OUT_PP_BLOC")==NULL)  ier=putenv((char *)"DIR_OUT_PP_BLOC=./");
  if (getenv("OUTPUT_Z_DIR")==NULL)  ier=putenv((char *)"OUTPUT_Z_DIR=./");
  if (getenv("OUTPUT_RAIN_DIR")==NULL)  ier=putenv((char *)"OUTPUT_RAIN_DIR=./");
  if (getenv("BACINI_DIR")==NULL)  ier=putenv((char *)"BACINI_DIR=./");
  if (getenv("OUTPUT_Z_LOWRIS_DIR")==NULL)  ier=putenv((char *)"OUTPUT_Z_LOWRIS_DIR=./");
  if (getenv("DIR_QUALITY")==NULL)  ier=putenv((char *)"DIR_QUALITY=./");
  if (getenv("LISTA_FILE")==NULL)  ier=putenv((char *)"LISTA_FILE=lista_dbp");
  if (getenv("LAST_FILE")==NULL)   ier=putenv((char *)"LAST_FILE=anap_last");
  if (getenv("LOG_FILE")==NULL)  ier=putenv((char *)"LOG_FILE=log_file");
  if (getenv("TEST_VPR")==NULL)  ier=putenv((char *)"TEST_VPR=test_vpr");
  if (getenv("ANAP_STAT_FILE")==NULL)  ier=putenv((char *)"ANAP_STAT_FILE=stat_anap");
  if (getenv("BLOC_STAT_FILE")==NULL)  ier=putenv((char *)"BLOC_STAT_FILE=stat_bloc");
  if (getenv("ELEV_STAT_FILE")==NULL)  ier=putenv((char *)"ELEV_STAT_FILE=stat_elev");
  if (getenv("NOME_SIDME")==NULL) ier=putenv((char *)"NOME_SIDME=sidme.list");
  if (getenv("MATRICE_BACINI")==NULL) ier=putenv((char *)"MATRICE_BACINI=grid_bac.dat");
  if (getenv("BACINI_AUS_FILE")==NULL) ier=putenv((char *)"BACINI_AUS_FILE=Bac_aus_file");
  if (getenv("BACINI_BOLOGNA")==NULL) ier=putenv((char *)"BACINI_BOLOGNA=bacini_xdr_tmp");
  if (getenv("BACINI_HISTORY_FILE")==NULL) ier=putenv((char *)"BACINI_HISTORY_FILE=Bacini.History");

  if (!strcmp(sito,"SPC")){
  lv=(char *)"LAST_VPR=last_vpr_SPC";
  lvad=(char *)"LAST_VPR_RES=last_vpr_GAT";
  v0=(char *)"VPR0_FILE=vpr_SPC";
  v0ad=(char *)"VPR0_FILE_RES=vpr_GAT";
  vh=(char *)"VPR_HEATING=vpr_heat_SPC";
  vhad=(char *)"VPR_HEATING_RES=vpr_heat_GAT";
  logv=(char *)"LOG_VPR=log_vpr_SPC";
  }
  if (!strcmp(sito,"GAT")){
  lv=(char *)"LAST_VPR=last_vpr_GAT";
  lvad=(char *)"LAST_VPR_RES=last_vpr_SPC";
  v0=(char *)"VPR0_FILE=vpr_GAT";
  v0ad=(char *)"VPR0_FILE_RES=vpr_SPC";
  vh=(char *)"VPR_HEATING=vpr_heat_GAT";
  vhad=(char *)"VPR_HEATING_RES=vpr_heat_SPC";
  logv=(char *)"LOG_VPR=log_vpr_GAT";
  }

  if (getenv("LAST_VPR")==NULL) ier=putenv(lv);
  if (getenv("LAST_VPR_RES")==NULL) ier=putenv(lvad);
  if (getenv("VPR0_FILE")==NULL) ier=putenv(v0);
  if (getenv("VPR0_FILE_RES ")==NULL) ier=putenv(v0ad);
  if (getenv("VPR_HEATING")==NULL) ier=putenv(vh);
  if (getenv("VPR_HEATING_RES")==NULL) ier=putenv(vhad);
  if (getenv("LOG_VPR")==NULL) ier=putenv(logv);

  return ier;
}

void unsetwork(){
unsetenv((char *)"FIRST_LEVEL_FILE");
unsetenv((char *)"FIRST_LEVEL_DIM_FILE");
unsetenv((char *)"FILE_T");
unsetenv((char *)"VPR_HMAX");;
unsetenv((char *)"FILE_ZERO_TERMICO");
unsetenv((char *)"DIR_STORE_VPR");
unsetenv((char *)"SITO");
unsetenv((char *)"VPR_ARCH");
unsetenv((char *)"FILE_DEM_SPC");
unsetenv((char *)"FILE_DEM_GAT");
                
unsetenv((char *)"TZ");
unsetenv((char *)"DIR_OUT_PP_BLOC");
unsetenv((char *)"OUTPUT_Z_DIR");
unsetenv((char *)"OUTPUT_RAIN_DIR");
unsetenv((char *)"BACINI_DIR");
unsetenv((char *)"OUTPUT_Z_LOWRIS_DIR");
unsetenv((char *)"DIR_QUALITY");
unsetenv((char *)"LISTA_FILE");
unsetenv((char *)"LAST_FILE");
unsetenv((char *)"LOG_FILE");
unsetenv((char *)"TEST_VPR");
unsetenv((char *)"ANAP_STAT_FILE");
unsetenv((char *)"BLOC_STAT_FILE");
unsetenv((char *)"ELEV_STAT_FILE");
unsetenv((char *)"NOME_SIDME");
unsetenv((char *)"MATRICE_BACINI");
unsetenv((char *)"BACINI_AUS_FILE");
unsetenv((char *)"BACINI_BOLOGNA");
unsetenv((char *)"BACINI_HISTORY_FILE");
unsetenv((char *)"LAST_VPR");
unsetenv((char *)"LAST_VPR_RES");
unsetenv((char *)"VPR0_FILE");
unsetenv((char *)"VPR0_FILE_RES ");
unsetenv((char *)"VPR_HEATING");
unsetenv((char *)"VPR_HEATING_RES");
unsetenv((char *)"LOG_VPR");

return;
}
void printwork(){
LOG_CATEGORY("Environment");
LOG_INFO("FIRST_LEVEL_FILE=%s",getenv((char *)"FIRST_LEVEL_FILE"));
LOG_INFO("FIRST_LEVEL_DIM_FILE=%s",getenv((char *)"FIRST_LEVEL_DIM_FILE"));
LOG_INFO("FILE_T=%s",getenv((char *)"FILE_T"));
LOG_INFO("VPR_HMAX=%s",getenv((char *)"VPR_HMAX"));
LOG_INFO("FILE_ZERO_TERMICO=%s",getenv((char *)"FILE_ZERO_TERMICO"));
LOG_INFO("DIR_STORE_VPR=%s",getenv((char *)"DIR_STORE_VPR"));
LOG_INFO("SITO=%s",getenv((char *)"SITO"));
LOG_INFO("VPR_ARCH=%s",getenv((char *)"VPR_ARCH"));
LOG_INFO("FILE_DEM_SPC=%s",getenv((char *)"FILE_DEM_SPC"));
LOG_INFO("FILE_DEM_GAT=%s",getenv((char *)"FILE_DEM_GAT"));
                
LOG_INFO("TZ=%s",getenv((char *)"TZ"));
LOG_INFO("DIR_OUT_PP_BLOC=%s",getenv((char *)"DIR_OUT_PP_BLOC"));
LOG_INFO("OUTPUT_Z_DIR=%s",getenv((char *)"OUTPUT_Z_DIR"));
LOG_INFO("OUTPUT_RAIN_DIR=%s",getenv((char *)"OUTPUT_RAIN_DIR"));
LOG_INFO("BACINI_DIR=%s",getenv((char *)"BACINI_DIR"));
LOG_INFO("OUTPUT_Z_LOWRIS_DIR=%s",getenv((char *)"OUTPUT_Z_LOWRIS_DIR"));
LOG_INFO("DIR_QUALITY=%s",getenv((char *)"DIR_QUALITY"));
LOG_INFO("LISTA_FILE=%s",getenv((char *)"LISTA_FILE"));
LOG_INFO("LAST_FILE=%s",getenv((char *)"LAST_FILE"));
LOG_INFO("LOG_FILE=%s",getenv((char *)"LOG_FILE"));
LOG_INFO("TEST_VPR=%s",getenv((char *)"TEST_VPR"));
LOG_INFO("ANAP_STAT_FILE=%s",getenv((char *)"ANAP_STAT_FILE"));
LOG_INFO("BLOC_STAT_FILE=%s",getenv((char *)"BLOC_STAT_FILE"));
LOG_INFO("ELEV_STAT_FILE=%s",getenv((char *)"ELEV_STAT_FILE"));
LOG_INFO("NOME_SIDME=%s",getenv((char *)"NOME_SIDME"));
LOG_INFO("MATRICE_BACINI=%s",getenv((char *)"MATRICE_BACINI"));
LOG_INFO("BACINI_AUS_FILE=%s",getenv((char *)"BACINI_AUS_FILE"));
LOG_INFO("BACINI_BOLOGNA=%s",getenv((char *)"BACINI_BOLOGNA"));
LOG_INFO("BACINI_HISTORY_FILE=%s",getenv((char *)"BACINI_HISTORY_FILE"));
LOG_INFO("LAST_VPR=%s",getenv((char *)"LAST_VPR"));
LOG_INFO("LAST_VPR_RES=%s",getenv((char *)"LAST_VPR_RES"));
LOG_INFO("VPR0_FILE=%s",getenv((char *)"VPR0_FILE"));
LOG_INFO("VPR0_FILE_RES=%s",getenv((char *)"VPR0_FILE_RES "));
LOG_INFO("VPR_HEATING=%s",getenv((char *)"VPR_HEATING"));
LOG_INFO("VPR_HEATING_RES=%s",getenv((char *)"VPR_HEATING_RES"));
LOG_INFO("LOG_VPR=%s",getenv((char *)"LOG_VPR"));

}

