#include <setwork.h>

#include <stdlib.h> 
#include <errno.h> 
#include <string.h> 
#include <stdio.h>

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
printf("FIRST_LEVEL_FILE=%s\n",getenv((char *)"FIRST_LEVEL_FILE"));
printf("FIRST_LEVEL_DIM_FILE=%s\n",getenv((char *)"FIRST_LEVEL_DIM_FILE"));
printf("FILE_T=%s\n",getenv((char *)"FILE_T"));
printf("VPR_HMAX=%s\n",getenv((char *)"VPR_HMAX"));
printf("FILE_ZERO_TERMICO=%s\n",getenv((char *)"FILE_ZERO_TERMICO"));
printf("DIR_STORE_VPR=%s\n",getenv((char *)"DIR_STORE_VPR"));
printf("SITO=%s\n",getenv((char *)"SITO"));
printf("VPR_ARCH=%s\n",getenv((char *)"VPR_ARCH"));
printf("FILE_DEM_SPC=%s\n",getenv((char *)"FILE_DEM_SPC"));
printf("FILE_DEM_GAT=%s\n",getenv((char *)"FILE_DEM_GAT"));
                
printf("TZ=%s\n",getenv((char *)"TZ"));
printf("DIR_OUT_PP_BLOC=%s\n",getenv((char *)"DIR_OUT_PP_BLOC"));
printf("OUTPUT_Z_DIR=%s\n",getenv((char *)"OUTPUT_Z_DIR"));
printf("OUTPUT_RAIN_DIR=%s\n",getenv((char *)"OUTPUT_RAIN_DIR"));
printf("BACINI_DIR=%s\n",getenv((char *)"BACINI_DIR"));
printf("OUTPUT_Z_LOWRIS_DIR=%s\n",getenv((char *)"OUTPUT_Z_LOWRIS_DIR"));
printf("DIR_QUALITY=%s\n",getenv((char *)"DIR_QUALITY"));
printf("LISTA_FILE=%s\n",getenv((char *)"LISTA_FILE"));
printf("LAST_FILE=%s\n",getenv((char *)"LAST_FILE"));
printf("LOG_FILE=%s\n",getenv((char *)"LOG_FILE"));
printf("TEST_VPR=%s\n",getenv((char *)"TEST_VPR"));
printf("ANAP_STAT_FILE=%s\n",getenv((char *)"ANAP_STAT_FILE"));
printf("BLOC_STAT_FILE=%s\n",getenv((char *)"BLOC_STAT_FILE"));
printf("ELEV_STAT_FILE=%s\n",getenv((char *)"ELEV_STAT_FILE"));
printf("NOME_SIDME=%s\n",getenv((char *)"NOME_SIDME"));
printf("MATRICE_BACINI=%s\n",getenv((char *)"MATRICE_BACINI"));
printf("BACINI_AUS_FILE=%s\n",getenv((char *)"BACINI_AUS_FILE"));
printf("BACINI_BOLOGNA=%s\n",getenv((char *)"BACINI_BOLOGNA"));
printf("BACINI_HISTORY_FILE=%s\n",getenv((char *)"BACINI_HISTORY_FILE"));
printf("LAST_VPR=%s\n",getenv((char *)"LAST_VPR"));
printf("LAST_VPR_RES=%s\n",getenv((char *)"LAST_VPR_RES"));
printf("VPR0_FILE=%s\n",getenv((char *)"VPR0_FILE"));
printf("VPR0_FILE_RES=%s\n",getenv((char *)"VPR0_FILE_RES "));
printf("VPR_HEATING=%s\n",getenv((char *)"VPR_HEATING"));
printf("VPR_HEATING_RES=%s\n",getenv((char *)"VPR_HEATING_RES"));
printf("LOG_VPR=%s\n",getenv((char *)"LOG_VPR"));

}

