#include <stdio.h> 
#include <stdlib.h> 
#include <errno.h> 
#include <string.h> 

int setstat(char *sito, int mese, char *n_dem, char *n_fl)
{  
  char *fl,*dem;
  int ier;

  /*settaggio files statici: dem e first level*/  

  if (strcmp(sito,"SPC") && strcmp(sito,"GAT")) {
      printf("errore, il sito è sconosciuto");
      exit(1);
    }
  if (!strcmp(sito,"SPC")){    
    if (getenv("FILE_DEM_SPC")==NULL) strcpy(n_dem,"../../PP+BLOC/dati/dem_SanPi.txt");
    else strcpy(n_dem,getenv("FILE_DEM_SPC"));
    if (mese > 0 && mese < 4){ 
      if (getenv("FIRST_LEVEL_FILE")==NULL) strcpy(n_fl,"../dati/FIRST_LEVEL_SPC_2006_INV");
      else strcpy(n_fl,getenv("FIRST_LEVEL_FILE"));
    }
    if (mese > 9 && mese <=12){ 
      if (getenv("FIRST_LEVEL_FILE")==NULL) strcpy(n_fl,"../dati/FIRST_LEVEL_SPC_2006_AUT");
      else strcpy(n_fl,getenv("FIRST_LEVEL_FILE")); 
    }
    if (mese > 3 && mese < 10){  
      if (getenv("FIRST_LEVEL_FILE")==NULL) strcpy(n_fl,"../dati/FIRST_LEVEL_SPC_2006_PRI-EST");
      else strcpy(n_fl,getenv("FIRST_LEVEL_FILE"));
    }
  }
  if (!strcmp(sito,"GAT")){
    if (getenv("FILE_DEM_GAT")==NULL) strcpy(n_dem,"../../PP+BLOC/dati/dem_Gatta.txt");
    else strcpy(n_dem,getenv("FILE_DEM_GAT"));
    if (mese > 0 && mese < 4) 
      if (getenv("FIRST_LEVEL_FILE")==NULL) strcpy(n_fl,"../dati/FIRST_LEVEL_GAT_2006_INV");
      else strcpy(n_fl,getenv("FIRST_LEVEL_FILE"));
    if (mese > 9 && mese <=12) 
      if (getenv("FIRST_LEVEL_FILE")==NULL) strcpy(n_fl,"../dati/FIRST_LEVEL_GAT_2006_AUT");
      else strcpy(n_fl,getenv("FIRST_LEVEL_FILE")); 
    if (mese > 3 && mese < 10) 
      if (getenv("FIRST_LEVEL_FILE")==NULL) strcpy(n_fl,"../dati/FIRST_LEVEL_GAT_2006_PRI-EST");
      else strcpy(n_fl,getenv("FIRST_LEVEL_FILE"));
  }
 
  // if (getenv("FILE_DEM")==NULL) ier=putenv(dem);
  //if (getenv("FIRST_LEVEL_FILE")==NULL) ier=putenv(fl);

  return 0;
}
