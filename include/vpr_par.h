
//PARAMETRI VPR

#define  TCK_VPR 200  // M ALTEZZA STRATO CAMPIONAMENTO
#define  RMIN_VPR 500 // M DISTANZA MINIMA CONSIDERATA
#define  RMAX_VPR 70000 // M DISTANZA MASSIMA CONSIDERATA
#define  VEXTMIN_VPR 2000  // M SPESSORE MINIMO DATI USATI
#define  THR_VPR 13 //DBZ MINIMO PER CONSIDERARE DATO
#define  THR_CORR 0 //DBZ MINIMO PER CORREGGERE DATO
#define  CL_VPR  0  // USO SOLO DATI PRIVI DI CLUTTER
#define  BBMAX_VPR  20  // USO SOLO DATI CON VISIBILITA' > 80%
#define  QMIN_VPR  20  // USO SOLO DATI CON QUALITA' > 20
#define  CT_MIN  0.01    // FRAZIONE DI VOLUME MINIMO PIENO PER CALCOLO PROF.
#define  K_C0  2 // COSTANTE MOLTIPLICATIVA PER RICAVARE C0 DA CV
#define  MIN_AREA   100000 //area minima per strato in km^2/1000 cambiata a 100000
#define  NMAXLAYER  70 // NUMERO MASSIMO STRATI (AGGIUNTA ANNA)
#define  NODATAVPR  -9999. // 
#define  INODATA     -9999
#define  IAZ_MIN_SPC  -100   // spc -90° (-100)--> 135°(150) gat -45°(-50)--> 99°(110)
#define  IAZ_MIN_GAT  -50
#define  IAZ_MAX_SPC   150
#define  IAZ_MAX_GAT   110
#define  MOD_VPR 0 // modalità calcolo VPR: se 0 vpr mediato, se 1 istantaneo
#define  WARM      3 //no profili per 'riscaldamento' VPR
#define  MEMORY      14 // no di istanti di 'raffreddamento'
#define HALF_BB 500 // M spessore indicativo del MELTING LAYER
#define T_MIN_ML 0. //  °C temperatura minima indicativa della bright band
#define T_MAX_ML 4. // °C  temperatura massima indicativa della Bright band
#define THR_STDEVR  0.03  // soglia sulla standard deviation di vol R  a 1100 m normalizzato rispetto vol R totale
#define DCHISQTHR 0.00001 //minima differenza di chisquare tra una iterazione e l' altra
#define DIST_MAX 10.  // (old 5) massima distanza tra i punti e il profilo di prima approssimazione affinchè i punti siano considerati
#define CHISQ_MAX 15. // (old 30) massimo chisquare affinchè l'interpolazione sia accettata
#define MAX_CORR 20. // massima correzione in dB
#define THR_SN 1.36 //differenza massima tra vprmax e vpr a hmax+400m per neve (in rapporto precipitazione
#define T_MAX_SN 2. // t max neve
#define MIN_PEAK_VPR 3.5 // picco minim VPR per cercare massimo in dbR
# define maxdlon 0.1  // massimo scarto lon  per usare temperature per calcolo t radar 
# define maxdlat 0.1 // massimo scarto lat  per usare temperature per calcolo t radar 
