#include "vpr.h"
#include "radarelab/vpr_par.h"

// TODO: should aplitude be computed rather than hardcoded?
#define AMPLITUDE 0.9
#define DTOR  (M_PI/180.)

using namespace std;

namespace radarelab {
namespace algo {

InstantaneousVPR::InstantaneousVPR(const Volume<double>& volume, const Volume<unsigned char>& qual, Volume<unsigned char>& flag_vpr, int az_min, int az_max)
    : volume(volume), qual(qual), flag_vpr(flag_vpr), az_min(az_min), az_max(az_max), dbz(volume)
{
}

/*
idx funzione calcolo VPR istantaneo
calcola il VPR istantaneo secondo il metodo di Germann e Joss (2003)
Per il calcolo si considerano i punti con Z>THR_VPR, qualità>QMIN_VPR, BeamBlocking<20% e clutter free all'interno del volume scelto.
Il punto del vpr corrispondente ad un livello è calcolato come la somma dei volumi di pioggia di ciascuna cella polare presente nell'intervallo di quota, diviso l'area totale. Per ogni livello devo avere un'area minima coperta di pioggia,se non si verifica, ricopio il valore del primo livello che raggiunge questa copertura.
La parte del profilo alta, che sta cioè al di sopra dell'ultimo livello calcolato, viene riempite tirando una retta con coeff. amgolre negativo e costante, pari a -0.03 (valore passibile di discussione)
Il profilo è poi soggetto a quality check e viene rigettato (return(1)) se:
- estensione verticale troppo bassa (devo avere almeno 2 km di profilo da un punto che sta non più alto di 700 m dal suolo)
- cv<CT_MIN*ct, cioè frazione di volume precipitante piccola
- la deviazione standard del volume di R a 1100 m normalizzata rispetto al volume totale precipitante supera soglia prefissata


COSTANTI
THR_VPR:soglia valore Z per calcolo vpr
QMIN_VPR: qualità minima necessaria per entrare nel calcolo vpr
MIN_AREA: area minima a un livello per avere un punto del vpr
IAZ_MIN_SPC,IAZ_MAX_SPC,IAZ_MIN_GAT,I_AZ_MAX_GAT,R_MIN_VPR,R_MAX_VPR:limiti dell'area considerata per il calcolo VPR
TCK_VPR: spessore dei livelli
NMAXLAYER: n max. livelli
VEXTMIN_VPR: estensione verticale minima del profilo
THR_TDEVR: soglia di massima dev. standard di R per considerare buono il profilo
CT_MIN: frazione di volume minimo precipitante per considerare buono il profilo

VARIABILI
int flag_vpr[][][]: flag che indica l'usabilità della cella polare in termini di posizione del bin, beam blocking e clutter
long int *cv,*ct:  volume del settore libero,volume precipitante
float vpr1[]: vpr istantaneo
long int vert_ext,vol_rain: estensione verticale profilo, volume pioggia del singolo bin
long int area_vpr[NMAXLAYER]; area totale usata per calcolo vpr

*/
int InstantaneousVPR::compute(long int *cv, long int *ct, vector<float>& vpr1, vector<long int>& area_vpr)
{
    LOG_CATEGORY("radar.vpr");
    int i,iA,ilay,il,ilast;
    long int dist,dist_plain,vol_rain;
    float area,quota_true_st;
    float grad;

    //----------------inizializzazioni----------------
    long int vert_ext=0;

    *cv=0;
    *ct=0;


    /*SECTION A: cv and ct retrieval and calculation of current area- weighted vpr*/


    //------------riconoscimento sito per definizione limiti azimut---------
    int iaz_min = az_min;
    int iaz_max = az_max;
    LOG_DEBUG(" Iaz_min %d   iaz_max %d",iaz_min,iaz_max);

    for (unsigned l=0; l < volume.size(); ++l)//ciclo elevazioni
    {
        const PolarScan<double>& scan = volume[l];

        for (unsigned k=0; k < scan.beam_size; k++)/*ciclo range*/
        {
            //-------------calcolo distanza-----------
            // TODO: why the casts to int?
            dist=k*(long int)(scan.cell_size)+(int)(scan.cell_size)/2.;

            //-----ciclo settore azimut???
            for (iA=iaz_min; iA<iaz_max; iA++)//ciclo sulle unità di azimut  (0,9°)    ---------
            {
                //----------ottengo il valore dell'indice in unità di azimut (0,9°) ---------
                i = (iA + volume.beam_count) % volume.beam_count;

                //--------calcolo elevazione e quota---------
                quota_true_st = scan.sample_height_real(i, k);

                //--------trovo ilay---------
                ilay=floor(quota_true_st/TCK_VPR);//  in teoria quota indipendente da azimuth  , in realtà no (se c'è vento)

                if (ilay <0 || ilay >= NMAXLAYER) {
                    //fprintf(log_vpr,"ilay %d errore\n",ilay);
                    break;
                }


                vol_rain=0;
                // dist=(long int)(dist*cos((float)(cum_bac.volume_at_elev_preci(i, k).teta_true)*CONV_RAD));

                /* //---------calcolo la distanza proiettata sul piano-------------  */

                const double elevaz = scan.elevations_real(i) * M_PI / 180.;
                dist_plain=(long int)(dist*cos(elevaz));
                if (dist_plain <RMIN_VPR || dist_plain > RMAX_VPR )
                    flag_vpr.scan(l).set(i, k, 0);
//      if (iA == iaz_min) LOG_DEBUG(" k %3d dist %6d  dist_plain %6d quota_true_st %8.2f ilay %3d elevaz %5.2f %f", k, dist, dist_plain,quota_true_st, ilay,scan.elevation,elevaz);

                if (qual.scan(l).get(i, k) < QMIN_VPR) flag_vpr.scan(l).set(i, k, 0);

                //AGGIUNTA PER CLASS
#if 0
                if(cum_bac.do_class){
                    if(conv(i,k)>= CONV_VAL){
                        flag_vpr.scan(l).set(i, k, 0);
                    }
                }
#endif

                // ------per calcolare l'area del pixel lo considero un rettangolo dim bin x ampiezzamediafascio x flag vpr/1000 per evitare problemi di memoria?
                area = scan.cell_size * dist_plain * AMPLITUDE * DTOR * flag_vpr.scan(l).get(i, k)/1000.; // divido per  mille per evitare nr troppo esagerato

                // ------incremento il volume totale di area

                *cv=*cv+(long int)(area);


                //---------------------condizione per incrementare VPR contributo: valore sopra 13dbz, qualità sopra 20 flag>0 (no clutter e dentro settore)------------------
                double sample = scan.get(i, k);
                if (sample > THR_VPR &&  flag_vpr.scan(l).get(i, k) > 0 )
                {
                    //-------incremento il volume di pioggia = pioggia x area
                    vol_rain=(long int)(dbz.DBZ_to_mp_func(sample)*area);//peso ogni cella con la sua area

                    //-------incremento l'area precipitante totale ct,aggiungendo però,cosa che avevo messo male una THR solo per ct, cioè per il peso
                    if (sample > THR_PDF)
                        *ct=*ct+(long int)(area);

                    //------se l'area in quello strato è già maggiore di 0 allora incremento il volume dello strato altrimenti lo scrivo ex novo. poi vpr1 andrà diviso per l'area
                    if (area_vpr[ilay]> 0) vpr1[ilay]=vpr1[ilay]+(float)(vol_rain);
                    else vpr1[ilay]=(float)(vol_rain);

                    //------incremento l'area dello strato----------
                    area_vpr[ilay]=area_vpr[ilay]+area;
                }
            }
        }
    }

    LOG_INFO("calcolati ct e cv ct= %li cv= %li",*ct,*cv);

    /*SECTION B: vpr quality checks and re-normalisation of vpr*/

    //--------------CONTROLLO DI QUALITA' E NORMALIZZAZIONE DEL PROFILO ISTANTANEO CALCOLATO
    //-------- se il volume supera quello minimo------
    if ((*ct) > CT_MIN*(*cv)) {

        ilast=0;
        vert_ext=0;


        //----- calcolo 'estensione verticale del profilo , negli strati dove l'area è troppo piccola assegno NODATAVPR,  NORMALIZZO il profilo, e se l'estensione è minore di VEXTMIN_VPR esco-


        for (ilay=0; ilay<NMAXLAYER; ilay++){


            LOG_INFO("  ilay %d area_vpr= %ld  ct= %ld  cv= %ld", ilay, area_vpr[ilay],*ct,*cv );

            if (area_vpr[ilay]>=MIN_AREA) {
                vert_ext=vert_ext+TCK_VPR;
                vpr1[ilay]=vpr1[ilay]/(float)(area_vpr[ilay]);

            }
            else
            {
                vpr1[ilay]=NODATAVPR;


                //----  se incontro un punto vuoto oltre 700 m ( o se sono arrivata alla fine) assegno ilast ed esco dal ciclo

                /*   if (ilast > 0 && vert_ext>VEXTMIN_VPR){ */

                if (ilay*TCK_VPR+TCK_VPR/2 > MINPOINTBASE || ilay== NMAXLAYER -1){  //cambio il criterio per calcolare estensione minima. devo avere almeno 2 km consecutivi a partire da un punto che stia almeno a 700 m (MINPOINTBASE),
                    LOG_INFO("raggiunta cima profilo");
                    ilast=ilay-1;// c'era errore!!!

                    //---------- raggiunta la cima profilo faccio check immediato sull'estensione verticale
                    if (vert_ext<VEXTMIN_VPR ){
                        LOG_INFO("estensione profilo verticale troppo bassa");
                        *ct=0;
                        ilast=0;
                        for  (il=0; il<ilast; il++) vpr1[il]=NODATAVPR;
                        return(1);
                    }

                    break; // esco dal ciclo..modifica
                }
            }
        }
    }// fine se volumeprecipitante sufficiente

    // ---------se il volume non supera quello minimo esco---------
    else {
        LOG_INFO("volume precipitante troppo piccolo");
        *ct=0;
        ilast=0;
        for  (il=0; il<NMAXLAYER; il++) vpr1[il]=NODATAVPR; //--devo riassegnare o mi rimane 'sporco' forse si potrebbe usare una ver diversa
        return(1);
    }


    //------calcolo il gradiente del profilo come media del gradiente negli ultimi 3 strati per assegnare la parte 'alta' (novità)

    grad=((vpr1[ilast]-vpr1[ilast-1]) + (vpr1[ilast-1]-vpr1[ilast-2])/(2.)+ (vpr1[ilast-2]-vpr1[ilast-3])/(3.) ) /3.;
    if (grad > 0.002)
        grad=-0.03 ;
    LOG_INFO(" %f", grad);

    //--riempio la parte alta del profilo decrementando di grad il profilo in ogni strato fino a raggiunere 0, SI PUÒ TOGLIERE E METTERE NODATA

    for (ilay=ilast+1; ilay<NMAXLAYER; ilay++) {
        if (vpr1[ilay-1] + grad > 0.002)
            vpr1[ilay]= vpr1[ilay-1]+grad;
        else
            vpr1[ilay]=0;
    }

    // HO CAMBIATO DA GRADIENTE FISSO PARI A (V. VECCHIO) A GRADIENTE RICAVATO DAL PROFILO PER LA PARTE ALTA
    //---HO TOLTO TUTTA LA PARTE CHE FA IL CHECK SULLA STDEV A 1100 M SI PUO' RIVEDERE SE METTERLA SE SERVE.


    return(0);
}


}
}
