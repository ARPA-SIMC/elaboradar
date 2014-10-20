#include "algo/anaprop.h"
#include "algo/utils.h"

// Soglie algoritmi
#define MAX_DIF_OR 30.            /* differenzio limiti controllo anap      */
#define MIN_VALUE_OR -10.         /* a seconda che sia alla prima o success.*/
#define MAX_DIF_NEXT_OR 15.       /*/ elevazione                             */
#define MIN_VALUE_NEXT_OR 0.
#define THR_CONT_ANAP 1 /* limite in numero occorrenze anaprop sul raggio dopo i 30 km per non togliere =*/

namespace elaboradar {
namespace algo {

using namespace std;

namespace anaprop {

GridStats::GridStats()
{
}

GridStats::~GridStats()
{
    if (stat_anap) delete[] stat_anap;
    if (stat_tot) delete[] stat_tot;
    if (stat_bloc) delete[] stat_bloc;
    if (stat_elev) delete[] stat_elev;
}

void GridStats::init(const Volume<double>& volume)
{
    size_az = volume[0].beam_count / step_stat_az + 1;
    size_beam = volume[0].beam_size / step_stat_range + 1;

    stat_anap = new unsigned[size_az * size_beam];
    stat_tot = new unsigned[size_az * size_beam];
    stat_bloc = new unsigned[size_az * size_beam];
    stat_elev = new unsigned[size_az * size_beam];

    for (unsigned i = 0; i < size_az * size_beam; ++i)
        stat_anap[i] = stat_tot[i] = stat_bloc[i] = stat_elev[i] = 0;
}

}

template<class T>
Anaprop<T>::Anaprop()
{
    logging_category = log4c_category_get("radar.anaprop");
}

template<class T>
void Anaprop<T>::init(const Volume<T>& volume)
{
    elev_fin.init(volume);
    grid_stats.init(volume);

    dato_corrotto.resize(400, volume.max_beam_size());
    dato_corrotto.fill(0);
    quota.resize(400, volume.max_beam_size());
    quota.fill(0);
}

template<class T>
void Anaprop<T>::init_elev_fin_static(const Volume<T>& volume, const PolarScan<unsigned char>& first_level_static)
{
    init(volume);
    for(unsigned i=0; i < volume[0].beam_count; ++i)
        for(unsigned k=0; k < volume[0].beam_size; ++k)
            //---assegno el_inf a mappa statica
            elev_fin[i][k] = first_level_static(i, k);
}

template<class T>
void Anaprop<T>::remove(
        Volume<T>& volume,
        PolarScan<unsigned char>& beam_blocking,
        const PolarScan<unsigned char>& first_level,
        const PolarScan<unsigned char>& first_level_static,
        const Volume<double>& SD)
{
    const double fondo_scala = volume[0].undetect;

    init(volume);

    //for(unsigned i=200; i<201; i++)
    for(unsigned i=0; i<NUM_AZ_X_PPI; i++) 
    {
        //------------assegno le soglie per anaprop : se sono oltre 60 km e se la differenza tra il bin sotto il base e quello sopra <10 non applico test (cambio i limiti per renderli inefficaci)
        /* differenza massima tra le due elevazioni successive perchè non sia clutter e valore minimo a quella superiore pe il primo e per i successivi (NEXT) bins*/

        bool do_test_AP = true;
        double MAX_DIF=MAX_DIF_OR;
        double MAX_DIF_NEXT=MAX_DIF_NEXT_OR;
        double MIN_VALUE=MIN_VALUE_OR;
        double MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;

        double MIN_VALUE_USED, MAX_DIF_USED;

        bool flag_anap = false;
        unsigned cont_anap=0;// aggiunto per risolvere problema di uso con preci shallow
        unsigned count_first_elev=0;

        //for(unsigned k=150 ; k<volume[0].beam_size; k++)
        for(unsigned k=0; k<volume[0].beam_size; k++)
        {
            //------------- incremento statistica tot ------------------
            grid_stats.incr_tot(i, k);
            // ------------assegno l'elevazione el_inf a first_level e elev_fin a el_inf---------
            int loc_el_inf = first_level(i, k);
            while ( k >= volume[loc_el_inf].beam_size)
            {
                LOG_INFO("Decremento el_inf per k fuori range (i,k,beam_size,el_inf_dec) (%d,%d,%d,%d)",i,k,volume[loc_el_inf].beam_size,loc_el_inf-1);
                loc_el_inf--;
                if (loc_el_inf < 0) throw std::runtime_error("loc_el_inf < 0");
            }
            while (loc_el_inf > 0 && SD[loc_el_inf-1].get(i,k) < conf_texture_threshold &&  SD[loc_el_inf-1].get(i,k)>=0.01 && volume[loc_el_inf-1].get(i,k) > volume[loc_el_inf].get(i,k)){
                //  LOG_WARN("Decremento el_inf Sotto esiste qualcosa %2d %3d %3d %6.2f %6.2f %6.2f",loc_el_inf, i, k , SD[loc_el_inf-1].get(i,k),volume[loc_el_inf-1].get(i,k),volume[loc_el_inf].get(i,k));
                loc_el_inf--;
            }
            const unsigned el_inf = loc_el_inf;
            if (el_inf == 0) count_first_elev++;
            if (do_quality)
                elev_fin[i][k]=el_inf;

            // ------------assegno a el_up il successivo di el_inf e se >=NEL metto bin_high=fondo_scala
            const unsigned el_up = el_inf +1;

            // ------------assegno bin_low e bin_high 

            float bin_low  = volume[el_inf].get(i, k);
            float bin_high;
            if (el_up >= volume.size() || k >= volume[el_up].beam_size){
                bin_high=fondo_scala;
            } else{
                bin_high = volume[el_up].get(i, k);
            }

            //----------questo serviva per evitare di tagliare la precipitazione shallow ma si dovrebbe trovare un metodo migliore p.es. v. prove su soglia
            if(bin_high == fondo_scala && SD[el_inf].get(i,k)<= conf_texture_threshold && SD[el_inf].get(i,k) > 0.01)                     //-----------ANNULLO EFFETTO TEST ANAP
            {
		do_test_AP=false;
                MAX_DIF_NEXT=BYTEtoDB(255);
                MAX_DIF=BYTEtoDB(255);
                MIN_VALUE=fondo_scala;
                MIN_VALUE_NEXT= fondo_scala;
            }
            else
            {
                do_test_AP=true;
                MAX_DIF=MAX_DIF_OR;
                MAX_DIF_NEXT=MAX_DIF_NEXT_OR;
                MIN_VALUE=MIN_VALUE_OR;
                MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;
            }
            // ------------separo i diversi casi x analisi anaprop: ho dati sia al livello base che sopra o no  e ho trovato anaprop in precedenza sul raggio o no
            bool test_an;
            if (cont_anap> THR_CONT_ANAP ||count_first_elev < 80  )
                test_an=(bin_low > fondo_scala && bin_high >= fondo_scala );
            else
                test_an=(bin_low > fondo_scala && bin_high > fondo_scala );
            //  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f - cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --  %6.2f %6.2f %6.2f",i,k,el_inf,el_up,bin_low,bin_high,   cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT, SD[el_inf].get(i,k),SD[el_inf].get((i+1)%NUM_AZ_X_PPI,k) ,SD[el_inf].get((i-1+NUM_AZ_X_PPI)%NUM_AZ_X_PPI,k));
            //------------------ se ho qualcosa sia al livello base che sopra allora effettuo il confronto-----------------
            if(test_an )
            {
                //------------------ se ho trovato anap prima nel raggio cambio le soglie le abbasso)-----------------
                if(flag_anap)
                {
                    //-----------caso di propagazione anomala presente nella cella precedente ---------
                    MAX_DIF_USED   = MAX_DIF_NEXT;
                    MIN_VALUE_USED = MIN_VALUE_NEXT;
                }
                else
                {
                    //-----------caso di propagazione anomala non presente nella cella precedente ---------
                    MAX_DIF_USED   = MAX_DIF;
                    MIN_VALUE_USED = MIN_VALUE;
                }

                if( do_test_AP &&
                   ( 
                    (
                      bin_low-bin_high >= MAX_DIF_USED 
                           || 
                      (
                          bin_high <= MIN_VALUE_USED 
                                  && 
                          bin_low > MIN_VALUE + 5 
                      )
                    ) 
                                  || 
                    (
                      SD[el_inf].get(i,k) > conf_texture_threshold
//                                  &&  
//                    (
//                          SD[el_inf].get((i+1)%NUM_AZ_X_PPI,k) < 3
//                                 || 
//                          SD[el_inf].get((i-1+NUM_AZ_X_PPI)%NUM_AZ_X_PPI,k) < 3 
//                    )
                    )
                   )
                   )
               {
                   //--------ricopio valore a el_up su tutte elev inferiori--------------
                   if(k < volume[el_up].beam_size && SD[el_up].get(i,k) <= conf_texture_threshold  && SD[el_up].get(i,k) >= 0.01){
                       for(unsigned l=0; l<el_up; l++){
                           volume[l].set(i, k, bin_high);  //ALTERN
                       }
                       //  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --  %6.2f %6.2f %6.2f %6.2f TA-AN",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT, SD[el_inf].get(i,k),SD[el_inf].get((i+1)%NUM_AZ_X_PPI,k) ,SD[el_inf].get((i-1+NUM_AZ_X_PPI)%NUM_AZ_X_PPI,k), SD[el_up].get(i,k)) ;
                   } else {
                       for(unsigned l=0; l<el_up; l++){
                           volume[l].set(i, k, fondo_scala);  //ALTERN
                       }
                       // if (k < volume[el_up].beam_size) LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --   %6.2f %6.2f %6.2f TA-AN set to fondo_scala",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT,SD[el_up].get(i,k),SD[el_up].get((i+1)%NUM_AZ_X_PPI,k) ,SD[el_up].get((i-1+NUM_AZ_X_PPI)%NUM_AZ_X_PPI,k));
                   }
                   //---------assegno l'indicatore di presenza anap nel raggio e incremento statistica anaprop, assegno matrici che memorizzano anaprop e elevazione_finale e azzero beam blocking perchè ho cambiato elevazione
                   flag_anap = true;
                   cont_anap=cont_anap+1;
                   grid_stats.incr_anap(i, k);
                   if (do_quality)
                   {
                       dato_corrotto(i, k)=ANAP_YES;/*matrice risultato test: propagazione anomala*/
                       elev_fin[i][k]=el_up;
                   }
                   if (el_up > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                   if (do_beamblocking)
                       beam_blocking(i, k)=0;
               }

              //-----non c'è propagazione anomala:ricopio su tutte e elevazioni il valore di el_inf e correggo il beam blocking,  incremento la statistica beam_blocking, assegno matrice anaprop a 0 nel punto , assegno a 0 indicatore anap nel raggio, assegno elevazione finale e incremento statisica cambio elevazione se el_inf > first_level_static(i, k)-----------
               else
               {
                   unsigned count_low =0;
                   unsigned count_high=0;
                   for (unsigned ii=0; ii<7; ii++){
                       int iaz_prox= (i+ii-3+NUM_AZ_X_PPI)%NUM_AZ_X_PPI;
                       if( SD[el_inf].get(iaz_prox,k) < conf_texture_threshold && SD[el_inf].get(iaz_prox,k) > 0.01) count_low++;
                       if( k < SD[el_up].beam_size && SD[el_up].get(iaz_prox,k) < 1.7*conf_texture_threshold&& SD[el_up].get(iaz_prox,k) > 0.01) count_high++;
                   }
                   if ( !(SD[el_inf].get(i,k) < 1.3*conf_texture_threshold && SD[el_inf].get(i,k) > 0.01 && count_low >=5)) {
                       if ( k >= SD[el_up].beam_size || !(SD[el_up].get(i,k) < 1.7* conf_texture_threshold && SD[el_up].get(i,k) > 0.01 && count_high >=3)){
                           bin_low = fondo_scala;
                           //		        if ( k < SD[el_up].beam_size )  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --  %6.2f %1d TA-NO_AN low=fondo",i,k,el_inf,el_up,bin_low,bin_high,  volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT, SD[el_up].get(i,k),count_high );  
                           //		        else   LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  -- %1d TA-NO_AN low=fondo",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT,count_high );
                       }
                       else {
                           bin_low=bin_high;
                           if (do_quality)
                           {
                               elev_fin[i][k]=el_up;
                           }
                           if (el_up > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                           if (do_beamblocking)
                               beam_blocking(i, k)=0;
                       }
                   }
                   else {	
                       if (do_beamblocking && do_bloccorr)
                       {
                           bin_low = beam_blocking_correction(bin_low, beam_blocking(i, k));
                           grid_stats.incr_bloc(i, k, beam_blocking(i, k));
                       }
		  }
                  for(unsigned l=0; l<=el_inf; l++)
                      volume[l].set(i, k, bin_low);
//  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --  %6.2f %1d TA-NO_AN",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT, SD[el_inf].get(i,k),count_low );

                  if (do_quality)
                  {
                      dato_corrotto(i, k)=ANAP_OK;
                      elev_fin[i][k]=el_inf;
                  }
                  if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                  flag_anap = false;
                }
            }/* test_anap */
            //----------------se al livello base non ho dato riempio con i valori di el_up tutte le elevazioni sotto (ricostruisco il volume) e assegno beam_blocking 0
            else if (bin_low < fondo_scala)
            {
                for(unsigned l=0; l<el_up; l++)
                {
#warning Here the .get(i, k) on level el_up may find a smaller beam size than the one k can reach, causing a read out of bound
                    if (volume[l].beam_size > k && volume[el_up].beam_size > k)
                        volume[l].set(i, k, bin_high);
                    else if (volume[l].beam_size > k )
                        volume[l].set(i, k, fondo_scala);
		
                }
//  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --  %6.2f %6.2f %6.2f NO_TA-low <fondo",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT, SD[el_inf].get(i,k),SD[el_inf].get((i+1)%NUM_AZ_X_PPI,k) ,SD[el_inf].get((i-1+NUM_AZ_X_PPI)%NUM_AZ_X_PPI,k));
                //----------------controlli su bin_high nel caso in cui bin_low sia un no data per assegnare matrice anap  (dato_corrotto(i, k))
                if (do_quality)
                {
                    if (bin_high<fondo_scala)   dato_corrotto(i, k)=ANAP_NODAT;/*manca dato sotto e sopra*/
                    bool test_an1;
                    if (cont_anap< THR_CONT_ANAP )
                        test_an1=(bin_high>=fondo_scala); //modificato per contemplare > o >=
                    else
                        test_an1=(bin_high>fondo_scala);

                    if (test_an1) dato_corrotto(i, k)=ANAP_NOCONTROL;
                    if (bin_high==fondo_scala) dato_corrotto(i, k)=ANAP_OK;/*non piove (oppure sono sopra livello preci...)*/
                }

                if (do_beamblocking)
                    beam_blocking(i, k)=0;
            }

            //----------------se bin_low == fondo_scala riempio matrice volume.vol_pol con dato a el_inf (mi resta il dubbio di quest'if se seve o basti un else ) azzero matrice anap (dato ok)
            else if (bin_low == fondo_scala || bin_high <= fondo_scala)/* quel che resta da (bin_low > fondo_scala && bin_high >= fondo_scala) e (bin_low < fondo_scala) ; messo =per bin_high*/

            {
                unsigned count =0;
                for (unsigned ii=0; ii<7; ii++){
                    int iaz=(i+ii-3+NUM_AZ_X_PPI)%NUM_AZ_X_PPI;
                    if( SD[el_inf].get(iaz,k) < conf_texture_threshold && SD[el_inf].get(iaz,k) > 0.01) count++;
                }
                if ( !(SD[el_inf].get(i,k) < conf_texture_threshold && SD[el_inf].get(i,k) >0.01 && count >=5 )) 
                    bin_low = fondo_scala;

                for(unsigned l=0; l<=el_inf; l++)//riempio con i valori di el_inf tutte le elevazioni sotto (ricostruisco il volume)
                {
                    if (volume[l].beam_size > k)
                        volume[l].set(i, k, bin_low);
                }
//  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  --  %6.2f %6.2f %6.2f NO_TA-low ==fondo",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT, SD[el_inf].get(i,k),SD[el_inf].get((i+1)%NUM_AZ_X_PPI,k) ,SD[el_inf].get((i-1+NUM_AZ_X_PPI)%NUM_AZ_X_PPI,k));

                if (do_quality)
                {
                    dato_corrotto(i, k)=ANAP_OK; // dubbio
                    elev_fin[i][k]=el_inf;
                }

                if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k);
            }
            /*-----------------------------------------------------------fine di tutti gli if-----------*/
            //-----finiti tutti i controlli assegno le varibili di qualita definitive: elevazione, quota calcolata sull'elevazione reale con propagazione standard , e quota relativa al suolo calcolata con elevazione nominale e propagazione da radiosondaggio.

            if (do_quality)
                quota(i, k)=(unsigned short)PolarScanBase::sample_height(
                        elev_fin.elevation_rad_at_elev_preci(i, k),
                        (k + 0.5) * volume[0].cell_size);
        }
     }    //  end for over beam_count
}

template<class T>
void Anaprop<T>::remove(
        Volume<T>& volume,
        PolarScan<unsigned char>& beam_blocking,
        const PolarScan<unsigned char>& first_level,
        const PolarScan<unsigned char>& first_level_static)
{
    const double fondo_scala = volume[0].undetect;

    init(volume);

    //for(unsigned i=200; i<201; i++)
    for(unsigned i=0; i<volume[0].beam_count; i++) 
    {
        //------------assegno le soglie per anaprop : se sono oltre 60 km e se la differenza tra il bin sotto il base e quello sopra <10 non applico test (cambio i limiti per renderli inefficaci)
        /* differenza massima tra le due elevazioni successive perchè non sia clutter e valore minimo a quella superiore pe il primo e per i successivi (NEXT) bins*/

        bool do_test_AP = true;
        double MAX_DIF=MAX_DIF_OR;
        double MAX_DIF_NEXT=MAX_DIF_NEXT_OR;
        double MIN_VALUE=MIN_VALUE_OR;
        double MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;

        double MIN_VALUE_USED, MAX_DIF_USED;

        bool flag_anap = false;
        unsigned cont_anap=0;// aggiunto per risolvere problema di uso con preci shallow
        unsigned count_first_elev=0;

        //for(unsigned k=150 ; k<volume[0].beam_size; k++)
        for(unsigned k=0; k<volume[0].beam_size; k++)
        {
            //------------- incremento statistica tot ------------------
            grid_stats.incr_tot(i, k);
            // ------------assegno l'elevazione el_inf a first_level e elev_fin a el_inf---------
            int loc_el_inf = first_level(i, k);
            while ( k >= volume[loc_el_inf].beam_size)
            {
                LOG_INFO("Decremento el_inf per k fuori range (i,k,beam_size,el_inf_dec) (%d,%d,%d,%d)",i,k,volume[loc_el_inf].beam_size,loc_el_inf-1);
                loc_el_inf--;
                if (loc_el_inf < 0) throw std::runtime_error("loc_el_inf < 0");
            }
/* ---------------------------------
 * PER IL MOMENTO NON BUTTO ANCORA QUESTO CODICE CI DEVO PENSARE
 * while (loc_el_inf > 0 && SD_Z6[loc_el_inf-1].get(i,k) < conf_texture_threshold &&  SD_Z6[loc_el_inf-1].get(i,k)>=0.01 && volume[loc_el_inf-1].get(i,k) > volume[loc_el_inf].get(i,k)){
                //  LOG_WARN("Decremento el_inf Sotto esiste qualcosa %2d %3d %3d %6.2f %6.2f %6.2f",loc_el_inf, i, k , SD_Z6[loc_el_inf-1].get(i,k),volume[loc_el_inf-1].get(i,k),volume[loc_el_inf].get(i,k));
                loc_el_inf--;
            }
*/   
            const unsigned el_inf = loc_el_inf;
            if (el_inf == 0) count_first_elev++;
            if (do_quality)
                elev_fin[i][k]=el_inf;

            // ------------assegno a el_up il successivo di el_inf e se >=NEL metto bin_high=fondo_scala
            const unsigned el_up = el_inf +1;

            // ------------assegno bin_low e bin_high 

            float bin_low  = volume[el_inf].get(i, k);
            float bin_high;
            if (el_up >= volume.size() || k >= volume[el_up].beam_size){
                bin_high=fondo_scala;
            } else{
                bin_high = volume[el_up].get(i, k);
            }

            //----------questo serviva per evitare di tagliare la precipitazione shallow ma si dovrebbe trovare un metodo migliore p.es. v. prove su soglia
            if(bin_high == fondo_scala )                     //-----------ANNULLO EFFETTO TEST ANAP
            {
		do_test_AP=false;
                MAX_DIF_NEXT=BYTEtoDB(255);
                MAX_DIF=BYTEtoDB(255);
                MIN_VALUE=fondo_scala;
                MIN_VALUE_NEXT= fondo_scala;
            }
            else
            {
                do_test_AP=true;
                MAX_DIF=MAX_DIF_OR;
                MAX_DIF_NEXT=MAX_DIF_NEXT_OR;
                MIN_VALUE=MIN_VALUE_OR;
                MIN_VALUE_NEXT=MIN_VALUE_NEXT_OR;
            }
            // ------------separo i diversi casi x analisi anaprop: ho dati sia al livello base che sopra o no  e ho trovato anaprop in precedenza sul raggio o no
            bool test_an;
            if (cont_anap> THR_CONT_ANAP ||count_first_elev < 80  )
                test_an=(bin_low > fondo_scala && bin_high >= fondo_scala );
            else
                test_an=(bin_low > fondo_scala && bin_high > fondo_scala );
            //  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f - cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  ",i,k,el_inf,el_up,bin_low,bin_high,   cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT);
            //------------------ se ho qualcosa sia al livello base che sopra allora effettuo il confronto-----------------
            if(test_an )
            {
                //------------------ se ho trovato anap prima nel raggio cambio le soglie le abbasso)-----------------
                if(flag_anap)
                {
                    //-----------caso di propagazione anomala presente nella cella precedente ---------
                    MAX_DIF_USED   = MAX_DIF_NEXT;
                    MIN_VALUE_USED = MIN_VALUE_NEXT;
                }
                else
                {
                    //-----------caso di propagazione anomala non presente nella cella precedente ---------
                    MAX_DIF_USED   = MAX_DIF;
                    MIN_VALUE_USED = MIN_VALUE;
                }

                if( do_test_AP &&
                      bin_low-bin_high >= MAX_DIF_USED 
                           || 
                      (
                          bin_high <= MIN_VALUE_USED 
                                  && 
                          bin_low > MIN_VALUE + 5 
                      )
                  )
               {
                   //--------ricopio valore a el_up su tutte elev inferiori--------------
                   if(k < volume[el_up].beam_size ){
                       for(unsigned l=0; l<el_up; l++){
                           volume[l].set(i, k, bin_high);  //ALTERN
                       }
                       //  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  TA-AN",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT ) ;
                   } else {
                       for(unsigned l=0; l<el_up; l++){
                           volume[l].set(i, k, fondo_scala);  //ALTERN
                       }
                       // if (k < volume[el_up].beam_size) LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f   TA-AN set to fondo_scala",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT);
                   }
                   //---------assegno l'indicatore di presenza anap nel raggio e incremento statistica anaprop, assegno matrici che memorizzano anaprop e elevazione_finale e azzero beam blocking perchè ho cambiato elevazione
                   flag_anap = true;
                   cont_anap=cont_anap+1;
                   grid_stats.incr_anap(i, k);
                   if (do_quality)
                   {
                       dato_corrotto(i, k)=ANAP_YES;/*matrice risultato test: propagazione anomala*/
                       elev_fin[i][k]=el_up;
                   }
                   if (el_up > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                   if (do_beamblocking)
                       beam_blocking(i, k)=0;
               }

              //-----non c'è propagazione anomala:ricopio su tutte e elevazioni il valore di el_inf e correggo il beam blocking,  incremento la statistica beam_blocking, assegno matrice anaprop a 0 nel punto , assegno a 0 indicatore anap nel raggio, assegno elevazione finale e incremento statisica cambio elevazione se el_inf > first_level_static(i, k)-----------
               else
               {
                       if (do_beamblocking && do_bloccorr)
                       {
                           bin_low = beam_blocking_correction(bin_low, beam_blocking(i, k));
                           grid_stats.incr_bloc(i, k, beam_blocking(i, k));
                       }
		  
                  for(unsigned l=0; l<=el_inf; l++)
                      volume[l].set(i, k, bin_low);
//  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f   TA-NO_AN",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT );

                  if (do_quality)
                  {
                      dato_corrotto(i, k)=ANAP_OK;
                      elev_fin[i][k]=el_inf;
                  }
                  if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k);//incremento la statistica cambio elevazione
                  flag_anap = false;
                }
            }/* test_anap */
            //----------------se al livello base non ho dato riempio con i valori di el_up tutte le elevazioni sotto (ricostruisco il volume) e assegno beam_blocking 0
            else if (bin_low < fondo_scala)
            {
                for(unsigned l=0; l<el_up; l++)
                {
#warning Here the .get(i, k) on level el_up may find a smaller beam size than the one k can reach, causing a read out of bound
                    if (volume[l].beam_size > k && volume[el_up].beam_size > k)
                        volume[l].set(i, k, bin_high);
                    else if (volume[l].beam_size > k )
                        volume[l].set(i, k, fondo_scala);
		
                }
//  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  NO_TA-low <fondo",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT);
                //----------------controlli su bin_high nel caso in cui bin_low sia un no data per assegnare matrice anap  (dato_corrotto(i, k))
                if (do_quality)
                {
                    if (bin_high<fondo_scala)   dato_corrotto(i, k)=ANAP_NODAT;/*manca dato sotto e sopra*/
                    bool test_an1;
                    if (cont_anap< THR_CONT_ANAP )
                        test_an1=(bin_high>=fondo_scala); //modificato per contemplare > o >=
                    else
                        test_an1=(bin_high>fondo_scala);

                    if (test_an1) dato_corrotto(i, k)=ANAP_NOCONTROL;
                    if (bin_high==fondo_scala) dato_corrotto(i, k)=ANAP_OK;/*non piove (oppure sono sopra livello preci...)*/
                }

                if (do_beamblocking)
                    beam_blocking(i, k)=0;
            }

            //----------------se bin_low == fondo_scala riempio matrice volume.vol_pol con dato a el_inf (mi resta il dubbio di quest'if se seve o basti un else ) azzero matrice anap (dato ok)
            else if (bin_low == fondo_scala || bin_high <= fondo_scala)/* quel che resta da (bin_low > fondo_scala && bin_high >= fondo_scala) e (bin_low < fondo_scala) ; messo =per bin_high*/

            {
                unsigned count =0;
                    bin_low = fondo_scala;

                for(unsigned l=0; l<=el_inf; l++)//riempio con i valori di el_inf tutte le elevazioni sotto (ricostruisco il volume)
                {
                    if (volume[l].beam_size > k)
                        volume[l].set(i, k, bin_low);
                }
//  LOG_WARN("b@(%3d,%3d) - el_inf %2d  - el_up %2d -low %6.2f - up %6.2f fin %6.2f- cont %3d %1d %1d %6.2f %6.2f %6.2f %6.2f  NO_TA-low ==fondo",i,k,el_inf,el_up,bin_low,bin_high, volume[0].get(i,k),cont_anap,test_an, flag_anap, MAX_DIF, MIN_VALUE, MAX_DIF_NEXT, MIN_VALUE_NEXT);

                if (do_quality)
                {
                    dato_corrotto(i, k)=ANAP_OK; // dubbio
                    elev_fin[i][k]=el_inf;
                }

                if (el_inf > first_level_static(i, k)) grid_stats.incr_elev(i, k);
            }
            /*-----------------------------------------------------------fine di tutti gli if-----------*/
            //-----finiti tutti i controlli assegno le varibili di qualita definitive: elevazione, quota calcolata sull'elevazione reale con propagazione standard , e quota relativa al suolo calcolata con elevazione nominale e propagazione da radiosondaggio.

            if (do_quality)
                quota(i, k)=(unsigned short)PolarScanBase::sample_height(
                        elev_fin.elevation_rad_at_elev_preci(i, k),
                        (k + 0.5) * volume[0].cell_size);
        }
     }    //  end for over beam_count
}

// Explicit template instantiation for T=double
template class Anaprop<double>;

}
}
