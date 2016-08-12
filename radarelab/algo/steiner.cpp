#include "steiner.h"
#include "radarelab/par_class.h"
#include "radarelab/algo/dbz.h"

//#ifdef __cplusplus
//extern "C" {
//#endif
//#include <func_Z_R.h>
//#ifdef __cplusplus
//}
//#endif

// valore mancante
#define MISSING 0

// fattore conversione gradi-radianti
#define DTOR  M_PI/180.

//Definizioni geometriche
static const double AMPLITUDE = 0.9; /* esternalizzo?*/ // ampiezza fascio radar

using namespace std;

namespace radarelab {
namespace algo {

namespace steiner {

void Point::add_sample(double sample)
{
    if (sample <= MINVAL_DB) return;
    Z_bckgr += DBZ::BYTEtoZ(DBZ::DBtoBYTE(sample));
    bckgr += sample;
    ++npoints;
}

void Point::finalize()
{
    if (npoints > 0){
        Z_bckgr = Z_bckgr / npoints;
        //bckgr[i]=bckgr[i]/npoints; //no
        if (Z_bckgr > 0) bckgr = 10 * (log10(Z_bckgr));
    }
    //il valore del raggio convettivo varia a seconda del background, da 1 a 5 km
    if (bckgr < 25.)
        convective_radius = 1.;
    else if (bckgr >= 25. && bckgr < 30.)
        convective_radius = 2.;
    else if (bckgr >= 30. && bckgr < 35.)
        convective_radius = 3.;
    else if (bckgr >= 35. && bckgr < 40.)
        convective_radius = 4.;
    else if (bckgr > 40.)
        convective_radius = 5.;
}

}

CalcoloSteiner::CalcoloSteiner(
        const Volume<double>& volume,
        const volume::ElevFin<double>& elev_fin,
        unsigned max_bin)
    // TODO: evaluate cell_size scan by scan?
    : volume(volume), elev_fin(elev_fin), max_bin(max_bin), size_cell(volume.scan(0).cell_size),
      conv_STEINER(Matrix2D<unsigned char>::Constant(NUM_AZ_X_PPI, max_bin, MISSING))
{
    using namespace steiner;

    logging_category = log4c_category_get("radar.vpr");

    // Traccio una lista dei punti che hanno valore non nullo e sotto base
    // bright band (lista_bckg) contenente iaz e irange e conto i punti
    for (unsigned i=0; i < NUM_AZ_X_PPI; ++i)
        for (unsigned j=0; j < max_bin; ++j)  // propongo max_bin visto che risoluzione è la stessa
            //if ( volume.scan(0)[i][j] > 1 &&  (float)(quota[i][j])/1000. < hbbb ) //verifico che il dato usato per la ZLR cioè la Z al lowest level sia > soglia e la sua quota sia sotto bright band o sopra bright band
            if (j < volume.scan(0).beam_size && volume.scan(0).get(i, j) > MINVAL_DB)
                lista_bckg.push_back(Point(i, j)); // IAZIMUT, IRANGE
}

void CalcoloSteiner::calcolo_background() // sui punti precipitanti calcolo bckgr . nb LA CLASSIFICAZIONE DI STEINER NON HA BISOGNO DI RICAMPIONAMENTO CILINDRICO PERCIÒ uso direttamente la matrice polare
    // definisco una lista di punti utili per l'analisi, cioè con valore non nullo e sotto la bright band o sopra (NB. per l'analisi considero i punti usati per la ZLR che si suppongono non affetti da clutter e beam blocking<50%. questo ha tutta una serie di implicazioni..... tra cui la proiezione, il fatto che nel confronto entrano dei 'buchi'...cioè il confronto è fatto su una matrice pseudo-orizzontale bucata e quote variabili tuttavia, considerato che i gradienti verticali in caso convettivo e fuori dalla bright band non sono altissimi spero funzioni ( andro' a verificare che l'ipotesi sia soddisfatta)
    // per trovare il background uso uno pseudo quadrato anzichè cerchio, mantenendo come semi-lato il raggio di steiner (11KM); devo perciò  definire una semi-ampiezza delle finestre in azimut e range corrispondente al raggio di 11 km
    // quella in azimut  dipende dal range

{
    using namespace steiner;

    if (lista_bckg.size() < 2)
        return;

    // Per il calcolo della finestra range su cui calcolare il background
    // divido il raggio di Steiner (11km) per la dimensione della cella.
    // Definisco ampiezza semi-finestra range corrispondente al raggio di
    // steiner (11km), in unità matrice polare.
    unsigned delta_nr = round(STEINER_RADIUS * 1000. / size_cell);
    LOG_DEBUG("delta_n range per analisi Steiner = %u  --- dimensione lista_bckg %d", delta_nr,lista_bckg.size());

    for (vector<Point>::iterator i = lista_bckg.begin(); i != lista_bckg.end(); ++i)
    {
        // estremi della finestra in range
        int kmin = i->range - delta_nr;
        unsigned kmax = min(i->range + delta_nr, max_bin);

        if (kmin>0)
        {
            //definisco ampiezza semi finestra nazimut  corrispondente al raggio di steiner (11km)  (11/distanzacentrocella)(ampiezzaangoloscansione)
            unsigned delta_naz=ceil(STEINER_RADIUS/((i->range * size_cell/1000. + size_cell/2000.)/(AMPLITUDE*DTOR)));
            if (delta_naz > NUM_AZ_X_PPI / 2)
                delta_naz = NUM_AZ_X_PPI / 2;

            int jmin = i->azimut - delta_naz;
            int jmax = i->azimut + delta_naz;

            for (int j = jmin; j < jmax; ++j)
                for (unsigned k = kmin; k < kmax; ++k)
                    i->add_sample(elev_fin.db_at_elev_preci((j + NUM_AZ_X_PPI) % NUM_AZ_X_PPI, k));
        } else {
            // FIXME: questo fa mezzo scan tra 0 e kmax, e mezzo scan tra 0 e
            // -kmin. Sempre gli stessi mezzi scan a prescindere dalla
            // posizione di i. Ha senso?
            for (unsigned j=0   ; j<NUM_AZ_X_PPI/2  ; j++)
                for (unsigned k=0  ; k<kmax   ; k++)
                    i->add_sample(elev_fin.db_at_elev_preci(j, k));

            for (unsigned j= NUM_AZ_X_PPI/2  ; j<NUM_AZ_X_PPI  ; j++)
                for (int k=0  ; k<-kmin   ; k++)
                    i->add_sample(elev_fin.db_at_elev_preci(j, k));
        }
        i->finalize();
    }
}

void CalcoloSteiner::ingrasso_nuclei(float cr,int ja,int kr)
{
    int dr=(int)(cr*1000./size_cell);//definisco ampiezza semi-finestra range corrispondente al raggio di steiner (11km), unità matrice polare

    int kmin=kr-dr;
    unsigned kmax=kr+dr;

    int daz=ceil(cr/((kr*size_cell/1000.+size_cell/2000.)/(AMPLITUDE*DTOR)));
    int jmin=ja-daz;
    unsigned jmax=ja+daz;

    LOG_DEBUG("dr cr kmin kmax  %d %f %d %d %d %d", dr,cr, kmin,kmax,jmin,jmax);

    if (kmin>0) {
        if (kmax > max_bin) kmax = max_bin;

        if (jmin<0) {
            jmin=NUM_AZ_X_PPI-jmin%NUM_AZ_X_PPI;
            for (unsigned j=jmin; j< NUM_AZ_X_PPI ; j++) {
                for (unsigned k=kmin ; k<kmax  ; k++) {
                    conv_STEINER(j, k)=CONV_VAL;
                }
            }
            LOG_DEBUG("jmin %d", jmin);
            jmin=0;

        }

        if (jmax>=NUM_AZ_X_PPI) {
            jmax=jmax%NUM_AZ_X_PPI;
            for (unsigned j=0; j<jmax ; j++) {
                for (unsigned k=kmin; k<kmax  ; k++) {
                    conv_STEINER(j, k)=CONV_VAL;
                }
            }
            LOG_DEBUG("jmax %d", jmax);
            jmax=NUM_AZ_X_PPI;
        }
        for (unsigned j=jmin; j<jmax ; j++) {
            for (unsigned k=kmin; k<kmax  ; k++) {
                conv_STEINER(j%NUM_AZ_X_PPI, k)=CONV_VAL;
            }
        }
    }
    else
    {
        for (unsigned j=0   ; j<NUM_AZ_X_PPI/2  ; j++)
            for (unsigned k=0  ; k<kmax   ; k++){
                conv_STEINER(j, k)=CONV_VAL;
            }
        for (unsigned j= NUM_AZ_X_PPI/2  ; j<NUM_AZ_X_PPI  ; j++)
            for (unsigned k=0  ; k < (unsigned)-kmin   ; k++){
                conv_STEINER(j, k)=CONV_VAL;
            }
    }
    LOG_DEBUG ("Finita ingrasso nuclei");
    return;
}

void CalcoloSteiner::classifico_STEINER()
{
    using namespace steiner;

    for (vector<Point>::const_iterator i = lista_bckg.begin(); i != lista_bckg.end(); ++i)
    {
        int j = i->azimut;
        int k = i->range;
        if (j < 0 || k < 0) continue;

        double db = elev_fin.db_at_elev_preci(j, k);
        // calcolo diff col background
        float diff_bckgr = db - i->bckgr;
        // test su differenza con bckground , se soddisfatto e simultaneamente il VIZ non ha dato class convettiva (?)
        if ((db > 40.) ||
            (i->bckgr < 0 && diff_bckgr > 10) ||
            (i->bckgr < 42.43 && i->bckgr > 0 && diff_bckgr > 10. - i->bckgr * i->bckgr / 180.) ||
            (i->bckgr > 42.43 && diff_bckgr > 0))
        {
            // assegno il punto nucleo di Steiner
            conv_STEINER(j, k) = CONV_VAL;

            // ingrasso il nucleo
            float cr = i->convective_radius;
            LOG_DEBUG(" %f cr", cr);
            ingrasso_nuclei(cr, j, k);
        }
    }
    return;
}

}
}
