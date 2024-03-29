#include "dbz.h"
#include "volume.h"

//#define  aMP 316.
//#define  bMP 1.5
#define  aMP_conv 500.0 
#define  bMP_conv 1.5
#define  aMP_strat 250.
#define  bMP_strat 1.5
#define  aMP_SNOW 400.0
#define  bMP_SNOW 2.0
#define  aMP_class 200.
#define  bMP_class 1.6

#define THRES_ATT 0 /* minimo valore di Z in dBZ per calcolare att rate */

namespace radarelab {
namespace algo {

using namespace std;

DBZ::DBZ(const Volume<double>& volume)
{
    // --- ricavo il mese x definizione first_level e  aMP bMP ---------
    time_t Time = volume.load_info->acq_date;
    struct tm* tempo = gmtime(&Time);
    init(tempo->tm_mon+1, volume[0].cell_size);
}

DBZ::DBZ(int month, double base_cell_size)
{
    init(month, base_cell_size);
}

void DBZ::init(int month, double base_cell_size)
{
    this->base_cell_size = base_cell_size;

    if (month > 4 && month < 10)
    {
        aMP = aMP_conv;
        bMP = bMP_conv;
    }
    else
    {
        aMP = aMP_strat;
        bMP = bMP_strat;
    }
}

double DBZ::attenuation(unsigned char DBZbyte, double  PIA)  /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/
{
    double Zhh,att_rate,R;/* PIA diventa att_tot devo decidere infatti se PIA sarà 3d percio' temp. uso  nomi diversi*/
    double att_tot;

    //---ricevo in ingresso il dato e l'attenuazione fino  quel punto
    //---la formula recita che l'attenuazione è pari una funzione di Z reale (quindi corretta dell'attenuazione precedente). ovviamente devo avere un segnale per correggere.
    //--------- CALCOL
    att_tot=PIA;
    Zhh= (BYTEtoZ(DBZbyte));
    if (10*log10(Zhh) > THRES_ATT )
    {
        Zhh=pow(10., (log10(Zhh)+ 0.1*att_tot));
        R=pow((Zhh/aMP),(1.0/bMP));
        att_rate=0.0018*pow(R,1.05);
        // TODO: to compute scan by scan?
        att_tot=att_tot+2.*att_rate*0.001 * base_cell_size;
        if (att_tot>BYTEtoDB(254)) att_tot=BYTEtoDB(254);
    }
    return att_tot;
}
double DBZ::attenuation(double DBZvalue, double  PIA)  /* Doviak,Zrnic,1984 for rain as reported in cost 717 final document*/
{
    double Zhh,att_rate,R;/* PIA diventa att_tot devo decidere infatti se PIA sarà 3d percio' temp. uso  nomi diversi*/
    double att_tot;

    //---ricevo in ingresso il dato e l'attenuazione fino  quel punto
    //---la formula recita che l'attenuazione è pari una funzione di Z reale (quindi corretta dell'attenuazione precedente). ovviamente devo avere un segnale per correggere.
    //--------- CALCOL
    att_tot=PIA;
    Zhh=(DBZtoZ(DBZvalue));
    if (DBZvalue > THRES_ATT )
    {
        //Zhh=pow(10., (log10(Zhh)+ 0.1*att_tot));
        Zhh=DBZtoZ(DBZvalue+ att_tot);
        R=pow((Zhh/aMP),(1.0/bMP));
        att_rate=0.0018*pow(R,1.05);
        // TODO: to compute scan by scan?
        att_tot=att_tot+2.*att_rate*0.001 * base_cell_size;
        if (att_tot>BYTEtoDB(254)) att_tot=BYTEtoDB(254);
    }
    return att_tot;
}

double DBZ::RtoDBZ(double rain) const
{
    return RtoDBZ(rain, aMP, bMP);
}

double DBZ::DBZtoR(double dbz) const
{
    return DBZtoR(dbz, aMP, bMP);
}

double DBZ::DBZ_snow(double dbz) const
{
    return RtoDBZ(DBZtoR(dbz, aMP_SNOW, bMP_SNOW), aMP_class, bMP_class);
}

double DBZ::DBZ_conv(double dbz) const
{
    return RtoDBZ(DBZtoR(dbz, aMP_conv, bMP_conv), aMP_class, bMP_class);
}

double DBZ::RtoDBZ_class(double R) const
{
    return RtoDBZ(R, aMP_class, bMP_class);
}

double DBZ::DBZ_to_mp_func(double sample) const
{
    return DBZtoR(sample, aMP, bMP);
}

double DBZ::BYTEtoZ(unsigned char byte)
{
    const double gain = 80. / 255.;
    const double offset = -20.;
    static bool precomputed = false;
    static double Z[256];

    if (!precomputed)
    {
        for (unsigned i=0; i < 256; ++i)
            Z[i] = pow(10., (i * gain + offset) * 0.1);
        precomputed = true;
    }

    return Z[byte];
}

}
}
