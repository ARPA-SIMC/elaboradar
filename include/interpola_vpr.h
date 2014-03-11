#include "vpr_par.h"

namespace cumbac {

struct CalcoloVPR;

struct InterpolaVPR
{
    // Output parameters
    double B, E, G, C, F;
    double chisqfin;
    double rmsefin;
    double vpr_int[NMAXLAYER];

    InterpolaVPR();
    ~InterpolaVPR();

    /**
     *  @brief   funzione che esegue interpolazione del profilo
     *  @details interpola profilo usando una funzione gaussiana+lineare  y= B*exp(-((x-E)/G)^2)+C+Fx 
     *  @param[out] a[] vettore dei parametri della funzione
     *  @param[out] ma  dimensione vettore parametri
     *  @return ier_int codice di uscita (0=ok 1=fallito)
     */
    virtual int interpola_VPR(const float* vpr, int hvprmax, int livmin) = 0;
};

struct InterpolaVPR_NR : public InterpolaVPR
{
    virtual int interpola_VPR(const float* vpr, int hvprmax, int livmin);
};

struct InterpolaVPR_GSL : public InterpolaVPR
{
    virtual int interpola_VPR(const float* vpr, int hvprmax, int livmin);
};

}
