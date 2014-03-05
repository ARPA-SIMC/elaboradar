namespace cumbac {

struct CalcoloVPR;

struct InterpolaVPR
{
    double B, E, G, C, F;
    double chisqfin;
    double rmsefin;

    InterpolaVPR();
    ~InterpolaVPR();

    /**
     *  @brief   funzione che esegue interpolazione del profilo
     *  @details interpola profilo usando una funzione gaussiana+lineare  y= B*exp(-((x-E)/G)^2)+C+Fx 
     *  @param[out] a[] vettore dei parametri della funzione
     *  @param[out] ma  dimensione vettore parametri
     *  @return ier_int codice di uscita (0=ok 1=fallito)
     */
    virtual int interpola_VPR(const CalcoloVPR& cv) = 0;
};

struct InterpolaVPR_NR : public InterpolaVPR
{
    virtual int interpola_VPR(const CalcoloVPR& cv);
};

struct InterpolaVPR_GSL : public InterpolaVPR
{
    virtual int interpola_VPR(const CalcoloVPR& cv);
};

}
