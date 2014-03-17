#include "interpola_vpr.h"
#include "cum_bac.h"
#include "logging.h"

#ifdef __cplusplus
extern "C" {
#endif
// nr
#include <nrutil.h>
#include <nr.h>
#ifdef __cplusplus
}
#endif

# include <vpr_par.h>

namespace {

/**
 *  restituisce il valore della combinazione tra gaussiana e lineare
 *  @brief   funzione che  calcola derivate della gaussiana + lineare rispetto ai parametri e valore (*y) in un punto x 
 *  @details    y(x,a) is the sum of a gaussian and a linear function with amplitude B=a[1], center E=a[2] and width G=a[3] and a linear function with coefficients C(shift)=a[4] and F(slope)=a[5] 
 *  a is the parameters vector and dyda is the vector of derivatives respect to the different parameters
 *  @param[in] x vettore delle quote
 *  @param[in] a vettore dei parametri
 *  @param[in] y vettore dei valori della funzione 
 *  @param[out]  dyda derivate
 *  @return 0 codice di uscita 0
 */
void lineargauss(float x, float a[], float *y, float dyda[], int na)
{
  float fac, ex, arg;

  *y=0.0;
  arg=(x-a[2])/a[3];
  ex=exp(-arg*arg);
  fac=a[1]*ex*2.0*arg;
  *y+=a[1]*ex+a[4]+a[5]*x;
  dyda[1]=ex;
  dyda[2]=fac/a[3];
  dyda[3]=fac*arg/a[3];
  dyda[4]=1.;
  dyda[5]=x;
}

/**
 *  testa i parametri del fit in modo che abbiano significato fisico 
 *
 *  @brief   funzione che testa il fit dell'interpolazione del profilo
 *  @details verifica che i parametri del fit del profilo abbiano senso
 *  @param[in] a[] vettore dei parametri della funzione
 *  @param[in] chisq  chiquare
 *  @return codice di uscita 0
 *
 */
int testfit(float a[], float chisq, float chisqin)
{
    if (a[1]<0. || a[1] >15.) return 1;
    if (a[2] >10.) return 1;
    if (a[3]<0.2 || a[3] > 0.6 ) return 1; //da analisi set dati
    if (a[4]<0. ) return 1;
    if (a[5]>0 ) return 1;
    if (chisq>chisqin ) return 1;
    return 0;
}

}

namespace cumbac {

/*
   comstart interpola_VPR
   idx interpola il profilo verticale tramite una funzione lingauss
   interpola il profilo verticale tramite una funzione gaussiana + lineare del tipo

   y= B*exp(-((x-E)/G)^2)+C+Fx

   usa la funzione mrqmin delle numerical recipes in C: tutti i vettori passati a mrqmin devono essere allocati e deallcocati usando le funzioni di NR (vector, matrix, free_vector, free_matrix.. etc) che definiscono vettori con indice a partire da 1.
   NB gli ndata dati considerati partono da 1000 m sotto il massimo (in caso il massimo sia più basso di 1000 m partono da 0 m)
   A ogni iterazione si esegue un test sui parametri. Se ritorna 1 si torna ai valori dell'iterazione precedente.
   A fine interpolazione si verifica che il chisquare non superi una soglia prefissata, in tal caso ritorna 1 e interpol. fallisce.

   INIZIALIZZAZIONE PARAMETRI:
   a[1]=B=vpr(liv del massimo)-vpr(liv. del massimo+500m);
   a[2]=E= quota liv. massimo vpr (in KM);
   a[3]=G=semiampiezza BB (quota liv .massimo- quota massimo decremento vpr nei 600 m sopra il massimo ) in KM;
   a[4]=C=vpr(liv massimo + 700 m);
   a[5]=F=coeff. angolare segmento con estremi nel vpr ai livelli max+900m e max+1700m, se negativo =0.;


   float a[ma], int ma: vettore parametri e n0 parametri
   float *x, *y:  quote in KM e valori del vpr usati per l'interpolazione
   float *sig,alamda : vettore dev st. e variabile che decrementa al convergere delle iterazioni
   float *dyda: vettore derivate rispetto ai parametri
   float B,E,C,G,F:  parametri da ottimizzare, first guess
   float chisq; scarto quadratico
   int i,in1,in2,in3,in4,*ia,ifit,ii,ndati_ok,k;
   int ndata=15;  numero di dati considerati
   float **covar,**alpha; matrice ovarianze, matrice alpha

   comend
*/
int InterpolaVPR_NR::interpola_VPR(const float* vpr, int hvprmax, int livmin)
{
    LOG_CATEGORY("radar.vpr");
    static const unsigned npar=5;
    float *x, *y,*sig,alamda,y1=0,*dyda,xint,qdist,*abest;
    float chisq=100.;
    float chisqold=0.0;
    float chisqin=0.0;
    int i,in1,in2,in3,in4,*ia,ifit,ii,ndati_nok,k,ier_int;
    //int ma=5;
    int ndata=10;
    float **covar;
    float **alpha;
    FILE *file;

    float *a=vector(1,npar);
    for (i=1;i<=npar;i++){
        a[i]=NODATAVPR;
    }

    LOG_INFO("sono in interpola_vpr");
    ier_int=0;

    in1=(hvprmax-TCK_VPR/2)/TCK_VPR; //indice del massimo
    in2=(hvprmax+HALF_BB)/TCK_VPR; //indice del massimo + 500 m
    in3=in2+1;
    in4=in2+5; //indice del massimo + 1000 m
    LOG_INFO("in1 in2 %i %i %f %f",in1,in2,vpr[in1],vpr[in2]);

    if (in4 > NMAXLAYER-1) {
        ier_int=1;
        return ier_int;
    }

    /* inizializzazione vettore parametri */
    abest=vector(1,npar);
    x=vector(1,ndata);
    y=vector(1,ndata);
    sig=vector(1,ndata);
    ia=ivector(1,npar);
    covar=matrix(1,npar,1,npar);
    alpha=matrix(1,npar,1,npar);

    for (k=in1+2; k<=in3; k++)
    {
        ier_int=0;

        dyda=vector(1,npar);
        a[1]=B=vpr[in1]-vpr[in2];
        a[2]=E=hvprmax/1000.;
          a[3]=G=(k-in1-0.5)*TCK_VPR/1000.;
          //  a[3]= G=0.25;
        a[4]=C=vpr[in2];
        a[5]=F=vpr[in4]<vpr[in3]?(vpr[in4]-vpr[in3])/((in4-in3)*TCK_VPR/1000.):0.;
        //fprintf(stderr, "k:%d, a1:%f a2:%f a3:%f a4:%f a5:%f\n", k, a[1], a[2], a[3], a[4], a[5]);

        alamda=-0.01;

        for (i=1;i<=npar;i++) ia[i]=1;
        qdist=0;
        ii=1;
        ndati_nok=0;

        for (i=1; i<=ndata; i++)
        {
            sig[ii]=0.5;
            x[ii]= ((hvprmax-1000.)>livmin)? (i*TCK_VPR+(hvprmax-800)-TCK_VPR)/1000. : (livmin+(i-1)*TCK_VPR)/1000.;
            y[ii]= ((hvprmax-1000.)>livmin)? vpr[i+((hvprmax-800)-TCK_VPR)/TCK_VPR] : vpr[i-1+livmin/TCK_VPR];
            // x[ii]= ((hvprmax-800.)>livmin)? (i*TCK_VPR+(hvprmax-600)-TCK_VPR)/1000. : (livmin+(i-1)*TCK_VPR)/1000.;
            //y[ii]= ((hvprmax-800.)>livmin)? vpr[i+((hvprmax-600)-TCK_VPR)/TCK_VPR] : vpr[i-1+livmin/TCK_VPR];
            lineargauss(x[ii], a, &y1, dyda, ndata);
            qdist=(y1-y[ii])*(y1-y[ii]);
            //fprintf(stderr, "i:%d, ii:%d, xii:%f, yii:%f, y1:%f, qdist:%f, chisqin:%f\n", i, ii, x[ii], y[ii], y1, qdist, chisqin);
            if (sqrt(qdist) < DIST_MAX)
            {
                ii+=1;
                chisqin=qdist+chisqin;
            }
            else
                ndati_nok=ndati_nok+1;

            if    ( ndati_nok > 2  )
            {
                LOG_WARN("  first guess troppo lontano dai punti , interpolazione fallisce");
                ier_int=1;
                break;
            }
        }

        if (!ier_int)
        {
            LOG_INFO("\n alamda %f   chisqin % f a[1]  % f a[2] % f a[3]  % f a[4]  % f a[5]  % f", alamda,chisqin,a[1], a[2], a[3],a[4],a[5]);
            //fprintf(stderr, "i    t    y    sigma\n");
            //for (unsigned i = 1; i <= ndata; ++i)
                //fprintf(stderr, "%2d %.2f %.2f %.2f\n", i, x[i], y[i], sig[i]);
            ifit=0;
            while (fabs(chisq-chisqold) > DCHISQTHR && ifit < 1)
            {
                chisqold=chisq;
                B=a[1];
                E=a[2];
                G=a[3];
                C=a[4];
                F=a[5];
                mrqmin(x, y, sig, ndata, a, ia, npar, covar, alpha, &chisq, &lineargauss, &alamda);
                LOG_INFO("alamda %f   chisq % f a[1]  % f a[2] % f a[3]  % f a[4]  % f a[5]  % f", alamda,chisq,a[1], a[2], a[3],a[4],a[5]);
                ifit=testfit(a,chisq,chisqin);  /*test sul risultato del fit */
                if (ifit)
                {/*test sul risultato del fit */
                    a[1]=B; /*test sul risultato del fit */
                    a[2]=E; /*test sul risultato del fit */
                    a[3]=G; /*test sul risultato del fit */
                    a[4]=C; /*test sul risultato del fit */
                    a[5]=F; /*test sul risultato del fit */
                    chisq=chisqin;
                }
            }


            if (chisq < chisqfin)
            {
                chisqfin=chisq;
                for (i=1;i<=npar;i++) abest[i]=a[i];
            }
        }
    }

    for (i=1; i<=ndata-ndati_nok; i++)
    {
        lineargauss(x[i], abest, &y1, dyda, ndata);
        rmsefin=rmsefin+  (y[i]-y1)*(y[i]-y1) ;
    }
    rmsefin=sqrt(rmsefin/(float)((ndata-ndati_nok)*(ndata-ndati_nok)));
    LOG_INFO("RMSEFIN %f", rmsefin );


    if (chisqfin>CHISQ_MAX)
    {
        ier_int=1;
    }
    else {
        // Calcola il profilo interpolato
        for (i=1;i<=npar;i++) a[i]=abest[i];
        for (i=1; i<=NMAXLAYER; i++)
        {
            xint=(i*TCK_VPR-TCK_VPR/2)/1000.;
            lineargauss(xint, a, &y1, dyda, ndata);
            vpr_int[i-1] = y1;
        }
    }
    B=a[1];
    E=a[2];
    G=a[3];
    C=a[4];
    F=a[5];
    free_vector(dyda,1,npar);
    free_vector(abest,1,npar);
    free_ivector(ia,1,npar);
    free_vector(x,1,ndata);
    free_vector(y,1,ndata);
    free_vector(sig,1,ndata);
    free_matrix(alpha,1,ndata,1,ndata);
    free_matrix(covar,1,ndata,1,ndata);
    free_vector(a,1,npar);

    return ier_int;
}

}
