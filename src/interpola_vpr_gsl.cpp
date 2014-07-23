#include "interpola_vpr.h"
#include "cum_bac.h"
#include "logging.h"
#include <vpr_par.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

// Funzioni di supporto a interpola_vpr
namespace {

struct data {
  const unsigned n;
  double *t;
  double *y;
  double *sigma;
  data(unsigned n)
    : n(n), t(new double[n]), y(new double[n]), sigma(new double[n])
  {
  }
  ~data()
  {
    delete[] t;
    delete[] y;
    delete[] sigma;
  }
  void print()
  {
      fprintf(stderr, "i    t    y    sigma\n");
      for (unsigned i = 0; i < n; ++i)
          fprintf(stderr, "%2d %.2f %.2f %.2f\n", i, t[i], y[i], sigma[i]);
  }
};

int
expb_f (const gsl_vector * x, void *data, 
        gsl_vector * f)
{
  const struct data& args = *(struct data*)data;
  double B = gsl_vector_get (x, 0);
  double E = gsl_vector_get (x, 1);
  double G = gsl_vector_get (x, 2);
  double C = gsl_vector_get (x, 3);
  double F = gsl_vector_get (x, 4);
  for (unsigned i = 0; i < args.n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      //double t = i;
      double Yi = B * exp (- (args.t[i]-E)/G*(args.t[i]-E)/G ) + C+F*args.t[i];
      gsl_vector_set (f, i, (Yi - args.y[i])/args.sigma[i]);
    }
  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *data, 
         gsl_matrix * J)
{
  const struct data& args = *(struct data*)data;
  double B = gsl_vector_get (x, 0);
  double E = gsl_vector_get (x, 1);
  double G=  gsl_vector_get (x, 2);
  for (unsigned i = 0; i < args.n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      // double t = i;
      //    double e = exp(-lambda * t);
      double arg =((args.t[i]-E)/G);
      double ex=exp(-arg*arg);
      double fac=B*ex*2.0*arg;

      gsl_matrix_set (J, i, 0, ex); 
      gsl_matrix_set (J, i, 1, fac/G);
      gsl_matrix_set (J, i, 2, fac*arg/G);
      gsl_matrix_set (J, i, 3, 1.); 
      gsl_matrix_set (J, i, 4, args.t[i]); 

    }
  return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *data,
          gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, data, f);
  expb_df (x, data, J);

  return GSL_SUCCESS;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  fprintf(stderr, "iter: %3zu x = % 15.8f % 15.8f % 15.8f  % 15.8f  % 15.8f  "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->x, 4), 
          gsl_blas_dnrm2 (s->f));
}

int testfit(double a[])
{
    if (a[0]<0. || a[0] >15.) return 1;
    if (a[1] >10.) return 1;
    if (a[2]<0.2 || a[2] > 0.6 ) return 1; //da analisi set dati
    if (a[3]<0. ) return 1;
    if (a[4]>0 ) return 1;
    return 0;
}

float lineargauss(double xint , double a[])
{
 return a[0] * exp (- (xint-a[1])/a[2]*(xint-a[1])/a[2] ) + a[3]+a[4]*xint;

}

}


namespace elaboradar {

int InterpolaVPR_GSL::interpola_VPR(const float* vpr, int hvprmax, int livmin)
{
    LOG_CATEGORY("radar.vpr");
    static const unsigned N = 10;
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    unsigned int i;
    const size_t n = N;
    const size_t p = 5;
    char file_vprint[512];
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double a[5];
    struct data d(N);
    gsl_multifit_function_fdf f;
    double x_init[5] = { 4, 0.2, 3. , 1.4, -0.4 };
    gsl_vector_view x = gsl_vector_view_array (x_init, p);

    //////////////////////////////////////////////////////////////////////////////
    int ier_int=0;
    double xint,yint;
    /* punti interessanti per inizializzare parametri*/
    int  in1=(int)((hvprmax-TCK_VPR/2)/TCK_VPR); //indice del massimo
    int  in2=(int)((hvprmax+HALF_BB)/TCK_VPR); //indice del massimo + 500 m
    int  in3=in2+1;
    int  in4=in2+5; //indice del massimo + 1000 m
    if (in4 > NMAXLAYER-1) {
        ier_int=1;
        return ier_int;
    }

    B=vpr[in1]-vpr[in2];
    E=hvprmax/1000.;
    G=0.25;
    C=vpr[in2-1];
    F=vpr[in4]<vpr[in3]?(vpr[in4]-vpr[in3])/((in4-in3)*TCK_VPR/1000.):0.;
    // fprintf(stderr, "const unsigned NMAXLAYER=%d;\n", NMAXLAYER);
    // fprintf(stderr, "float vpr[] = {");
    // for (unsigned i = 0; i < NMAXLAYER; ++i)
    //     fprintf(stderr, "%s%f", i==0?"":",", (double)vpr[i]);
    // fprintf(stderr, "};\n");

    x_init[0]= a[0]=B;
    x_init[1]= a[1]=E;
    x_init[2]= a[2]=G;
    x_init[3]= a[3]=C;
    x_init[4]= a[4]=F;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    /* This is the data to be fitted */

    for (i = 0; i < n; i++)
    {
        d.t[i]= ((hvprmax-1000.)>livmin)? (i*TCK_VPR+(hvprmax-800)-TCK_VPR)/1000. : (livmin+i*TCK_VPR)/1000.;
        d.y[i]= ((hvprmax-1000.)>livmin)? vpr[i+(int)(((hvprmax-800)-TCK_VPR)/TCK_VPR)] : vpr[i+(int)(livmin/TCK_VPR)];
        d.sigma[i] = 0.5;
    };

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    //print_state (0, s);
    bool found = false;
    for (unsigned iter = 0; !found && iter < 500; ++iter)
    {
        //fprintf(stderr, "Iter %d\n", iter);
        //d.print();
        int status = gsl_multifit_fdfsolver_iterate (s);
        if (status != 0)
        {
            LOG_ERROR("gsl_multifit_fdfsolver_iterate: %s", gsl_strerror(status));
            return 1;
        }

        //print_state (iter, s);

        status = gsl_multifit_test_delta (s->dx, s->x,
                1e-4, 1e-4);
        switch (status)
        {
            case GSL_SUCCESS: found = true; break;
            case GSL_CONTINUE: break;
            default:
                LOG_ERROR("gsl_multifit_test_delta: %s", gsl_strerror(status));
                return 1;
        }
    }

    gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    { 
        double chi = gsl_blas_dnrm2(s->f);
        double dof = n - p;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

        // printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

        // printf ("B      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        // printf ("E = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
        // printf ("G     = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
        // printf ("C = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
        // printf ("F     = %.5f +/- %.5f\n", FIT(4), c*ERR(4));
    }

    B = a[0] = FIT(0);
    E = a[1] = FIT(1);
    G = a[2] = FIT(2);
    C = a[3] = FIT(3);
    F = a[4] = FIT(4);

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);

    /////////////////////////////////////////////////////////

    if (testfit(a) == 1)
        return 1;

    for (i=1; i<=N; i++)
    {
        xint=(i*TCK_VPR-TCK_VPR/2)/1000.;
        yint= lineargauss(xint, a);
        vpr_int[i-1] = yint;
    }

    return 0;
}

}
