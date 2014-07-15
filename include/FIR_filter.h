#include <fftw3.h>
#include <cmath>
#include <iostream>

using namespace std;

class FIR_filter
{
public:
	double *in,*re_in;//,*re_in2;
	fftw_complex *out;
	fftw_plan plan_forw,plan_back;//plan_back2;
	unsigned n;

	FIR_filter(unsigned N,double* ret)
	{
		n=N;	
		in = (double*) fftw_malloc(sizeof(double)*N);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(1+floor(0.5*N)));
		re_in = ret; // punto re_in (che è il risultato) su un puntatore passato al costruttore che deve avere n double allocati

		plan_forw = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE); // sign is not need for real trasforms
		plan_back = fftw_plan_dft_c2r_1d(N, out, re_in, FFTW_MEASURE);
	}

	void feed(double* input)
	{
		copy(input,input+n,in);	// la copia è necessaria, in deve essere costante per plan_forw ed fftw_MEASURE ne distrugge il contenuto
					// magari esiste un metodo per evitare la copia, ma non lo conosco
					// si suppone che input abbia n double. Un numero diverso di n richiede l'istanza di un nuovo oggetto FIR
					// per poter avere i dovuti plan.
	}

	~FIR_filter()
	{
		fftw_destroy_plan(plan_forw);
		fftw_destroy_plan(plan_back);

		fftw_free(in);
		fftw_free(out);
	}	
		
	void perform()
	{
		fftw_execute(plan_forw);
		frequency_filter();
		frequency_filter();
		fftw_execute(plan_back);		
	}

	void dump()
	{
		cout<<endl<<endl;
		for(unsigned i=0;i<20;i++)
		{
			cout<<fixed<<in[i]<<"\t"<<out[i][0]<<"\t"<<out[i][1]<<"\t\t"<<re_in[i]/n<<endl;
		}
		cout<<endl<<endl;
	}
private:
	void frequency_filter()
	{
		double hn[20];
		hn[0]=hn[20]=1.625807356e-2;
		hn[1]=hn[19]=2.230852545e-2;
		hn[2]=hn[18]=2.896372364e-2;
		hn[3]=hn[17]=3.595993808e-2;
		hn[4]=hn[16]=4.298744446e-2;
		hn[5]=hn[15]=4.971005447e-2;
		hn[6]=hn[14]=5.578764970e-2;
		hn[7]=hn[13]=6.089991897e-2;
		hn[8]=hn[12]=6.476934523e-2;
		hn[9]=hn[11]=6.718151185e-2;
		hn[10]=6.8001e-2;
		fftw_complex* Hz=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(1+floor(0.5*n)));
		double omega;
		double real,imag;
		for(unsigned i=0;i<1+floor(0.5*n);i++)
		{
			Hz[i][0]=0;
			Hz[i][1]=0;
			omega=i*M_PI/(1+floor(0.5*n));
			for(unsigned j=0;j<20;j++)
			{
				Hz[i][0]+=hn[j]*cos(-omega*j);
				Hz[i][1]+=hn[j]*sin(-omega*j);
			}
			real=out[i][0]*Hz[i][0]-out[i][1]*Hz[i][1];
			imag=out[i][0]*Hz[i][1]+out[i][1]*Hz[i][0];
			out[i][0]=real;
			out[i][1]=imag;
		}
	}
};
