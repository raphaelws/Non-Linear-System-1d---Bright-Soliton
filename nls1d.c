#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define PI 3.141592653589793
#define x_min -64
#define x_max 64
#define dt 0.01*I
#define dt2 0.01
#define N 512 //power of 2, good for FFT
#define x0 0
#define time 100

//run with gcc nls1d.c -o nls1d -lfftw3 -lm

//if imaginary time dt = Idt, if real time dt = dt

double complex Gaussian(double);
void normalize();
void initialize();
void write_file(FILE *);
void evolve();
double V(double);

double complex V_exp[N];
double complex T_exp[N];
double complex psi[N];   
double complex psip[N]; 
double p0;
double x_step;
fftw_plan fft, ifft;    // fftw plan

int main()
{
	x_step = (double)(x_max-x_min)/(double)N;
	FILE *f;
	int i;
	double T, pl, pr,k=0;
	
	f = fopen("groundstate.txt", "w+");
	
	fft = fftw_plan_dft_1d(N, psi, psip, FFTW_FORWARD,  FFTW_ESTIMATE);
	ifft = fftw_plan_dft_1d(N, psip, psip, FFTW_BACKWARD,  FFTW_ESTIMATE);
	initialize();
	while(k<=time){
		k=k+dt2;
		evolve();
		printf("Steps take = %lf", k);
		printf("\n");
	}
		
	
	fftw_destroy_plan(fft);
	fftw_destroy_plan(ifft);
	
	write_file(f);
	fclose(f);
	return 0;
}

double V(double x)
{
	
		return pow(x,2);
}

void initialize()
{
	int i;
	double x, p;
	for(i=0; i<N; i++){
		x = x_min + x_step*(i+1);
		//V_exp[i]=cexp(-I*V(x)*dt);
		psi[i] = Gaussian(x);
		V_exp[i]=cexp(I*psi[i]*conj(psi[i]*dt/2.0));
		
	}

	for(i=0; i<N; i++){
		// Considering Sampling Theorem && Cyclic displacement
		if(i<N/2) 
			p = i*2.0*PI/x_step/N;
		else      
			p = (i-N)*2.0*PI/x_step/N;
		T_exp[i] = cexp(-I*p*p*dt/2.0);
	}
	normalize();

}

void write_file(FILE *f)
{
	int i;
	double wr,wi,x;
	for(i=0; i<N; i++){
	    	x = x_min + x_step*(i+1);
		wr = creal(psi[i]);
		wi = cimag(psi[i]);
		//writes position x and psi square
		fprintf(f,"%lf " ,x);
		fprintf(f," ");
		fprintf(f,"%lf " ,wr*wr+wi*wi);
		fprintf(f, "\n");
	}
}

double complex Gaussian(double x)
{
	return exp(-(x*x));
}

void normalize()
// This function nomalizes the wave packet numerically
{
	double ww[N];
	double prod=0;
	double A;
	int i;
	for(i=0; i<N; i++){
		ww[i] = psi[i]*conj(psi[i]);
	}
	//trapezoidal integral
	for(i=0; i<N-1; i++){
		prod += (ww[i]+ww[i+1])*x_step/2.0;
	}
	A=sqrt(prod);
	for(i=0; i<N; i++){
		psi[i] = psi[i]/A;
	}
}

void evolve()

{
	//envolve half nabla 2
	int i;
	for(i=0; i<N; i++){
		psi[i] *= V_exp[i];
	}
	//fft
	fftw_execute(fft);
	//envolve the knetic therm
	for(i=0; i<N; i++){
		psip[i] *= T_exp[i];
	}
	//inverse fft
	fftw_execute(ifft);
	
	//envolve half nabla 2
	for(i=0; i<N; i++){
		psi[i] *= V_exp[i];
	}
	normalize();
}





