#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define PI 3.141592653589793
#define x_min -40
#define x_max 40
#define dt 0.01*I
#define dt2 0.01
#define N 1024
#define x0 0
#define time 100
//time 0.17

//run with gcc nls1d.c -o nls1d -lfftw3 -lm

//if imaginary time dt = Idt, if real time dt = dt


//i*dpsi/dt = (-d^2/dx^2 + U(x) + g|psi|^2)psi
double complex psi0(double);
void normalize();
void attvexp();
void initialize();
void write_file(FILE *);
void evolve();
double V(double);

double complex V_exp[N];
double complex T_exp[N];
double complex psi[N];   
double complex psip[N]; 
double x_step;
fftw_plan fft, ifft;    // fftw plan

int main()
{
	x_step = (double)(x_max-x_min)/(double)N;
	FILE *f;
	int i;
	double k=0;
	
	f = fopen("groundstate.txt", "w+");;
	
	fft = fftw_plan_dft_1d(N, psi, psip, FFTW_FORWARD,  FFTW_ESTIMATE);
	ifft = fftw_plan_dft_1d(N, psip, psi, FFTW_BACKWARD,  FFTW_ESTIMATE);
	initialize();
	while(k<=time){
		k=k+dt2;
		evolve();
		printf("time envolved = %lf", k);
		printf("\n");
		attvexp();
		}
	
	fftw_destroy_plan(fft);
	fftw_destroy_plan(ifft);
	
	write_file(f);
	fclose(f);
	return 0;
}

double V(double x)
{
	
		return x*x;
}

void initialize()
{
	int i;
	double wr,wi,x, p;
	for(i=0; i<N; i++){
		x = x_min + x_step*(i+1);
		//V_exp[i]=cexp(-I*V(x)*dt);
		psi[i] = psi0(x);
		//+ signal for dark soliton, g<0 repulsive
		//- signal for bright soliton, g>0 atractive
		//wr = creal(psi[i]);
		//wi = cimag(psi[i]);
		//V_exp[i]=cexp(-100*I*(wr*wr+wi*wi)*dt/2.0);
		attvexp();
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

void attvexp(){
	int i;
	double wr,wi,x;
	for(i=0; i<N; i++){
		//x = x_min + x_step*(i+1);
		wr = creal(psi[i]);
		wi = cimag(psi[i]);
		V_exp[i]=cexp(-200*I*(wr*wr+wi*wi)*dt/2.0);
		}
		

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

double complex psi0(double x)
{
	return exp(-(x*x));
	return  1/cosh(x);
	
}

void normalize()
{
	double ww[N];
	double prod=0;
	double A;
	int i;

	//my normalization
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
	double wr,wi,x;
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
	
	normalize();
	attvexp();
	
	//envolve half nabla 2
	for(i=0; i<N; i++){
		psi[i] *= V_exp[i];
	}
	normalize();
	
}





