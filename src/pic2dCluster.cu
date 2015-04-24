#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <ctime>
#include <cstring>
#include <fstream>
#include <complex.h>
#include <fftw3.h>
#include <curand_kernel.h>

//Magnitudes f.

#define m_e			9.10938291e-31
#define e			1.6021e-19
#define m_e         9.10938291e-31	// Masa del Electrón
#define e           1.6021e-19		  // Carga del Electrón
#define k_b         1.3806504e-23	  // Constante de Boltzmann
#define epsilon_0	8.854187e-12	  // Permitividad eléctrica del vacío

#define max_SPe     10000        	// Limite (computacional) de Superpartículas electrónicas
#define max_SPi     10000      	// Limite (computacional) de Superpartículas iónicas
#define J_X         4096      		// Número de puntos de malla X
#define J_Y			1024			// Número de puntos de malla Y

int     le=0, li=0,kt;


//Distribución números aleatorios para las coordenadas x.
__device__  double create_Velocities_X(double fmax, double vphi, double  aleatorio, curandState *states) // función para generar distribución semi-maxwelliana de velocidades de las particulas
                                             // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
{
  double sigma=vphi;                           // sigma=vflujo=vth    ( "dispersión" de la distribución Maxweliana)
  double vmin= 0. ;                            // Rapidez mínima
  double vmax= 4.*sigma;                       // Rapidez máxima
  double v,f,f_random;

  int Idx = blockIdx.x * blockDim.x + threadIdx.x;

  while (true) {

  v=vmin+(vmax-vmin)*aleatorio; // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
  f =fmax*exp(-(1.0/M_PI)*pow(v/vphi,2));     //
  f_random = fmax*aleatorio;    // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]

  if (f_random > f)

	  aleatorio = curand_uniform(states + Idx);
  else
	  return  v;
  }
}

// Funcion de distribución para la coordenadas en y.

__device__  double create_Velocities_Y(double fmax1, double vphi1, double  aleatorio, curandState *states) // función para generar distribución semi-maxwelliana de velocidades de las particulas
                                             // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
{
  double sigma=vphi1;                           // sigma=vflujo=vth    ( "dispersión" de la distribución Maxweliana)
  double vmin= -3.*sigma;                            // Rapidez mínima
  double vmax=  3.*sigma;                       // Rapidez máxima
  double v,f,f_random;

  int Idx = blockIdx.x * blockDim.x + threadIdx.x;

  while (true) {

  v=vmin+(vmax-vmin)*aleatorio; // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
  f =fmax1*exp(-(1.0/M_PI)*pow(v/vphi1,2));     //
  f_random = fmax1*aleatorio;    // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]

  if (f_random > f)

	  aleatorio = curand_uniform(states + Idx);
  else
	  return  v;
  }
}


__global__  void distribucionVelocidadX(double *vel, double fmax, double vphi, curandState *states, int seed){

	int Idx = blockIdx.x * blockDim.x + threadIdx.x;

	seed = (unsigned int) (clock() * Idx);
	curand_init(seed, 0, 0, states + Idx);

	if (Idx < max_SPe) {
		vel[Idx] = create_Velocities_X(fmax, vphi, curand_uniform(states + Idx), states); // Distribucion_X

	}

}

__global__  void distribucionVelocidadY(double *vel1, double fmax1, double vphi1, curandState *states, int seed){

	int Idx = blockIdx.x * blockDim.x + threadIdx.x;

	seed = (unsigned int) (clock() * Idx);
	curand_init(seed, 0, 0, states + Idx);

	if (Idx < max_SPe) {
		vel1[Idx] = create_Velocities_Y(fmax1, vphi1, curand_uniform(states + Idx), states); // Distribucion_X
		//vel_1[Idx] = distrib_vel_X(fmax, vphi, curand_uniform(states + Idx), states, N); // DIstribucion_Y
	}

}




using namespace std;

int main(void) {

	double razon_masas = 1.98e5;  //m_i/m_e (plata)
	double vphi_i_0;    		  // masa Ion en la posicion 0;
	double vphi_i_1;			  // masa Ion en la posicion 1;
	double vphi_e_0;			  // masa electrones en la posición 0;
	double vphi_e_1;			  // masa electrones en la posición 1;
	double fi_Maxwell_0;		  //
	double fi_Maxwell_1;
	double fe_Maxwell_0;
	double fe_Maxwell_1;
	double	vflux_i_0 = 1e3;
	double  vflux_i_1 = 1e3;
	double	vflux_e_0 =(sqrt(razon_masas)*vflux_i_0);
	double  vflux_e_1 =(sqrt(razon_masas)*vflux_i_1);
	double  vflux_i_magnitud=sqrt(vflux_i_0*vflux_i_0+vflux_i_1*vflux_i_1); // Velocidad de flujo iónico (m/s) = sqrt(2*k_b*Te/(M_PI*m_i))
	double  vflux_e_magnitud=sqrt(vflux_e_0*vflux_e_0+vflux_e_1*vflux_e_1);
	vphi_i_0=vflux_i_0/vflux_i_magnitud;    // Velocidad térmica Iónica (X)
	vphi_e_0=vflux_e_0/vflux_i_magnitud;    // Velocidad térmica Electrónica (X)
	vphi_i_1=vflux_i_1/vflux_i_magnitud;    // Velocidad térmica Iónica (Y)
	vphi_e_1=vflux_e_1/vflux_i_magnitud;    // Velocidad térmica Electrónica (Y)
	fi_Maxwell_0=  (2./(M_PI*vphi_i_0));    // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica (X)
	fe_Maxwell_0=  (2./(M_PI*vphi_e_0));    // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
	fi_Maxwell_1=  (1./(M_PI*vphi_i_0));    // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica
	fe_Maxwell_1=  (1./(M_PI*vphi_e_0));    // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
	int     NTSPe, NTSPI, max_SPe_dt, max_SPi_dt;
	int     NTe = 1e5, NTI = 1e5;
	int     Factor_carga_e=10, Factor_carga_i=10;

	NTSPe=NTe/Factor_carga_e;
	NTSPI=NTI/Factor_carga_i;
	int Kemision=10;
	double dt=1.e-5;
    int  dt_emision=Kemision*dt;
    max_SPe_dt= round((double)NTSPe*dt_emision);
    max_SPi_dt=max_SPe_dt;



	//////////////////////////////////////////////////////
	int size = max_SPe* sizeof(double);

	//Declaración de las variables en el host.
	double *vel_e_0;
	double *vel_e_1;
	double *vel_i_0;
	double *vel_i_1;
	vel_e_0 = (double *) malloc(size);
	vel_e_1 = (double *) malloc(size);
	vel_i_0 = (double *) malloc(size);
	vel_i_1 = (double *) malloc(size);

	//Declaración de las variables en el Device
	double *vel_e_0_d;
	double *vel_e_1_d;
	double *vel_i_0_d;
	double *vel_i_1_d;

	cudaMalloc((void **) &vel_e_0_d, size);
	cudaMalloc((void **) &vel_e_1_d, size);
	cudaMalloc((void **) &vel_i_0_d, size);
	cudaMalloc((void **) &vel_i_1_d, size);

	// crear la semilla para generar el número aleatorio
	curandState *devStates;
	cudaMalloc((void **) &devStates, max_SPe * sizeof(curandState));
	int seed = time(NULL);

	//Cantidad de Hilos a correr en cada bloque.
	float blockSize = 1024;
	dim3 dimBlock(ceil(max_SPe/ blockSize), 1, 1);
	dim3 dimGrid(blockSize, 1, 1);

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	distribucionVelocidadX<<<blockSize, dimBlock>>>(vel_e_0_d, fe_Maxwell_0,  vphi_e_0, devStates, seed);
	cudaDeviceSynchronize();
	distribucionVelocidadX<<<blockSize, dimBlock>>>(vel_i_0_d, fi_Maxwell_0,  vphi_i_0, devStates, seed);
	cudaDeviceSynchronize();
	distribucionVelocidadY<<<blockSize, dimBlock>>>(vel_e_1_d, fe_Maxwell_1,  vphi_e_1, devStates, seed);
	cudaDeviceSynchronize();
	distribucionVelocidadY<<<blockSize, dimBlock>>>(vel_i_1_d, fi_Maxwell_1,  vphi_i_1, devStates, seed);
	cudaDeviceSynchronize();


	cudaMemcpy(vel_e_0, vel_e_0_d, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_i_0, vel_i_0_d, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_e_1, vel_e_1_d, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_i_1, vel_i_1_d, size, cudaMemcpyDeviceToHost);

	ofstream init;
	init.open("velocidad_X");//se escribe un archivo de salida para analizar los datos. la salida corresponde al potencial electrostatico en cada celda conocido como phi.
		for (int i = 0; i < max_SPe; i++){
			init<<vel_e_0[i]<<" "<<vel_i_0[i]<<" "<<vel_e_1[i]<<" "<<vel_i_1[i]<<"\n";
		}
		init<<endl;
	init.close();

	printf("hola");


return 0;

}
