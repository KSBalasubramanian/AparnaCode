#include <iostream>
#include <boost/random.hpp>
#include <time.h>
#include <math.h>

using namespace std;
using namespace boost;


const double pi = 3.14159265358979323846;

int main()
 { 
//Initialize random class etc.
	 mt19937    generator;    // Mersenne Twister
	 generator.seed(time(0)); // seed
	 normal_distribution<>norm_dist;
	 variate_generator<mt19937&,normal_distribution<> >norm_rand(generator,norm_dist);//distribution & generator
//done with it, Use norm_rand()

int ensemble=10000;
int runtime=10;
double dt=0.001;
double sqdt=sqrt(dt);
//Numerical run constants

double theta0=pi/4,phi0=pi/4;
double Dr=1;
double theta=theta0;
double phi=phi0;
//Physics constants

double* nx= new double[int(runtime/dt)]();
double* ny= new double[int(runtime/dt)]();
double* nz= new double[int(runtime/dt)](); 


	
//double nx[floor(int(floor(runtime/dt)))]={},ny[int(runtime/dt)]={},nz[int(floor(runtime/dt))]={};

nx[0]=sin(theta0)*cos(phi0);
ny[0] = sin(theta0)*sin(phi0); 
nz[0]=cos(theta0);


//cout<<nx[0]<<endl<<nz[0]<<endl<<sum[0];

    for (int i = 0; i < ensemble; i++)
	{

		theta=theta0;
		phi=phi0;

		for(int j=0;j<(int(floor(runtime/dt)));j++)
		{
			

			theta+=dt*Dr*Dr/tan(theta)+sqrt(2.0)*Dr*sqdt*norm_rand();
			phi+=sqrt(2.0)*Dr*sqdt*norm_rand()/sin(theta);
			nx[j]+=sin(theta)*cos(phi);
			ny[j]+=sin(theta)*sin(phi);
			nz[j]+=cos(theta);
			
		}
	        
 }

for(int j=0;j<(int(floor(runtime/dt)));j++)
{
	cout<<dt*j<<" " <<(nx[j]*nx[0])/(ensemble*ensemble)<<" "<<(ny[j]*ny[0])/(ensemble*ensemble)<<" "<<(nz[j]*nz[0])/(ensemble*ensemble)<<" "<<sum[j]/ensemble<<endl;

}

    return 0;
}
