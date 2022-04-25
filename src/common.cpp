#include "common.h"
#include <gsl/gsl_rng.h>

int N=500000;//number of clusters to be sedimented (Note: overwritten in main function)
int br=0;//global event counter

int nx=200,ny=100,nz=1000;
long double dx=100,dy=100,dz=100;
long double xo=0,yo=0,zo=-2*dz;

long double lx=nx*dx;
long double ly=ny*dy;
long double pxb=0;
long double pxe=lx;
long double pyb=0;
long double pye=ly;

string cmsfile="cms/cms";
string positionsfile="positions/positions";
string detailsfile="details/state";

long double error=0;
long double error1=0;
long double eps=60*1e-10;


long double phi=0.01;//delta phi value
rotinftype rotinfphi;

long double Rmax=34;//Maximal radius of a sphere in a system
long double h0=60000000;//initial dropping height of particles
long double epsz=60*1e-10;//particle is allowed to increase its height by maximum this value, to prevent issues when the structure is crystilline
long double epstang=1e-6;//allows nearly grazing collisions to be classified as collisions.
int particle;

long double l_pi=3.14159265359;

const gsl_rng_type * T;//random number generation global variables
gsl_rng * r;

vtype rand_orientation()//returns vectors uniformly distributed on a unit sphere
{
	vtype res;
	do
	{
		res.x()=ranfs();
		res.y()=ranfs();
		res.z()=ranfs();
	}
	while(norm(res)>1);

	res=(1./norm(res))*res;
	return res;
}

