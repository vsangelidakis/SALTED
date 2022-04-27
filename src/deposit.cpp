#include "common.h"
#include "surround.h"
#include "functions1.h"
#include "cluster.h"
#include "numerics.h"



int main(int argc,char** argv)
{	
	initrand(); // initialisation of random number generator (apart from seed)
	load_parameters("input"); // reads file with input parameters

	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW); // simulation exits if NaN etc occur

	phi=atof(argv[1]); // reads parameter dphi from command line;
	rotinfphi.cc=cos(phi); // this and next line: assign cos and sin of dphi to global variable rotinfphi 
	rotinfphi.ss=sin(phi);
	ds("phi",phi);

	int seedgsl=atoi(argv[2]);
	gsl_rng_set(r,seedgsl); // sets the random number seed

	int nspheres=atoi(argv[3]); // cluster size
	particle=atoi(argv[4]); // particles shape to be deposited
	int N=atoi(argv[5]); // number of particles to be deposited
	bool doredrop=atoi(argv[6]); // should the particle be redropped?
	string file_identifier=argv[7]; // string appended to output filenames

	if(particle==0) ds("Spheres");
	if(particle==1) ds("Spiky particles");
	if(particle==2) ds("Cubes");
	if(particle!=0&&particle!=1&&particle!=2) {ds("Not known particle shape");exit(1);}

	vector<clustertype> co;
	co.resize(N);

	for(int i=0;i<N;i++) // filling up the vector containing clusters, represented with object cp. The particles from this vector will be deposited
	{
		clustertype cp;
		if(particle==0)	cp.init_sphere();
		if(particle==1)	cp.init_axes_lines_mshuffle(nspheres);
		if(particle==2)	cp.init_cube_mshuffle(nspheres);

		long double ddx=lx/2+lx/2*ranfs(),ddy=ly/2+ly/2*ranfs();
		for(int k=0;k<cp.r.size();k++) 
		{
			cp.r[k].x()+=ddx;
			cp.r[k].y()+=ddy;
		} 
		cp.cm=cp.cmass();
		cp.randrot();

		co[i]=cp;	
	}

	surroundtype su;
	su.init_grid(); // initializations of the object containing fixed particles

	for(int i=0;i<co.size();i++) // loop in which particles are sedimented 
	{
		if(i%100==0) ds("i",i);
		clustertype c=co[i];
		if(!doredrop) // if no redropping leave particle in (rare) unstable positions 
		{
			c.settle(su); // sediments the particle (cluster)		
		}
		else // if there is a problem redropping of a cluster is allowed
		{
			bool collapsed;
			clustertype cplr=c;
			do // loop until the particle reaches a proper stable position (needed very rarely)
			{
				collapsed=c.settle(su); // sediments the particle (cluster)		
				if(collapsed) // choose new random position and orientation
				{
					c=cplr;
					long double ddx=-c.cm.x()+lx/2.+lx/2.*ranfs();
					long double ddy=-c.cm.y()+ly/2.+ly/2.*ranfs();

					for(int q=0;q<c.r.size();q++) 
					{
						c.r[q].x()+=ddx;c.r[q].y()+=ddy;
					} 
					c.cm=c.cmass();
					c.randrot();	
				}
			}
			while(collapsed);
		}
		unload(c,su); // the cluster has been sedimented, add it to the fixed particles
	}

	// output information about the simulation

	string positionsfile_str=positionsfile+file_identifier;
	output(su,positionsfile_str.c_str()); // positions

	string cmsfile_str=cmsfile+file_identifier;
	write2cm(cmass(su.r,su.R).z(),cmsfile_str.c_str()); // filling height

	int intersect=checkinters(su); // check for intersections
	string detailsfile_str=detailsfile+file_identifier;
	writedetails(co.size(),intersect,detailsfile_str.c_str()); // other information

	return 0;
}


