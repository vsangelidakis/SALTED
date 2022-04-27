#include "functions1.h"
#include "common.h"

// initalises random number generator (apart from seed)
void initrand()
{
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
}

// reads file with input parameters
void load_parameters(const char* inputfile)
{
	string s1,s2,s3,s4,s5,s6;
	ifstream myfile(inputfile);
	myfile.precision(15);

	myfile>>s1>>N;
	assert(s1=="N:");

	myfile>>s1>>nx>>s2>>ny>>s3>>nz;
	assert(s1=="nx:"&&s2=="ny:"&&s3=="nz:");

	myfile>>s1>>dx>>s2>>dy>>s3>>dz;
	assert(s1=="dx:"&&s2=="dy:"&&s3=="dz:");

	myfile>>s1>>xo>>s2>>yo>>s3>>zo;
	assert(s1=="xo:"&&s2=="yo:"&&s3=="zo:");

	lx=nx*dx;
	ly=ny*dy;
	pxb=0;
	pxe=lx;
	pyb=0;
	pye=ly;

	myfile>>s1>>phi;
	assert(s1=="phi:");

	// Output files

	myfile>>s1>>cmsfile>>s2>>positionsfile>>s4>>detailsfile;
	assert(s1=="cmsfile:"&&s2=="positonsfile:"&&s4=="detailsfile:");

	myfile.close();
}



// calculates the center of mass of the packing
vtype cmass(const vector<vtype>& r, const vector<long double>& R)
{
	int i;long double mc=0,r3;vtype res=null;
	for(i=0;i<r.size();i++)
	{
		r3=R[i]*R[i]*R[i];
		res=res+r3*(r[i]-r[0]);
		mc=mc+r3;
	}
	res=(1/mc)*res;res=res+r[0];
	return res;
}

// calculates the mass of the packing assuming that sphere with radius 30 has unit mass
long double mass(const vector<vtype>& r, const vector<long double>& R)
{
	int i;long double m=0;
	for(i=0;i<r.size();i++)
	{
		long double rr;
		rr=R[i]/30.;
		m=m+rr*rr*rr;
	}
	return m;
}

// calculates the number of overlaps in the sediment
int checkinters(surroundtype& su) 
{
	int count=0;
	for(int ix=0;ix<nx;ix++)
	for(int iy=0;iy<ny;iy++)
	for(int iz=0;iz<nz;iz++)
	for(int ii=0;ii<su.box[ix][iy][iz].size();ii++)
	{
		int i=su.box[ix][iy][iz][ii];
		vtype r=su.r[i];
		long double R=su.R[i];
		int ixl=floor((r.x()-R-Rmax-xo)/dx);
		int ixu=floor((r.x()+R+Rmax-xo)/dx);

		int iyl=floor((r.y()-R-Rmax-yo)/dy);
		int iyu=floor((r.y()+R+Rmax-yo)/dy);

		int izl=floor((r.z()-R-Rmax-zo)/dz);if(izl<0) izl=0;
		int izu=floor((r.z()+R+Rmax-zo)/dz);if(izu>=nz) izu=nz-1;

		for(int ixnt=ixl;ixnt<=ixu;ixnt++)
		for(int iynt=iyl;iynt<=iyu;iynt++)
		for(int iznt=izl;iznt<=izu;iznt++)
		{
			int ixn=imod(ixnt,nx);
			int iyn=imod(iynt,ny);
			int izn=iznt; 
			for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
			{
				int j=su.box[ixn][iyn][izn][jj];
				vtype rj=vmod(su.r[j],r,lx,ly);
				long double Rj=su.R[j];
				if(i!=j) if(norm(r-rj)-R-Rj<-60*1e-8) {ds(i,j,norm(r-rj)-R-Rj);assert(0);count++;}
			}
		}
	}
	return count/2;
}


// calculates the number of overlaps between the cluster spheres and the sediment
bool checkinters(const clustertype& c,const surroundtype& s)
{
	for(int i=0;i<c.r.size();i++)
	{
		vtype rp=c.r[i];
		long double R=c.R[i];
		int ixl=floor((rp.x()-(R+Rmax)-xo)/dx);
		int ixu=floor((rp.x()+(R+Rmax)-xo)/dx);

		int iyl=floor((rp.y()-(R+Rmax)-yo)/dy);
		int iyu=floor((rp.y()+(R+Rmax)-yo)/dy);

		int izl=floor((rp.z()-(R+Rmax)-zo)/dz);if(izl<0) izl=0;
		int izu=floor((rp.z()+(R+Rmax)-zo)/dz);if(izu>=nz) izu=nz-1;

		for(int ixnt=ixl;ixnt<=ixu;ixnt++)
		for(int iynt=iyl;iynt<=iyu;iynt++)
		for(int iznt=izl;iznt<=izu;iznt++)
		{
			int ixn=imod(ixnt,nx);
			int iyn=imod(iynt,ny);
			int izn=iznt; 
			for(int jj=0;jj<s.box[ixn][iyn][izn].size();jj++)
			{
				int j=s.box[ixn][iyn][izn][jj];
				{
					vtype rj=vmod(s.r[j],rp,lx,ly);
					if (norm(c.r[i]-rj)<c.R[i]+s.R[j]-60.*1e-6) {ds("Intersection: ",i,j,norm(c.r[i]-rj)-c.R[i]-s.R[j]);return true;}
				}
			}
		}
	}
	return false;
}


// finds number of contacts that the cluster has with its surroundings
void find_cluster_contacts(long double eps_contact,const clustertype& c,const surroundtype& su,vector<vtype>& rps,vector<vtype>& rds,vector<int>& p,vector<int>& d)
{
	for(int i=0;i<c.r.size();i++)
	{
		long double lcd=c.R[i]+Rmax;
		vtype r=c.r[i];
		int ixl=floor((r.x()-lcd-xo)/dx);
		int ixu=floor((r.x()+lcd-xo)/dx);

		int iyl=floor((r.y()-lcd-yo)/dy);
		int iyu=floor((r.y()+lcd-yo)/dy);

		int izl=floor((r.z()-lcd-zo)/dz);if(izl<0) izl=0;
		int izu=floor((r.z()+lcd-zo)/dz);if(izu>=nz) izu=nz-1;

		for(int ixnt=ixl;ixnt<=ixu;ixnt++)
		for(int iynt=iyl;iynt<=iyu;iynt++)
		for(int iznt=izl;iznt<=izu;iznt++)
		{
			int ixn=imod(ixnt,nx);
			int iyn=imod(iynt,ny);
			int izn=iznt; 
			for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
			{
				int j=su.box[ixn][iyn][izn][jj];
				{
					vtype rj=vmod(su.r[j],r,lx,ly);
					if (norm(r-rj)<c.R[i]+su.R[j]+eps_contact) 
					{
						rps.push_back(rj);	
						rds.push_back(r);	
						p.push_back(j);
						d.push_back(i);
					}
				}
			}
		}
	}
}


// adds the cluster c to the fixed particles
void unload(clustertype& c,surroundtype& su)
{
	vtype v;
	int susize=su.r.size();
	for(int i=0;i<c.r.size();i++) 
	{
		int ii=susize+i;
		unloadsphere(vmod(c.r[i],lx,ly),ii,c.R[i],su);
	}
}


void unloadsphere(const vtype& v,int ii,long double R,surroundtype& su)
{
	su.r.push_back(v);
	su.R.push_back(R);

	int ix=floor(v.x()-xo)/dx;
	int iy=floor(v.y()-yo)/dy;
	int iz=floor(v.z()-zo)/dz;
	assert(0<=ix&&ix<nx&&0<=iy&&iy<ny&&0<=iz&&iz<nz);

	su.box[ix][iy][iz].push_back(ii);
	if(iz>su.izmax[ix][iy]) su.izmax[ix][iy]=iz;
}

//outputs the fixed particles to a file
void output(surroundtype& su,string fname) 
{
	int i;ofstream myfile;
	myfile.open (fname.c_str());
	myfile.precision(15);
	for(i=0;i<su.r.size();i++) {myfile<<su.r[i].x()<<" "<<su.r[i].y()<<" "<<su.r[i].z()<<" "<<su.R[i]<<endl;}
	myfile.close();
}


//outputs a set of clusters to a single file
void output(vector<clustertype>& cd,string fname)
{
	int i;ofstream myfile;
	myfile.open(fname.c_str());
	myfile.precision(15);

	for(i=0;i<cd.size();i++) 
	for(int j=0;j<cd[i].r.size();j++) 
	myfile<<cd[i].r[j].x()<<" "<<cd[i].r[j].y()<<" "<<cd[i].r[j].z()<<" "<<cd[i].R[j]<<endl;

	//cm
	//for(i=0;i<cd.size();i++) 
	//myfile<<cd[i].cm<<endl;

	myfile.close();
}

//outputs to a file the number of clusters sedimented (input varibale size) and the number of intersections in the sediment (also input parameter)
void writedetails(int size,int intersect,string name)
{
	ofstream myfile;
	time_t seconds=time(NULL);
	myfile.open(name.c_str(),ios_base::app);
	myfile<<"  Number of complex particles:  "<<size<<"  Number of intersections  "<<intersect<<endl;
	myfile.close();
}


//outputs to a file the filling hegiht of the sediment
void write2cm(long double z,string name)
{
	ofstream myfile;
	myfile.open(name.c_str(),ios::app);
	myfile.precision(15);
	myfile<<2*z<<endl;
	myfile.close();
}
