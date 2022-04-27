#include "cluster.h"
#include "functions1.h"
#include "mod.h"

char name[100];

vtype clustertype::cmass()
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


void clustertype::erase()
{
	r.erase(r.begin(),r.end());
	R.erase(R.begin(),R.end());
}

void clustertype::init_sphere()
{
	contacts=0;
	stable=false;
	p1=-1;d1=-1;
	p2=-1;d2=-1;
	p3=-1;d3=-1;
	tp=-1;td=-1;
	f=0;

	// addition of a sphere
	r.clear();R.clear();
	vtype v(0,0,h0); //FIXME: Rename this, as v() is already a function in functions1.cpp. Not a severe problem, as the v() function is not used here and the v variable here is local, just to avoid future confusions. Or even better, we could give the v() function a more intuitive name, e.g. cross()
	r.push_back(v);R.push_back(30.);
	// end of sphere addition

	cm=cmass();
}


void clustertype::init_axes_lines_mshuffle(int nsph_per_side)
{
	r.clear();R.clear();
	vtype n;
	long double eps_sh=60*1e-4;
	long double Rp=30.-eps_sh;
	long double le=2*Rp*(nsph_per_side-1);
	long double D=2*30.;

	n=1*ex;
	for(long double l=0;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=-1*ex;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=1*ey;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=-1*ey;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=1*ez;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=-1*ez;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	contacts=0;
	stable=false;
	p1=-1;d1=-1;
	p2=-1;d2=-1;
	p3=-1;d3=-1;
	tp=-1;td=-1;
	f=0;
	cm=cmass();
}

void clustertype::init_cube_mshuffle(int nsph_per_side)
{
	r.clear();R.clear();
	vtype n;
	long double eps_sh=60*1e-4;
	long double Rp=30.-eps_sh;
	long double le=2*30.*(nsph_per_side-1);
	long double D=2*30.;

	// xy plane
	for(int i=0;i<nsph_per_side;i++)
	for(int j=0;j<nsph_per_side;j++)
	{vtype v=D*(i*ex+j*ey);v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// xy plane shifted
	for(int i=0;i<nsph_per_side;i++)
	for(int j=0;j<nsph_per_side;j++)
	{vtype v=D*(i*ex+j*ey);v.z()=v.z()+h0;v.z()=v.z()+(nsph_per_side-1)*D;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}


	// xz plane
	for(int i=0;i<nsph_per_side;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ex+j*ez);v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// xz plane shifted
	for(int i=0;i<nsph_per_side;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ex+j*ez);v.z()=v.z()+h0;v.y()=v.y()+(nsph_per_side-1)*D;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}


	// yz plane
	for(int i=1;i<nsph_per_side-1;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ey+j*ez);v.z()=v.z()+h0;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// yz plane shifted
	for(int i=1;i<nsph_per_side-1;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ey+j*ez);v.z()=v.z()+h0;v.x()=v.x()+(nsph_per_side-1)*D;vtype dv(eps_sh*ranfs(),eps_sh*ranfs(),eps_sh*ranfs());r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}


	contacts=0;
	stable=false;
	p1=-1;d1=-1;
	p2=-1;d2=-1;
	p3=-1;d3=-1;
	tp=-1;td=-1;
	f=0;
	cm=cmass();
}

bool clustertype::settle(surroundtype& su) // sediments a particle
{
	bool collapsed=false;
	int countev=1;
	vtype cmev(cm);

	long double dtotal=0;
	while(!stable)
	{
		br++;

		countev++;
		vtype cmp=cm;

		switch(contacts) // propagates cluster by one event
		{
			case 0: fall(su);break; // propagate when cluster is falling
			case 1: rotate1(su);break; // propagate when cluster is rotating on one contact
			case 2: rotate2(su);break; // propagate when cluster is rotating on two contacts
		}

		long double dcm1=norm(cm-cmp);

		long double eps_collapse=60*1e-10;
		long double eps_contact=60*1e-8;
		if(dcm1<eps_collapse&&particle!=0) // for spheres no need to detect collapse 
		{
			// collapse detected, going through collapse

			vector<vtype> rps;
			vector<vtype> rds;
			vector<int> p;
			vector<int> d;
			find_cluster_contacts(eps_contact,*this,su,rps,rds,p,d); // finds contacts by proximity

			// check if there are contacts where one cluster particle touches two fixed particles (duplicate)
			bool duplicate_exists=false;
			for(int i=0;i<p.size()-1;i++)
			for(int j=i+1;j<p.size();j++)
			if(p[i]==p[j]) duplicate_exists=true;

			if(!duplicate_exists) // if there are no duplicates find steepest descent trajectory
			{
				rotattempt att;
				find_steppest_descent(rps,rds,p,d,cm,att); // find steepest descent with closed contacts considered
				// update state
				contacts=att.contacts;
				stable=att.stable;
				p1=att.p1;
				d1=att.d1;
				p2=att.p2;
				d2=att.d2;
				p3=att.p3;
				d3=att.d3;
				td=att.td;
				tp=att.tp;
				f=att.f;
				// particle stays where it is so no need to change cm and r[i]
			}
			else
			{
				// long double contact exists, cannot proceed
				stable=true;
				collapsed=true; // considered as collapsed, in order to be redropped, if chosen
			}

		}

		if(countev%1000000==0) // parabolic collapse fixing
		{
			if(norm(cmev-cm)<60) 
			{
				stable=true;
				collapsed=true;
				ds("Moved less than a particle diameter in 1000000 events");
			}
			else cmev=cm;
		}

		/*****************  Important: periodic boundary conditions treatment  *****************/
		periodic(lx,ly);
	}
	return collapsed;
}

void clustertype::periodic(long double lx,long double ly) // shifts cluster back to the periodic cell
{
	vtype v;
	v.z()=0;

	if(cm.x()<pxb) v.x()=lx;
	if(cm.x()>pxb&&cm.x()<pxe) v.x()=0;
	if(cm.x()>pxe) v.x()=-lx;

	if(cm.y()<pyb) v.y()=ly;
	if(cm.y()>pyb&&cm.y()<pye) v.y()=0;
	if(cm.y()>pye) v.y()=-ly;

	movecluster(v);
}

void clustertype::movecluster(const vtype& v) // moves the whole cluster by vector v
{
	for(int i=0;i<r.size();i++) r[i]=r[i]+v;
	cm=cm+v;
}


void clustertype::fall(surroundtype& su)
{
	vector<fallattempt> atvec;fallattempt att;

	bfall(su,atvec); // add events for collisions with fixed spheres
	gfall(atvec); // add events for collisions with the ground

	if(atvec.size()==0) {ds("Error: fall: no events");exit(1);}//must be at least one future event
	att=*max_element(atvec.begin(),atvec.end()); // choose the highest future event

	step0(att,su); // move the cluster
}

void clustertype::bfall(surroundtype& su,vector<fallattempt>& atvec) // possible collisons with fixed sphere
{
	int j,k;long double dlt;fallattempt att;
	bool contact;long double zfc;
	vtype rp;
	for(k=0;k<r.size();k++) // loop through all cluster particles
	{
		zfc=h0;contact=false;
		int ixl,ixr,iyl,iyr,izu,izd;findindex(r[k],R[k],ixl,ixr,iyl,iyr,izu);
		for(int it=ixl;it<=ixr;it++) // loop through boxes containing candidates
		for(int nt=iyl;nt<=iyr;nt++)
		{
			int i=imod(it,nx);
			int n=imod(nt,ny);
			for(int m=su.izmax[i][n];m>=0&&(!contact||(contact&&m>=floor((zfc-R[k]-Rmax-zo)/dz)));m--)
			for(int s=0;s<su.box[i][n][m].size();s++)
			{
				j=su.box[i][n][m][s];
				rp=vmod(su.r[j],cm,lx,ly);
				if(k!=td||j!=tp) // to avoid recontact when conacts was lost by rolling to horizontal
				if(pcontact0(r[k],R[k],rp,su.R[j],dlt)) // test for collision by translation 
				if(r[k].z()-dlt>rp.z()) // gradient check
				{
					if(dlt>0&&!contact) {contact=true;zfc=r[k].z()-dlt;}
					if(dlt>0&&r[k].z()-dlt>zfc) zfc=r[k].z()-dlt;

					att.contacts=1;
					att.stable=false;
					att.p1=j;
					att.d1=k;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=0;
					att.cm=cm;
					att.cm.z()=att.cm.z()-dlt;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}
				}
			}
		}
	}
}

void clustertype::gfall(vector<fallattempt>& atvec) // finds events for collision with the ground
{
	int k;long double dlt;fallattempt att;
	for(k=0;k<r.size();k++) 
	{
		dlt=r[k].z()-R[k];
		att.stable=true;
		att.contacts=1;
		att.d1=k;
		att.p1=-2;
		att.d2=-1;
		att.p2=-1;
		att.d3=-1;
		att.p3=-1;
		att.td=-1;
		att.tp=-1;
		att.f=0;
		att.cm=cm;
		att.cm.z()=att.cm.z()-dlt;
		if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
	}
}

void clustertype::step0(fallattempt att,surroundtype& su) // updates the clusters state
{
	d1=att.d1;
	p1=att.p1;
	d2=att.d2;
	p2=att.p2;
	d3=att.d3;
	p3=att.p3;
	td=att.td;
	tp=att.tp;
	f=att.f;
	contacts=att.contacts;
	stable=att.stable;
	translate(cm.z()-att.cm.z(),d1,p1,su);
}

void clustertype::translate(long double dlt,int d1,int p1,surroundtype& su) // translates the cluster
{
	int i;
	for(i=0;i<r.size();i++) 
	{
		r[i].z()-=dlt;
	}
	cm.z()-=dlt;

	if(p1>=0&&d1>=0)
	{
		// correction due to loss of precision due to use of dlt
		vtype rp1=su.r[p1];
		vtype rd1l=vmod(r[d1],rp1,lx,ly);

		long double dl=su.R[p1]+R[d1];
		long double ds2=pow(rd1l.x()-rp1.x(),2)+pow(rd1l.y()-rp1.y(),2);
		long double underroot=1-ds2/(dl*dl);
		if(underroot>eps) 
		{
			long double ovrl=norm(rp1-rd1l)-dl;
			r[d1]=r[d1]-ovrl*1/sqrt(underroot)*ez;
			rd1l=vmod(r[d1],rp1,lx,ly);
		}
	}
}

void clustertype::rotate1(surroundtype& su) // moves the cluster by one event when cluster moves on one particle, structure is analogous to function fall
{
	vector<rotattempt> atvec;rotattempt att;
	vtype rp_temp=vmod(su.r[p1],cm,lx,ly);

	// correction of the case where the center of mass is above the center of fixed particle
	bool cm_on_top_of_fixed_particle=fabs(cm.x()-rp_temp.x())<60*1e-6&&fabs(cm.y()-rp_temp.y())<60*1e-6;
	vtype sh;

	if(cm_on_top_of_fixed_particle) // moves p1 particle slightly
	{	
		sh.z()=0;
		do 
		{
			sh.x()=ranfs();
			sh.y()=ranfs();
		}
		while(norm(sh)<0.5&&norm(sh)>1);
		sh=eps*sh;
		su.r[p1]=su.r[p1]+sh;
	}
	int p1old=p1;

	vtype rp=vmod(su.r[p1],cm,lx,ly);

	steeperdphi(rp,atvec); // dphi event
	brot1(su,atvec); // add events for collisions with fixed spheres
	hrot(rp,atvec); // add event for rotation to horizontal 
	bottomrot1(rp,atvec); // add event for rotation to hanging
	grot1(rp,atvec); // add event for contact with the ground
	if(atvec.size()==0) {ds("Error: rotate1: no events");exit(1);}

	att=*max_element(atvec.begin(),atvec.end());
	step1(att,rp); // moves cluster

	if(cm_on_top_of_fixed_particle) su.r[p1old]=su.r[p1old]-sh;
}

void clustertype::steeperdphi(const vtype& rp,vector<rotattempt>& atvec) // deals with "frustrated" situation when cluster has one contact
{
	if(r[d1].z()<rp.z()) // if cluster particle is below fixed one then cluster is frustrated situation
	{		
		rotinftype rotinf=rotinfphi;
		vtype cmr=cm-rp;
		vtype n=v(ez,cmr);n=(1/norm(n))*n;
		rotinf.n=n;
		
		vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n); // rotated center of mass
		vtype rdr=rot(r[d1],rp,rotinf.cc,rotinf.ss,rotinf.n); // rotated cluster particle d1

		if(rightside(cm-rp,cmc-rp,rotinf.n)) // checks whether the cluster is on the same side of the plane with normal ez x n, please see paper
		{
			if(rdr.z()<rp.z()) // check whether after rotation by dphi angle cluster sphere is still below fixed sphere
			{
				rotattempt att;
				att.contacts=0;
				att.stable=false;
				att.d1=-1;
				att.p1=-1;
				att.d2=-1;
				att.p2=-1;
				att.d3=-1;
				att.p3=-1;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc;
				att.rotinf=rotinf;
				atvec.push_back(att);
				assert(cmc.z()<cm.z()+epsz);
			}
		}
	}
}


void clustertype::brot1(surroundtype& su,vector<rotattempt>& atvec) // search for collision with fixed spheres while rotating
{
	rotattempt att;
	gradrotattempt gatt;
	int j,k;
	vtype rp=vmod(su.r[p1],cm,lx,ly);

//find
	long double rmax=0;
	long double lmax=0;
	for(int k=0;k<r.size();k++) // find the rotation angle between test positions of the cluster such that no events are missed
	{
		long double rt=distance_from_axis(r[k],rp,cm);
		if(rt>rmax) 
		{
			rmax=rt;
			lmax=2*(Rmax+R[k]);
		}
	}
	long double dephi;

	if(rmax<lmax) dephi=l_pi/3.;
	else dephi=2*asin(lmax/(2*rmax)); // angle increment

	long double phic=0;
	long double phidet=l_pi;
	do//do while loop until first contact has been detected without doubt
	{
		for(int k=0;k<r.size();k++) // loop through cluster particles
		{
			vtype point=rotpoint(r[k],rp,cm,phic); // position of the cluster particle according to current test angle
			long double axisr=distance_from_axis(r[k],rp,cm);
			long double l=2*axisr*sin(dephi/2);
			long double max_cillinder=R[k]+Rmax+axisr-sqrt(axisr*axisr-l*l/4.);
			long double radius=sqrt(l*l/4.+max_cillinder*max_cillinder); // radius of the search volume

			vtype rref=point; // center of spheres giving the neighbours

			// testing collision with candidate spheres 
			int ixl=floor((rref.x()-radius-xo)/dx);
			int ixu=floor((rref.x()+radius-xo)/dx);

			int iyl=floor((rref.y()-radius-yo)/dy);
			int iyu=floor((rref.y()+radius-yo)/dy);

			int izl=floor((rref.z()-radius-zo)/dz);if(izl<0) izl=0;
			int izu=floor((rref.z()+radius-zo)/dz);if(izu>=nz) izu=nz-1;

			for(int ixnt=ixl;ixnt<=ixu;ixnt++)
			for(int iynt=iyl;iynt<=iyu;iynt++)
			for(int iznt=izl;iznt<=izu;iznt++)
			{
				int ixn=imod(ixnt,nx);
				int iyn=imod(iynt,ny);
				int izn=iznt; 

				for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
				{
					j=su.box[ixn][iyn][izn][jj]; // index of the candidate sphere
					vtype rp2=vmod(su.r[j],cm,lx,ly); // periodic conditions applied to candidate sphere
					{
						rotinftype rotinf1,rotinf2,rotinf; // information about rotation
						if(j!=p1) // no need to test for contact 1
						if(!f||k!=td||j!=tp) // to avoid reconctact after grazing loss of contact
						if(pcontact1(rp,cm,r[k],R[k],rp2,su.R[j],rotinf1,rotinf2)) // test for collision
						{
							vector<rotinftype> rotv;
							rotv.erase(rotv.begin(),rotv.end());
							rotv.push_back(rotinf1);rotv.push_back(rotinf2);

							for(int ir=0;ir<2;ir++) // consider both solutions for collision
							{
								rotinf=rotv[ir];

								// rotated relevant quantities
								vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rdc=rot(r[d1],rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rcc=rot(r[k],rp,rotinf.cc,rotinf.ss,rotinf.n);
								if(rightside(cm-rp,cmc-rp,rotinf.n)&&g(rp,cmc,rcc,rp2)<0&&cmc.z()<cm.z()+epsz) // final check if collision is a valid future event
								{
									vector<vtype> rps;
									rps.push_back(rp);rps.push_back(rp2);
									vector<vtype> rds;
									rds.push_back(rdc);rds.push_back(rcc);
									vector<int> p;
									p.push_back(p1);p.push_back(j);
									vector<int> d;
									d.push_back(d1);d.push_back(k);
									find_steppest_descent(rps,rds,p,d,cmc,att); // consider all contacts and find steepest descent trajectory
									att.rotinf=rotinf; 
									atvec.push_back(att); // add event

									// if by numerical error |cos|>1
									if(att.rotinf.cc<-1) att.rotinf.cc=att.rotinf.cc+1e-12;
									if(att.rotinf.cc>1) att.rotinf.cc=att.rotinf.cc-1e-12;
									long double phidetc=acos(att.rotinf.cc); // rotation angle for the current collision
									if(phidetc<phidet) phidet=phidetc; // update smallest colision angle
								}//right side if
							}//+- solution for loop
						}
					}
				}
			}		
		}
		phic+=dephi; // increase the angle for the cluster test position
	}
	while((phic-dephi)<phidet); // checks whether the current test angle is too large compared to current collision angle, for earlier collision to be possible  
}

void clustertype::hrot(const vtype& rp,vector<rotattempt>& atvec) // adds event for case in Fig 4b of the paper
{
	rotattempt att;
	rotinftype rotinf;
	tohorizontal(cm,r[d1],rp,rotinf);
	att.contacts=0;
	att.stable=false;
	att.p1=-1;
	att.d1=-1;
	att.p2=-1;
	att.d2=-1;
	att.p3=-1;
	att.d3=-1;
	att.td=d1;
	att.tp=p1;
	att.f=0;

	att.cm=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	att.rotinf=rotinf;

	vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	if(rightside(cm-rp,cmc-rp,rotinf.n)&&att.cm.z()<cm.z()+epsz) atvec.push_back(att);
}

void clustertype::grot1(const vtype& rp,vector<rotattempt>& atvec) // contact with the ground
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	for(int i=0;i<r.size();i++)
	{
		if(toground1(rp,cm,r[i],R[i],rotinf1,rotinf2)) // finds position (if any) where cluster sphere i makes contact with the ground
		{
			vtype cmc1=rot(cm,rp,rotinf1.cc,rotinf1.ss,rotinf1.n);
			vtype cmc2=rot(cm,rp,rotinf2.cc,rotinf2.ss,rotinf2.n);
			if(rightside(cm-rp,cmc1-rp,rotinf1.n))
			{
				att.contacts=2;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=i;
				att.p2=-2;
				att.d3=-1;
				att.p3=-1;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}

			if(rightside(cm-rp,cmc2-rp,rotinf2.n))
			{
				att.contacts=2;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=i;
				att.p2=-2;
				att.d3=-1;
				att.p3=-1;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}
		}
	}
}

void clustertype::bottomrot1(const vtype& rp,vector<rotattempt>& atvec) // event where clusters center of mass is below the center of fixed particle
{
	rotattempt att;
	rotinftype rotinf;
	tobelow1(cm,rp,rotinf); // finds amount of rotation such that the center of mass is below the fixed particle

	att.contacts=30;
	att.stable=true;
	att.p1=-1;
	att.d1=-1;
	att.p2=-1;
	att.d2=-1;
	att.p3=-1;
	att.d3=-1;
	att.cm=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	att.rotinf=rotinf;
	atvec.push_back(att);
	assert(att.cm.z()<=cm.z()+epsz);
}

void clustertype::step1(const rotattempt& att,const vtype& rp) // updates the cluster state to the next state
{
	contacts=att.contacts;
	stable=att.stable;
	p1=att.p1;
	d1=att.d1;
	p2=att.p2;
	d2=att.d2;
	p3=att.p3;
	d3=att.d3;
	td=att.td;
	tp=att.tp;
	f=att.f;
	cm=rot(cm,rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
	for(int i=0;i<r.size();i++) r[i]=rot(r[i],rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
}

void clustertype::rotate2(surroundtype& su) // propagates cluster by one event when cluster has two contacts
{
	vector<rotattempt> atvec;rotattempt att;

	vtype rp1_temp=vmod(su.r[p1],cm,lx,ly);
	vtype rp2_temp=vmod(su.r[p2],cm,lx,ly);

	vtype n=rp2_temp-rp1_temp;n=(1/norm(n))*n;
	vtype cmr=cm-rp1_temp;cmr=(1/norm(cmr))*cmr;
	vtype vcmn=v(n,cmr);vcmn=(1/norm(vcmn))*vcmn;

	vtype sh;
	bool is_middle=fabs(s(ez,vcmn))<60*1e-9;
	if(is_middle) 
	{
		long double sh_val;
		do
		{
			sh_val=ranfs();
		}
		while(fabs(sh_val)<0.5);	
		sh=60*1e-9*sh_val*vcmn;
		su.r[p1]=su.r[p1]+sh; // moves fixed particle by small amount when cm is exactly on top of rotation axis
	}
	int p1old=p1;

	vtype rp1=vmod(su.r[p1],cm,lx,ly);
	vtype rp2=vmod(su.r[p2],cm,lx,ly);
	brot2(su,atvec); // add events for collisions with fixed spheres
	disconnect(rp1,rp2,atvec,su); // add events for the case Fig6a of the paper
	steeperdphi2(rp1,rp2,atvec); // add event for "frustrated" situation
	hrot2(rp1,rp2,atvec); // add events for the case Fig6b of the paper
	bottomrot2(rp1,rp2,atvec); // add events for the case Fig6c of the paper
	grot2(rp1,rp2,atvec); // add events for contact with the ground

	if(atvec.size()==0) {ds("Error: rotate2: no events");exit(1);}

	att=*max_element(atvec.begin(),atvec.end()); // choose earliest event
	step2(att,su,rp1); // propagate cluster by chosen event

	if(is_middle) su.r[p1old]=su.r[p1old]-sh;
}

void clustertype::steeperdphi2(const vtype& rp1,const vtype& rp2,vector<rotattempt>& atvec) // adds event for frustrated situation (if any) with two contacts 
{
	// it is guarantied that rp1 corresponds to p1 and rp2 to p2, see function rotate2()
	vector<vtype> rp;rp.push_back(rp1);rp.push_back(rp2);
	vector<vtype> rd;rd.push_back(r[d1]);rd.push_back(r[d2]);
	vector<int> p;p.push_back(p1);p.push_back(p2);
	vector<int> d;d.push_back(d1);d.push_back(d2);

	rotattempt att;
	find_steppest_descent(rp,rd,p,d,cm,att); // find steepest descent trajectory

	if(att.contacts<contacts) // is it possible to move with less contacts?
	{
		vtype cmr=cm-rp1;vtype n=rp2-rp1;n=(1/norm(n))*n;
		if(s(ez,v(n,cm-rp1))>0) n=-1*n;

		rotinftype rotinf=rotinfphi;
		rotinf.n=n;

		vtype cmc=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n); // rotate center of mass of the cluster by small angle

		if(rightside(cm-rp1,cmc-rp1,rotinf.n))
		{
			// vectors rp and rd
			vector<vtype> rp;rp.push_back(rp1);rp.push_back(rp2);

			vtype rd1r=rot(r[d1],rp1,rotinf.cc,rotinf.ss,rotinf.n);
			vtype rd2r=rot(r[d2],rp1,rotinf.cc,rotinf.ss,rotinf.n);
			vector<vtype> rd;rd.push_back(rd1r);rd.push_back(rd2r);
			
			// indices of contacting particles
			vector<int> p;p.push_back(p1);p.push_back(p2);
			vector<int> d;d.push_back(d1);d.push_back(d2);

			rotattempt att;
			find_steppest_descent(rp,rd,p,d,cmc,att); // find steepest descent after the cluster has been rotated
			
			if(att.contacts<contacts) 
			{
				att.rotinf=rotinf;
				atvec.push_back(att);
				assert(att.cm.z()<=cm.z());
			}
			else//consistency check that rolling on two has the same contacts
			{
				assert(att.stable==stable);
				assert(att.contacts==contacts);
				assert((att.p1==p1&&att.d1==d1)||(att.p2==p1&&att.d2==d1));
			}
		}
	}
	else
	{
		// check of consistency between search for steepest descent in previous and current step
		assert(att.stable==stable);
		assert(att.contacts==contacts);
		assert((att.p1==p1&&att.d1==d1)||(att.p2==p1&&att.d2==d1));
	}
}

void clustertype::brot2(surroundtype& su,vector<rotattempt>& atvec) // find contacts with fixed spheres while move on two contacts
{
	rotattempt att;
	gradrotattempt gatt;
	rotinftype rotinf1,rotinf2,rotinf;

	vtype rp=vmod(su.r[p1],cm,lx,ly);
	vtype rp2=vmod(su.r[p2],cm,lx,ly);

	long double rmax=0;
	long double lmax=0;
	for(int k=0;k<r.size();k++)
	{
		long double rt=distance_from_axis2(r[k],rp,rp2);
		if(rt>rmax) 
		{
			rmax=rt;
			lmax=2*(Rmax+R[k]);
		}
	}
	long double dephi;

	if(rmax<lmax) dephi=l_pi/3.;
	else dephi=2*asin(lmax/(2*rmax));


	long double phic=0;
	long double phidet=l_pi;
	do
	{
		for(int k=0;k<r.size();k++)
		{
			vtype point=rotpoint2(r[k],rp,rp2,cm,phic);
			long double axisr=distance_from_axis2(r[k],rp,rp2);
			long double l=2*axisr*sin(dephi/2);
			long double max_cillinder=R[k]+Rmax+axisr-sqrt(axisr*axisr-l*l/4.);
			long double radius=sqrt(l*l/4.+max_cillinder*max_cillinder);

			vtype rref=point; // center of spheres giving the neighbours
			int ixl=floor((rref.x()-radius-xo)/dx);
			int ixu=floor((rref.x()+radius-xo)/dx);

			int iyl=floor((rref.y()-radius-yo)/dy);
			int iyu=floor((rref.y()+radius-yo)/dy);

			int izl=floor((rref.z()-radius-zo)/dz);if(izl<0) izl=0;
			int izu=floor((rref.z()+radius-zo)/dz);if(izu>=nz) izu=nz-1;

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
					vtype rp3=vmod(su.r[j],cm,lx,ly);
					if(j!=p1&&j!=p2)
					{
						if(pcontact2(rp,rp2,cm,r[k],R[k],rp3,su.R[j],rotinf1,rotinf2))
						{
//							vector<gradrotattempt> gradevent;
//							gradevent.erase(gradevent.begin(),gradevent.end());

							vector<rotinftype> rotv;
							rotv.erase(rotv.begin(),rotv.end());
							rotv.push_back(rotinf1);rotv.push_back(rotinf2);

							for(int ir=0;ir<2;ir++)
							{
								rotinf=rotv[ir];
								vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rcr=rot(r[k],rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rd1c=rot(r[d1],rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rd2c=rot(r[d2],rp,rotinf.cc,rotinf.ss,rotinf.n);
								if(rightside(cm-rp,cmc-rp,rotinf.n)&&g2(rp,rp2,cmc,rcr,rp3)<0&&cmc.z()<cm.z()+epsz)
								{
									vector<vtype> rps;
									rps.push_back(rp);rps.push_back(rp2);rps.push_back(rp3);
									vector<vtype> rds;
									rds.push_back(rd1c);rds.push_back(rd2c);rds.push_back(rcr);
									vector<int> p;
									p.push_back(p1);p.push_back(p2);p.push_back(j);
									vector<int> d;
									d.push_back(d1);d.push_back(d2);d.push_back(k);
									find_steppest_descent(rps,rds,p,d,cmc,att);
									att.rotinf=rotinf;
									atvec.push_back(att);
									if(att.rotinf.cc<-1) att.rotinf.cc=att.rotinf.cc+1e-12;
	if(att.rotinf.cc>1) att.rotinf.cc=att.rotinf.cc-1e-12;
									long double phidetc=acos(att.rotinf.cc);
									if(phidetc<phidet) phidet=phidetc;
								}//right side if
							}//soulution for loop
						}
					}
				}
			}
		}
		phic+=dephi;
	}
	while((phic-dephi)<phidet);
}

void clustertype::bottomrot2(const vtype& rp1,const vtype& rp2,vector<rotattempt>& atvec) // add event where clusters center of mass is below the axis
{
	rotattempt att;
	rotinftype rotinf;
	tobelow2(cm,rp1,rp2,rotinf);
	vtype cmc=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);
	att.contacts=20; // cluster is stable so this value is irrelevant. Value 20 signals that cluster hangs
	att.stable=true;
	att.p1=p1;
	att.d1=d1;
	att.p2=p2;
	att.d2=d2;
	att.p3=-1;
	att.d3=-1;
	att.tp=-1;
	att.td=-1;
	att.f=0;
	att.cm=cmc;
	att.rotinf=rotinf;
	atvec.push_back(att);
	assert(att.cm.z()<=cm.z()+epsz); 
	return;
}

void clustertype::grot2(const vtype rp1,const vtype rp2,vector<rotattempt>& atvec) // add event for contact with the ground while rolling on two contacts
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	for(int i=0;i<r.size();i++)
	{
		if(toground2(rp1,rp2,cm,r[i],R[i],rotinf1,rotinf2))
		{
			vtype cmc1=rot(cm,rp1,rotinf1.cc,rotinf1.ss,rotinf1.n);
			vtype cmc2=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
			if(rightside(cm-rp1,cmc1-rp1,rotinf1.n))
			{
				att.contacts=3;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=d2;
				att.p2=p2;
				att.d3=i;
				att.p3=-2;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}

			if(rightside(cm-rp1,cmc2-rp1,rotinf2.n))
			{
				att.contacts=3;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=d2;
				att.p2=p2;
				att.d3=i;
				att.p3=-2;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}
		}
	}
}

void clustertype::hrot2(const vtype& rp1,const vtype& rp2,vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	if(r[d1].z()>rp1.z()) 
	{
		if(tohorizontal2(rp1,rp2,cm,r[d1],rotinf1,rotinf2))
		{

			vector<rotinftype> rotv;
			rotv.push_back(rotinf1);rotv.push_back(rotinf2);
			for(int i=0;i<rotv.size();i++)
			{
				rotinftype rotinf=rotv[i];
				vtype cmc=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);
				vtype rd1r=rot(r[d1],rp1,rotinf.cc,rotinf.ss,rotinf.n);			
				vtype rd2r=rot(r[d2],rp1,rotinf.cc,rotinf.ss,rotinf.n);
				if(rightside(cm-rp1,cmc-rp1,rotinf.n)&&rd2r.z()<rp2.z())
				{
					att.contacts=0;
					att.stable=false;
					att.p1=-1;
					att.d1=-1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=p1;
					att.td=d1;
					att.f=0;
					att.cm=cmc;
					att.rotinf=rotinf;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}			
				}
			}
		}
	}

	if(r[d2].z()>rp2.z()) 
	{
		if(tohorizontal2(rp2,rp1,cm,r[d2],rotinf1,rotinf2))
		{

			vector<rotinftype> rotv;
			rotv.push_back(rotinf1);rotv.push_back(rotinf2);
			for(int i=0;i<rotv.size();i++)
			{
				rotinftype rotinf=rotv[i];
				vtype cmc=rot(cm,rp2,rotinf.cc,rotinf.ss,rotinf.n);
				vtype rd1r=rot(r[d1],rp2,rotinf.cc,rotinf.ss,rotinf.n);			
				vtype rd2r=rot(r[d2],rp2,rotinf.cc,rotinf.ss,rotinf.n);
				if(rightside(cm-rp2,cmc-rp2,rotinf.n)&&rd1r.z()<rp1.z())
				{
					att.contacts=0;
					att.stable=false;
					att.p1=-1;
					att.d1=-1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=p2;
					att.td=d2;
					att.f=0;
					att.cm=cmc;
					att.rotinf=rotinf;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}			
				}
			}
		}
	}
}

void clustertype::disconnect(const vtype& rp1,const vtype& rp2, vector<rotattempt>& atvec,surroundtype& su) // adds events for the grazing loss of contacts, see Fig 6a of the paper
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	vtype rd1=r[d1];
	vtype rd2=r[d2];

	// treatment of contact 2 
	bool push12=g(rp1,cm,rd2,rp2)<0; // is contact 2 closing while rolling just on contact 1?
	if(push12&&pdisconnect(cm,rp1,rp2,rd1,rd2,rotinf1,rotinf2)) // if so find position of the cluster where d2 moves tangentially to p2, while rotation is on p1 (see paper Eq. 25)
	{
		// starting solution 1
		vtype cmc1=rot(cm,rp1,rotinf1.cc,rotinf1.ss,rotinf1.n);
		if(rightside(cm-rp1,cmc1-rp1,rotinf1.n))
		{
			if(positivecircle(rp1,rp2,rd2,rotinf1)) // checking the condition (31) of the paper
			{
				// condition (31) not satisfied. Standard treatment.
				att.contacts=1;
				att.stable=false;
				att.p1=p1;
				att.d1=d1;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p2;
				att.td=d2;
				att.f=1;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}
			else
			{
				// condition (31) satisfied. Cluster needs to be rotated by dphi before contact can be lost
				rotinfphi.n=rotinf1.n;
				rotinftype rotinfdphi=rotn(rotinf1,rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd2r=rot(rd2,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				if(rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&g(rp1,cmcdph,rd2r,rp2)>epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p1;
					att.d1=d1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=1;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}
				}
			}
		}

		// starting solution 2
		vtype cmc2=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
		if(rightside(cm-rp1,cmc2-rp1,rotinf2.n))
		{
			if(positivecircle(rp1,rp2,rd2,rotinf2))
			{
				att.contacts=1;
				att.stable=false;
				att.p1=p1;
				att.d1=d1;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p2;
				att.td=d2;
				att.f=1;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}
			else
			{
				rotinfphi.n=rotinf2.n;
				rotinftype rotinfdphi=rotn(rotinf2,rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd2r=rot(rd2,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				if(rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&g(rp1,cmcdph,rd2r,rp2)>epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p1;
					att.d1=d1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=1;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}
				}	
			}
		}
	}

	// end of treatment of contact2

	// treatment of contact 1
	bool push21=g(rp2,cm,rd1,rp1)<0;
	if(push21&&pdisconnect(cm,rp2,rp1,rd2,rd1,rotinf1,rotinf2))
	{
		// solution 1
		vtype cmc1=rot(cm,rp1,rotinf1.cc,rotinf1.ss,rotinf1.n);
		if(rightside(cm-rp1,cmc1-rp1,rotinf1.n))
		{
			if(positivecircle(rp2,rp1,rd1,rotinf1))
			{
				att.contacts=1;
				att.stable=false;
				att.p1=p2;
				att.d1=d2;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p1;
				att.td=d1;
				att.f=1;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}
			else
			{
				rotinfphi.n=rotinf1.n;
				rotinftype rotinfdphi=rotn(rotinf1,rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd1r=rot(rd1,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				if(rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&g(rp2,cmcdph,rd1r,rp1)>epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p2;
					att.d1=d2;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=0;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}			
				}
			}
		}

		// solution 2
		vtype cmc2=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
		vtype rd2c2=rot(rd2,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
		if(rightside(cm-rp1,cmc2-rp1,rotinf2.n))
		{
			if(positivecircle(rp2,rp1,rd1,rotinf2))
			{
				att.contacts=1;
				att.stable=false;
				att.p1=p2;
				att.d1=d2;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p1;
				att.td=d1;
				att.f=1;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+epsz) atvec.push_back(att);
			}
			else
			{
				rotinfphi.n=rotinf2.n;
				rotinftype rotinfdphi=rotn(rotinf2,rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd1r=rot(rd1,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);

				vtype cmcn=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
				vtype cmcnd=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);

				if(rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&g(rp2,cmcdph,rd1r,rp1)>epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p2;
					att.d1=d2;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=0;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+epsz) {atvec.push_back(att);}			
				}
			}
		}
	}
}

bool clustertype::positivecircle(const vtype& rp1,const vtype& rp2,const vtype& rd2,const rotinftype& rotinf) // check condition (31) of the paper
{
	vtype cmr=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);
	vtype rdr=rot(rd2,rp1,rotinf.cc,rotinf.ss,rotinf.n);

	rotinftype rotinf180;
	vtype n=v(ez,cmr-rp1);n=(1/norm(n))*n;
	rotinf180.cc=-1;
	rotinf180.ss=0;
	rotinf180.n=n;

	vtype rdr180=rot(rdr,rp1,rotinf180.cc,rotinf180.ss,rotinf180.n);
	return norm(rdr180-rp2)>norm(rdr-rp2);
}


void clustertype::step2(const rotattempt& att,surroundtype& su,const vtype& rp) // update cluster state for rotation on two contacts
{
	contacts=att.contacts;
	stable=att.stable;
	p1=att.p1;
	d1=att.d1;
	p2=att.p2;
	d2=att.d2;
	p3=att.p3;
	d3=att.d3;
	td=att.td;
	tp=att.tp;
	f=att.f;
	cm=rot(cm,rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
	for(int i=0;i<r.size();i++) r[i]=rot(r[i],rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
}


void clustertype::find_steppest_descent(vector<vtype> rp,vector<vtype> rd,vector<int> p,vector<int> d,const vtype& cmc,rotattempt& att) // find steepest descent trajectory of the center of mass while considering all contacts
{
	gradrotattempt gatt;
	vector<gradrotattempt> gradevent;
	assert(rp.size()==rd.size());
	assert(p.size()==d.size());
	assert(rp.size()==p.size());
	int contact_number=rp.size();
	// Add event for cluster falling, if allowed
	bool can_fall=true;
	for(int i=0;i<contact_number;i++) {can_fall=can_fall&&rd[i].z()<rp[i].z();}
	if(can_fall) 
	{
		gatt.contacts=0;
		gatt.stable=false;
		gatt.d1=-1;
		gatt.p1=-1;
		gatt.d2=-1;
		gatt.p2=-1;
		gatt.d3=-1;
		gatt.p3=-1;
		gatt.td=-1;
		gatt.tp=-1;
		gatt.f=0;
		gatt.cm=cmc;
		gatt.gz=-1;
		gradevent.push_back(gatt);
	}

	// Add trajectories for rolling on one contact, if they exist 
	for(int i=0;i<contact_number;i++)
	{
		bool can_rotate=true;
		for(int j=0;j<contact_number;j++) 
		if(i!=j) 
		{
			can_rotate=can_rotate&&g(rp[i],cmc,rd[j],rp[j])>epstang;
		}
		if(can_rotate)
		{
			gatt.contacts=1;
			gatt.stable=false;
			gatt.d1=d[i];
			gatt.p1=p[i];
			gatt.d2=-1;
			gatt.p2=-1;
			gatt.d3=-1;
			gatt.p3=-1;
			gatt.td=-1;
			gatt.tp=-1;
			gatt.f=0;
			gatt.cm=cmc;
			gatt.gz=s(ez,grad(rp[i],cmc));
			gradevent.push_back(gatt);
		}

	}

	// Add trajectories for rolling on two contacts, if they exist
	for(int i=1;i<contact_number;i++)
	for(int j=0;j<i;j++)
	{
		bool can_roll=true;
		for(int k=0;k<contact_number;k++) if(k!=i&&k!=j) {can_roll=can_roll&&g2(rp[i],rp[j],cmc,rd[k],rp[k])>epstang;}
		if(can_roll)
		{
			gatt.contacts=2;
			gatt.stable=false;
			gatt.d1=d[i];
			gatt.p1=p[i];
			gatt.d2=d[j];
			gatt.p2=p[j];
			gatt.d3=-1;
			gatt.p3=-1;
			gatt.td=-1;
			gatt.tp=-1;
			gatt.f=0;
			gatt.cm=cmc;
			gatt.gz=s(ez,grad2(rp[i],rp[j],cmc));
			gradevent.push_back(gatt);
		}
	}

	if(!gradevent.empty()) // if there are trajectories find the steepest
	{
		gatt=gradevent[0];
		for(int ig=0;ig<gradevent.size();ig++)
		{
			if(gradevent[ig].gz<gatt.gz) gatt=gradevent[ig];
		}
		att.contacts=gatt.contacts;
		att.stable=gatt.stable;
		att.d1=gatt.d1;	
		att.p1=gatt.p1;
		att.d2=gatt.d2;
		att.p2=gatt.p2;
		att.d3=gatt.d3;
		att.p3=gatt.p3;
		att.td=gatt.td;
		att.tp=gatt.tp;
		att.f=gatt.f;
		att.cm=gatt.cm;
	}
	else//if there are no trajectories cluster is stable
	{
		att.contacts=rp.size();
		att.stable=true;
		att.d1=d[0];	
		att.p1=p[0];
		att.d2=d[1];
		att.p2=p[1];
		att.d3=d[2];
		att.p3=p[2];
		att.td=-1;
		att.tp=-1;
		att.f=0;
		att.cm=cmc;
	}
}

void clustertype::randrot()
{
	// random basis way

	vtype nxp; // first vector
	do
	{
		nxp.x()=ranfs();
		nxp.y()=ranfs();
		nxp.z()=ranfs();
	}
	while(norm(nxp)>1||norm(nxp)<0.5);
	nxp=(1/norm(nxp))*nxp;

	vtype nyp; // second vector
	do
	{
		nyp.x()=ranfs();
		nyp.y()=ranfs();
		nyp.z()=ranfs();
	}
	while(norm(nyp)>1||norm(nyp)<0.5||fabs(s(nyp,nxp))<0.1);

	nyp=nyp-nxp*s(nyp,nxp); // nyp is now normal to nxp
	nyp=(1/norm(nyp))*nyp;

	vtype nzp=v(nxp,nyp); // nzp normal to both nxp and nyp

	for(int i=0;i<r.size();i++) 
	{
		r[i]=cm+(r[i].x()-cm.x())*nxp+(r[i].y()-cm.y())*nyp+(r[i].z()-cm.z())*nzp;
	}
}
