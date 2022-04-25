
#include "numerics.h"

bool pcontact0(vtype r1,long double R1,vtype r2,long double R2,long double& dtl)//finds translation length for collision of one sphere with another fixed sphere, section 3.1 of the paper
{
	long double ds,ds2,dl,zc;

	ds2=pow(r1.x()-r2.x(),2)+pow(r1.y()-r2.y(),2);dl=R1+R2;
	if(ds2<dl*dl) 
	{
		ds=sqrt(ds2);
		zc=r2.z()+sqrt(dl*dl-ds2);dtl=r1.z()-zc;
		return true;
	}
	else return false;
}

bool pcontact(vtype rp,const vtype& n,vtype rc,long double Rc,vtype rs,long double Rs,rotinftype& rotinf1,rotinftype& rotinf2)//solution for situation in Fig 5 of the paper.
{
	vtype rci=rc,rpi=rp,rsi=rs;
	long double R=Rc+Rs;
	rp=rp+s(rc-rp,n)*n;//move reference point to the plane
	vtype rsp=rs+s(rp-rs,n)*n;//projection point of surround sphere to the plain
	if(norm(rp-rsp)<eps) return false;//no contact if they are exactly on the line
	if(norm2(rsp-rs)<R*R)
	{
		rc=rc-rp;//moving to rp as coordinate system beginning
		rsp=rsp-rp;
		rs=rs-rp;

		long double Rsp=sqrt(R*R-norm2(rsp-rs));//diameters
		long double Rcd=norm(rc);
		if(norm2(rsp)<(Rsp+Rcd)*(Rsp+Rcd))
		{
			vtype ext=rsp;ext=1/norm(ext)*ext;//orts of new system
			vtype eyt=v(n,ext);eyt=1/norm(eyt)*eyt;
			vtype ezt=n;
			matrix m(ext,eyt,ezt);
			vtype v1,v2;
			if(circle_intersect(Rcd,Rsp,norm(rsp),v1,v2))
			{
				v1=ma(m,v1);
				v2=ma(m,v2);
				v1=1/norm(v1)*v1;
				v2=1/norm(v2)*v2;
				rc=1/norm(rc)*rc;//rc is not needed anymore except as an ort

				rotinf1.n=n;
				rotinf1.cc=s(rc,v1);
				rotinf1.ss=s(n,v(rc,v1));

				rotinf2.n=n;
				rotinf2.cc=s(rc,v2);
				rotinf2.ss=s(n,v(rc,v2));

				vtype rcr;
				long double dst;

				rcr=rot(rci,rpi,rotinf1.cc,rotinf1.ss,rotinf1.n);
				dst=fabs(norm(rcr-rsi)-Rc-Rs);
				if(dst>error)  {error=dst;}

				rcr=rot(rci,rpi,rotinf2.cc,rotinf2.ss,rotinf2.n);
				dst=fabs(norm(rcr-rsi)-Rc-Rs);
				if(dst>error)  {error=dst;}
				return true;
			}
			else return false;
		}
		else return false;						
	}
	else return false;
}

bool circle_intersect(long double Rcd,long double Rsp,long double d,vtype& v1,vtype& v2)//intersection between two circles one at (0,0) another at (d,0)
{
	if(fabs(Rcd-Rsp)<fabs(d)&&fabs(Rcd+Rsp)>fabs(d))
	{
		long double x=(d*d+(Rcd*Rcd-Rsp*Rsp))/(2*d);
		if(Rcd*Rcd-x*x<0) return false;

		long double y=sqrt(Rcd*Rcd-x*x);
		v1.x()=x;
		v1.y()=y;
		v1.z()=0;

		v2.x()=x;
		v2.y()=-y;
		v2.z()=0;
		return true;
	}
	else return false;
}


bool pcontact1(vtype rp,vtype cmo,vtype rcb,long double Rcb,vtype rsb,long double Rsb,rotinftype& rotinf1,rotinftype& rotinf2)//adapts the function pcontact for collision of one sphere with fixed, whilec cluster is rolling on one contact
{
	vtype ns;
	ns=v(ez,cmo-rp);ns=(1/norm(ns))*ns;

	bool happens=pcontact(rp,ns,rcb,Rcb,rsb,Rsb,rotinf1,rotinf2);
	return happens;
}

void tohorizontal(const vtype& cm,const vtype& rc,const vtype& rs,rotinftype& rotinf)//see section 3.2.3 of the paper
{
	vtype n=v(ez,cm-rs);n=(1/norm(n))*n;
	vtype rd=rc-rs;
	rd=rd-s(rd,n)*n;rd=(1/norm(rd))*rd;
	vtype rdn=v(n,ez);rdn=(1/norm(rdn))*rdn;
	rotinf.cc=s(rd,rdn);
	rotinf.ss=s(n,v(rd,rdn));
	rotinf.n=n;
}


void tobelow1(const vtype& cm,const vtype& rp,rotinftype& rotinf)//see section 3.2.4 of the paper
{
	assert(cm.x()!=rp.x()||cm.y()!=rp.y());
	vtype cmr=cm-rp;
	vtype n=v(ez,cmr);n=(1/norm(n))*n;
	cmr=(1/norm(cmr))*cmr;
	rotinf.cc=s(cmr,-1*ez);
	rotinf.ss=s(n,v(cmr,-1*ez));
	rotinf.n=n;
	assert(rotinf.ss>0);
}


bool toground1(const vtype&rs,vtype cm,vtype rc,const long double Rc,rotinftype& rotinf1,rotinftype& rotinf2)//see section 3.2.5 of the paper
{
	vtype n=v(ez,cm-rs);n=(1/norm(n))*n; 

	cm=cm-rs;
	rc=rc-rs;

	vtype rp=s(rc,n)*n;
	long double R=norm(rc-rp);

	vtype npl=ez;
	vtype rpl(0,0,-rs.z()+Rc);

	vtype rcr=rc-rp;rcr=(1/norm(rcr))*rcr;

	vtype sp,sm;

	if(planecontact(rp,n,R,rpl,npl,sp,sm))
	{
		vtype spr=sp-rp;spr=(1/norm(spr))*spr;
		vtype smr=sm-rp;smr=(1/norm(smr))*smr;

		rotinf1.n=n;
		rotinf1.cc=s(spr,rcr);
		rotinf1.ss=s(v(rcr,spr),n);

		rotinf2.n=n;
		rotinf2.cc=s(smr,rcr);
		rotinf2.ss=s(v(rcr,smr),n);

		return true;
	}
	return false;
}

long double g(const vtype& rp, const vtype& cm, const vtype& rc, const vtype& rs)//determines how "fast" rc is moving towards of away from rs, rotation on one contact
{
	vtype cmr=cm-rp;
	vtype n=v(ez,cmr);n=(1/norm(n))*n;
	return s(v(n,(rc-rp)),rc-rs)/(norm(rc-rp)*norm(rc-rs));
}

vtype grad(const vtype& rp,const vtype& cm)//direction of motion of center of mass when moving on one contact
{
	vtype cmr=cm-rp;
	vtype n=v(ez,cmr);
	vtype gr=v(n,cmr);
	gr=(1/norm(gr))*gr;
	return gr;
}

bool pcontact2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rc,const long double Rcb,const vtype& rs,const long double Rsb,rotinftype& rotinf1,rotinftype& rotinf2)//adapts the function pcontact for collision of one sphere with fixed, whilec cluster is rolling on two contacts, section 3.3.2 of the paper
{
	vtype n=rp2-rp1;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp1))>0) n=-1*n;
	return pcontact(rp1,n,rc,Rcb,rs,Rsb,rotinf1,rotinf2);
}

bool pdisconnect(vtype cm,vtype rp1,vtype rp2,vtype rd1,vtype rd2,rotinftype& rotinf1,rotinftype& rotinf2) //section 3.3.3 of the paper
{
	vtype n=rp2-rp1;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp1))>0) n=-1*n;//in order to have cos>0 and sin>0 for small rotatons

	//all relative to rp1
	rp2-=rp1;
	rd1-=rp1;
	rd2-=rp1;
	cm-=rp1;

	long double A=s(ez,rp2)*s(cm,rd2)/s(cm,rp2);

	vtype npl=ez;
	vtype rpl=A*ez;

	vtype rp=s(rd2,n)*n;
	long double R=norm(rd2-rp);
	vtype sp,sm;

	if(planecontact(rp,n,R,rpl,npl,sp,sm))
	{
		rotinf1.n=n;
		rotinf1.cc=s(sp-rp,rd2-rp)/norm2(rd2-rp);
		rotinf1.ss=s(v(rd2-rp,sp-rp),n)/norm2(rd2-rp);

		rotinf2.n=n;
		rotinf2.cc=s(sm-rp,rd2-rp)/norm2(rd2-rp);
		rotinf2.ss=s(v(rd2-rp,sm-rp),n)/norm2(rd2-rp);
		return true;
	}
	return false;
}

bool tohorizontal2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rd1,rotinftype& rotinf1,rotinftype& rotinf2)//section 3.3.4 of the paper
{
	vtype n=rp2-rp1;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp1))>0) n=-1*n;

	vtype rp=rp1+s(rd1-rp1,n)*n;
	long double R=norm(rd1-rp);
	vtype rpl=rp1;
	vtype npl=ez;
	vtype sp,sm;
	if(planecontact(rp,n,R,rpl,npl,sp,sm))
	{
		rotinf1.n=n;
		rotinf1.cc=s(sp-rp,rd1-rp)/norm2(rd1-rp);
		rotinf1.ss=s(v(rd1-rp,sp-rp),n)/norm2(rd1-rp);

		rotinf2.n=n;
		rotinf2.cc=s(sm-rp,rd1-rp)/norm2(rd1-rp);
		rotinf2.ss=s(v(rd1-rp,sm-rp),n)/norm2(rd1-rp);
		long double weps=1e-08;
		if(fabs(rp1.z()-sp.z())>weps) ds("Warning: tohorizontal2: fabs(rp1.z()-sp.z())>weps");
		if(fabs(rp1.z()-sm.z())>weps) ds("Warning: tohorizontal2: fabs(rp1.z()-sm.z())>weps");
		if(fabs(norm(sp-rp1)-norm(rd1-rp1))>weps) ds("Warning: tohorizontal2: fabs(norm(sp-rp1)-norm(rd1-rp1))>weps");
		if(fabs(norm(sm-rp1)-norm(rd1-rp1))>weps) ds("Warning: tohorizontal2: fabs(norm(sm-rp1)-norm(rd1-rp1))>weps");
		if(fabs(s(sp-rp1,n)-s(rd1-rp1,n))>weps) ds("Warning: tohorizontal2: fabs(s(sp-rp1,n)-s(rd1-rp1,n))>weps");
		if(fabs(s(sm-rp1,n)-s(rd1-rp1,n))>weps) ds("Warning: tohorizontal2: fabs(s(sm-rp1,n)-s(rd1-rp1,n))>weps");
		return true;
	}
	else return false;
}

void tobelow2(const vtype& cm,const vtype& rp1,const vtype& rp2,rotinftype& rotinf)///section 3.3.5 of the paper
{
	vtype cmr=cm-rp1;
	vtype n=rp2-rp1;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp1))>0) n=-1*n;

	vtype rp=rp1+s(cmr,n)*n;

	vtype ncmr=cm-rp;ncmr=(1/norm(ncmr))*ncmr;
	vtype cmrr=v(n,v(n,ez));
	vtype ncmrr=cmrr;ncmrr=(1/norm(ncmrr))*ncmrr;
	vtype nr=v(ncmr,ncmrr);nr=(1/norm(nr))*nr;

	rotinf.cc=s(ncmr,ncmrr);
	rotinf.ss=s(n,v(ncmr,ncmrr));
	rotinf.n=n;
	assert(rotinf.ss>0);

	long double weps=1e-08;

	vtype nvp=v(ncmr,ncmrr);nvp=(1/norm(nvp))*nvp;
	if(norm(nvp-n)>weps) ds("Warning: tobelow2: correct normal and recalculated normal do not match, nvp-n=",nvp-n);

	if(fabs(rotinf.cc*rotinf.cc+rotinf.ss*rotinf.ss-1)>weps) ds("Warning: tobelow2: sum cc^2+ss^2 differs from 1, cc^2+ss^2-1",rotinf.cc*rotinf.cc+rotinf.ss*rotinf.ss-1);

	vtype nim=v(ncmrr,n);nim=(1/norm(nim))*nim;
	if(fabs(s(nim,ez))>weps) ds("Warning: tobelow2: rotated position not in vertical plane, s(nim,ez)=",s(nim,ez));
}

bool toground2(const vtype& rp1,vtype rp2,const vtype& cm, vtype rc,const long double Rc,rotinftype& rotinf1,rotinftype& rotinf2)//section 3.3.6 of the paper
{
	vtype n=rp2-rp1;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp1))>0) n=-1*n;

	rc=rc-rp1;

	vtype rp=s(rc,n)*n;
	long double R=norm(rc-rp);

	vtype npl=ez;
	vtype rpl(0,0,-rp1.z()+Rc);

	vtype rcr=rc-rp;rcr=(1/norm(rcr))*rcr;

	vtype sp,sm;

	if(planecontact(rp,n,R,rpl,npl,sp,sm))
	{
		vtype spr=sp-rp;spr=(1/norm(spr))*spr;
		vtype smr=sm-rp;smr=(1/norm(smr))*smr;

		rotinf1.n=n;
		rotinf1.cc=s(spr,rcr);
		rotinf1.ss=s(v(rcr,spr),n);

		rotinf2.n=n;
		rotinf2.cc=s(smr,rcr);
		rotinf2.ss=s(v(rcr,smr),n);

		return true;
	}
	return false;
}


long double g2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rc,const vtype& rs)//analog of g for two contacts
{
	vtype n=rp2-rp1;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp1))>0) n=-1*n;
	vtype g=v(n,rc-rp1);g=(1/norm(g))*g;
	return s(g,rc-rs)/norm(rc-rs);
}

vtype grad2(const vtype& rp1,const vtype& rp2,const vtype& cm)//analog of grad for two contacts
{
	vtype n=rp2-rp1;
	vtype cmr=cm-rp1;
	if(s(ez,v(n,cmr))>0) n=-1*n;
	vtype g=v(n,cmr);g=(1/norm(g))*g;
	return g;
}

bool rightside(const vtype& cmr,const vtype& cmcr,const vtype n)//determines whether center of mass is on the same side of a vertical plane containg rotation axis, before and after the rotation
{
	vtype n1=v(ez,n);
	long double s1=s(n1,cmr);
	long double s2=s(n1,cmcr);
	return s1*s2>=0;
}

bool planecontact(vtype rp,vtype n,long double R,const vtype& rpl,vtype npl,vtype& sp,vtype& sm)//finds contact of a point with a horizontal plane
{
	vtype tr1=rpl;

	rp=rp-tr1;

	if(fabs(n.y())>fabs(n.x()))
	{
		long double A=1+(n.x()*n.x())/(n.y()*n.y());
		long double B=-2*rp.z()*n.z()*n.x()/(n.y()*n.y());
		long double C=(rp.z()*n.z()/n.y())*(rp.z()*n.z()/n.y())+rp.z()*rp.z()-R*R;
		long double dis=B*B-4*A*C;
		if(dis<0) return false;
		long double xp=(-B+sqrt(dis))/(2*A);
		long double yp=(rp.z()*n.z()-xp*n.x())/n.y();

		long double xm=(-B-sqrt(dis))/(2*A);
		long double ym=(rp.z()*n.z()-xm*n.x())/n.y();

		sp.x()=xp+rp.x();sp.y()=yp+rp.y();sp.z()=0;
		sm.x()=xm+rp.x();sm.y()=ym+rp.y();sm.z()=0;
	}
	else
	{
		long double A=1+(n.y()*n.y())/(n.x()*n.x());
		long double B=-2*rp.z()*n.z()*n.y()/(n.x()*n.x());
		long double C=(rp.z()*n.z()/n.x())*(rp.z()*n.z()/n.x())+rp.z()*rp.z()-R*R;
		long double dis=B*B-4*A*C;
		if(dis<0) return false;
		long double yp=(-B+sqrt(dis))/(2*A);
		long double xp=(rp.z()*n.z()-yp*n.y())/n.x();

		long double ym=(-B-sqrt(dis))/(2*A);
		long double xm=(rp.z()*n.z()-ym*n.y())/n.x();

		sp.x()=xp+rp.x();sp.y()=yp+rp.y();sp.z()=0;
		sm.x()=xm+rp.x();sm.y()=ym+rp.y();sm.z()=0;
	}

	sp=sp+tr1;
	sm=sm+tr1;

	return true;
}



/********** Grid operations ***********/

void findindex(const vtype& r,long double R,int& ixl,int& ixr,int& iyl,int& iyr,int& izu)//find adjecent grid indexes upper and lower bound
{
	ixl=floor((r.x()-R-Rmax-xo)/dx);
	ixr=floor((r.x()+R+Rmax-xo)/dx);
	iyl=floor((r.y()-R-Rmax-yo)/dy);
	iyr=floor((r.y()+R+Rmax-yo)/dy);
	izu=floor((r.z()-zo)/dy);if(izu>nz-1) izu=nz-1;
}


/********** Speedup functions ***********/

long double distance_from_axis(const vtype& rc,const vtype& rp,const vtype& cm)//finds distance of point rc from the axis when motion is accoring to rule from Fig 2b of the paper.
{
	vtype rcr=rc-rp;
	vtype cmr=cm-rp;
	vtype n=v(ez,cmr);n=(1/norm(n))*n;
	vtype normal=rcr-s(rcr,n)*n;
	return norm(normal);
}


void findpoints(const vtype& rp,const vtype& cm,const vtype& rc,long double d,vector<vtype>& points)//finds testing positons for search collision, movement on one contact
{
	points.erase(points.begin(),points.end());
	vtype cmr=cm-rp;
	vtype rcr=rc-rp;
	vtype n=v(ez,cmr);n=(1/norm(n))*n;
	long double normaldist=distance_from_axis(rc,rp,cm);
	long double ss2=d/normaldist;
	long double cc2=sqrt(1-ss2*ss2);
	long double cc=1-2*ss2*ss2;
	long double ss=2*ss2*cc2;
	rotinftype rotinf;
	rotinf.n=n;
	rotinf.cc=cc;
	rotinf.ss=ss;
	long double dalpha=2*asin(d/(2*normaldist));
	int nrot=ceil(3.141592653/dalpha);
	vtype rcm=rc;
	for(int i=0;i<nrot;i++)
	{
		points.push_back(rcm);
		rcm=rot(rcm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	}
}

vtype rotpoint(const vtype& r,const vtype& rp,const vtype& cm,long double phi)//rotates point around axis defined by rp and cm (Fig 2b of the papaer)
{
	vtype n=v(ez,cm-rp);n=(1/norm(n))*n;
	return rot(r,rp,cos(phi),sin(phi),n);
}


long double distance_from_axis2(const vtype& rc,const vtype& rp,const vtype& rp2)//finds distance of point rc from the axis when motion is accoring to rule from Fig 2c of the paper.
{
	vtype rcr=rc-rp;
	vtype n=rp2-rp;n=(1/norm(n))*n;
	vtype normal=rcr-s(rcr,n)*n;
	return norm(normal);
}


void findpoints2(const vtype& rp,const vtype& rp2,const vtype& cm,const vtype& rc,long double d,vector<vtype>& points)//finds testing positons for search collision, movement on two contact
{
	points.erase(points.begin(),points.end());
	vtype cmr=cm-rp;
	vtype rcr=rc-rp;
	vtype n=rp2-rp;n=(1/norm(n))*n;
	if(s(ez,v(n,cmr))>0) n=-1*n;
	long double normaldist=distance_from_axis2(rc,rp,rp2);
	long double ss2=d/normaldist;
	long double cc2=sqrt(1-ss2*ss2);
	long double cc=1-2*ss2*ss2;
	long double ss=2*ss2*cc2;
	rotinftype rotinf;
	rotinf.n=n;
	rotinf.cc=cc;
	rotinf.ss=ss;
	long double dalpha=2*asin(d/(2*normaldist));
	int nrot=ceil(3.141592653/dalpha);
	vtype rcm=rc;
	for(int i=0;i<nrot;i++)
	{
		points.push_back(rcm);
		rcm=rot(rcm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	}
}

vtype rotpoint2(const vtype& r,const vtype& rp,const vtype& rp2,const vtype& cm,long double phi)//rotates point around axis defined by rp and cm (Fig 2c of the papaer)
{
	vtype n=rp2-rp;n=(1/norm(n))*n;
	if(s(ez,v(n,cm-rp))>0) n=-1*n;
	return rot(r,rp,cos(phi),sin(phi),n);
}

