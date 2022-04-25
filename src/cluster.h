#ifndef _CLUSTER_
#define _CLUSTER_

#include "common.h"
#include "rotinf.h"
#include "surround.h"
#include "numerics.h"

class surroundtype;

class clusterbase
{
	public:
		int contacts;
		bool stable;
		int d1,p1,d2,p2,d3,p3,td,tp;
		bool f;
		vtype cm;
		bool operator<(const clusterbase& c) const {return (*this).cm.z()<c.cm.z();}; 
		bool operator>(const clusterbase& c) const {return (*this).cm.z()>c.cm.z();}; 
};


class fallattempt : public clusterbase
{

};


class rotattempt : public clusterbase
{
	public:
		rotinftype rotinf;
};


class gradrotattempt : public rotattempt
{
	public:
		long double gz;
};


class clustertype : public clusterbase//full cluster implementation
{
	public:
		vector<vtype> r;
		vector<long double> R;

	public:
		vtype cmass();
		void erase();
		void init_sphere();
		void init_axes_lines_mshuffle(int nsph_per_side);
		void init_cube_mshuffle(int nsph_per_side);
		bool settle(surroundtype& su);
		void periodic(long double lx,long double ly);
		void movecluster(const vtype& v);
		void randrot();

	private:
		void fall(surroundtype& su);
		void bfall(surroundtype& su,vector<fallattempt>& attempts);
		void gfall(vector<fallattempt>& attempts);
		void step0(fallattempt att,surroundtype& su);
		void translate(long double dlt,int d1,int p1,surroundtype& su);

		void rotate1(surroundtype& su);
		void steeperdphi(const vtype& rp,vector<rotattempt>& atvec);
		void brot1(surroundtype& s,vector<rotattempt>& attempts);
		void hrot(const vtype& rp,vector<rotattempt>& attempts);
		void grot1(const vtype& rp,vector<rotattempt>& attempts);
		void bottomrot1(const vtype& rp,vector<rotattempt>&  attempts);
		void step1(const rotattempt& att,const vtype& rp);

		void rotate2(surroundtype& su);
		void steeperdphi2(const vtype& rp1,const vtype& rp2,vector<rotattempt>& atvec);
		void brot2(surroundtype& su,vector<rotattempt>&  attempts);
		void bottomrot2(const vtype& rp1,const vtype& rp2,vector<rotattempt>& atvec);
		void grot2(const vtype rp1,const vtype rp2,vector<rotattempt>&  attempts);
		void hrot2(const vtype& rp1,const vtype& rp2,vector<rotattempt>& atvec);
		void disconnect(const vtype& rp1,const vtype& rp2, vector<rotattempt>&  attempts,surroundtype& su);
		bool positivecircle(const vtype& rp1,const vtype& rp2,const vtype& rd2,const rotinftype& rotinf);
		void step2(const rotattempt& att,surroundtype& su,const vtype& rp);
		void find_steppest_descent(vector<vtype> rp,vector<vtype> rd,vector<int> p,vector<int> d,const vtype& cmc,rotattempt& att);

};


#endif
