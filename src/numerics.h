#ifndef _NUMERICS_
#define _NUMERICS_
#include "common.h"
#include "rotinf.h"
#include "matrix.h"

/***** Contact numerics *****/
bool pcontact0(vtype r1,long double R1,vtype r2,long double R2,long double& dtl);

bool pcontact(vtype rp,const vtype& n,vtype rc,long double Rc,vtype rs,long double Rs,rotinftype& rotinf1,rotinftype& rotinf2);
bool circle_intersect(long double Rcd,long double Rsp,long double d,vtype& v1,vtype& v2);

bool pcontact1(vtype rp,vtype cmo,vtype rcb,long double Rcb,vtype rsb,long double Rsb,rotinftype& rotinf1,rotinftype& rotinf2);

void tohorizontal(const vtype& cm,const vtype& rc,const vtype& rs,rotinftype& rotinf);
void tobelow1(const vtype& cm,const vtype& rp,rotinftype& rotinf);
bool toground1(const vtype&rs,vtype cm,vtype rc,const long double Rc,rotinftype& rotinf1,rotinftype& rotinf2);

long double g(const vtype& rp, const vtype& cm, const vtype& rc, const vtype& rs);
vtype grad(const vtype& rp,const vtype& cm);

bool pcontact2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rc,const long double Rcb,const vtype& rs,const long double Rsb,rotinftype& rotinf1,rotinftype& rotinf2);
bool pdisconnect(vtype cm,vtype rp1,vtype rp2,vtype rd1,vtype rd2,rotinftype& rotinf1,rotinftype& rotinf2);
bool tohorizontal2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rd1,rotinftype& rotinf1,rotinftype& rotinf2);
void tobelow2(const vtype& cm,const vtype& rp1,const vtype& rp2,rotinftype& rotinf);
bool toground2(const vtype& rp1,vtype rp2,const vtype& cm, vtype rc,const long double Rc,rotinftype& rotinf1,rotinftype& rotinf2);

long double g2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rc,const vtype& rs);
vtype grad2(const vtype& rp1,const vtype& rp2,const vtype& cm);


bool rightside(const vtype& cmr,const vtype& cmcr,const vtype n);
bool planecontact(vtype rp,vtype n,long double R,const vtype& rpl,vtype npl,vtype& sp,vtype& sm);


/***** Grid function *****/
void findindex(const vtype& r,long double R,int& ixl,int& ixr,int& iyl,int& iyr,int& izu);

/********** Speedup functions ***********/

long double distance_from_axis(const vtype& rc,const vtype& rp,const vtype& cm);
void findpoints(const vtype& rp,const vtype& cm,const vtype& rc,long double d,vector<vtype>& points);
vtype rotpoint(const vtype& r,const vtype& rp,const vtype& cm,long double phi);

long double distance_from_axis2(const vtype& rc,const vtype& rp,const vtype& rp2);
void findpoints2(const vtype& rp,const vtype& rp2,const vtype& cm,const vtype& rc,long double d,vector<vtype>& points);
vtype rotpoint2(const vtype& r,const vtype& rp,const vtype& rp2,const vtype& cm,long double phi);

#endif

