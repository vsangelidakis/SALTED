#include "rotinf.h"
#include "common.h"
#include <cassert>

rotinftype rotn(const rotinftype& rotinf1,const rotinftype& rotinf2)
{
	assert(rotinf1.n.x()==rotinf2.n.x()&&rotinf1.n.y()==rotinf2.n.y()&&rotinf1.n.z()==rotinf2.n.z());
	rotinftype res;
	res.cc=rotinf1.cc*rotinf2.cc-rotinf1.ss*rotinf2.ss; // cos=cos1*cos2-sin1*sin2
	res.ss=rotinf1.cc*rotinf2.ss+rotinf1.ss*rotinf2.cc; // sin=cos1*sin2+sin1*cos2
	res.n=rotinf1.n; //n: normal vector / axis of rotation

	return res;
}

