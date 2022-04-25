#include "rotinf.h"
#include "common.h"
#include <cassert>

rotinftype rotn(const rotinftype& rotinf1,const rotinftype& rotinf2)
{
	assert(rotinf1.n.x()==rotinf2.n.x()&&rotinf1.n.y()==rotinf2.n.y()&&rotinf1.n.z()==rotinf2.n.z());
	rotinftype res;
	res.cc=rotinf1.cc*rotinf2.cc-rotinf1.ss*rotinf2.ss;	
	res.ss=rotinf1.cc*rotinf2.ss+rotinf1.ss*rotinf2.cc;
	res.n=rotinf1.n;

	return res;
}

