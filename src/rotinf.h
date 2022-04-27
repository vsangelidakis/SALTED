#ifndef _ROTINFTYPE_
#define _ROTINFTYPE_
#include "SimpleVector.h"

class rotinftype
{
	public:
	long double cc,ss; // cos, sin
	vtype n; // normal vector
};

rotinftype rotn(const rotinftype& rotinf1,const rotinftype& rotinf2);

#endif
