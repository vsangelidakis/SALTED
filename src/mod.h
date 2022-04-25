#ifndef _MOD_
#define _MOD_

#include "SimpleVector.h"
#include <cmath>

int imod(int i,int n);
long double rmod(long double x,long double l);
long double moddist(long double x, long double l);
vtype vmod(const vtype& r,long double lx,long double ly);
vtype vmod(const vtype& ri,const vtype& rb,long double lx,long double ly);

#endif
