#ifndef _COMMON_INC_
#define _COMMON_INC_

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <deque>
#include <ctime>
#include "SimpleVector.h"
#include "mod.h"
#include "rotinf.h"
#include "short_cout.h"
#include <gsl/gsl_rng.h>
#include <fenv.h>

using namespace std;

extern int N;
extern int br;

extern long double error;
extern long double error1;
extern long double eps;

extern int nx,ny,nz;
extern long double xo,yo,zo,dx,dy,dz;

extern  long double lx;
extern  long double ly;
extern  long double pxb;
extern  long double pxe;
extern  long double pyb;
extern  long double pye;

extern string cmsfile;
extern string positionsfile;
extern string detailsfile;

extern  long double phi;
extern rotinftype rotinfphi;

extern long double Rmax;
extern long double h0;
extern long double lcd;
extern long double epsz;
extern long double epstang;
extern int particle;

extern long double l_pi;

extern const gsl_rng_type * T;
extern gsl_rng * r;

inline long double ranf() {return gsl_rng_uniform (r);}//uniform random numbers in interval [0,1)
inline long double ranfs() {return -1+2*ranf();}//uniform random numbers in interval [-1,1)
inline long double sgn(long double x) {return (x>=0?1:-1);}//sign function

vtype rand_orientation();


#endif
