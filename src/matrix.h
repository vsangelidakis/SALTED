#ifndef _MATRIX_
#define _MATRIX_
#include "SimpleVector.h"

class matrix
{
	public:
	matrix(const vtype& v1,const vtype& v2,const vtype& v3);

	public:
	long double m11,m12,m13;
	long double m21,m22,m23;
	long double m31,m32,m33;
};

vtype ma(const matrix& m, const vtype& v);

#endif



