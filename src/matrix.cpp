#include "matrix.h"

matrix::matrix(const vtype& v1,const vtype& v2,const vtype& v3)
{ // [3x3] matrix class: v1,v2,v3 are [1x3] vertices which form each row
	m11=v1.x();m12=v2.x();m13=v3.x();
	m21=v1.y();m22=v2.y();m23=v3.y();
	m31=v1.z();m32=v2.z();m33=v3.z();
}

vtype ma(const matrix& m, const vtype& v)
{ // multiplies matrix with vector: [3x3] * [3x1] = [3x1] 
	vtype res;
	res.x()=m.m11*v.x()+m.m12*v.y()+m.m13*v.z();
	res.y()=m.m21*v.x()+m.m22*v.y()+m.m23*v.z();
	res.z()=m.m31*v.x()+m.m32*v.y()+m.m33*v.z();
	return res;
}

// Overload a function for the mirror calculation [1x3] * [3x3]
vtype ma(const vtype& v, const matrix& m)
{ // multiplies vector with matrix: [1x3] * [3x3] = [1x3]
	vtype res;
	res.x()=m.m11*v.x()+m.m12*v.y()+m.m13*v.z();
	res.y()=m.m21*v.x()+m.m22*v.y()+m.m23*v.z();
	res.z()=m.m31*v.x()+m.m32*v.y()+m.m33*v.z();
	return res;
}

//TODO: Create also matrix multiplications [3x3] * [3x3] if needed later on
//TODO: Create functions to extract rows/columns
//TODO: Create functions to diagonalise










