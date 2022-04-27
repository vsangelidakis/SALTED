#ifndef _vtype_
#define _vtype_
#include <iostream>
#include <cmath>
using namespace std;

class vtype { // Vector class
  // Add operators for input >> and output << stream for vectors
  friend istream & operator >> (istream & is, vtype & v) {
    is >> v._x >> v._y >> v._z;
    return is;
  }
  
  friend ostream & operator << (ostream & os, const vtype & v) {
    os << v._x << " " << v._y << " " << v._z;
    return os;
  }
  
  // Add + - operators for two vectors
  friend vtype operator + (const vtype & v1, const vtype & v2) {
    vtype res(v1);
    res+=v2;
    return res;
  }
  
  friend vtype operator - (const vtype & v1, const vtype & v2) {
    vtype res(v1);
    res-=v2;
    return res;
  }
  
    // Add * operator for number-vector and vector-number multiplications
  friend vtype operator * (long double c, const vtype & p) {
    vtype res=p;
    res*=c;
    return res;
  }
  
  friend vtype operator * (const vtype & p, long double c) {
    return c*p;
  }

  // norm2: calculates the squared magnitude of a vector  
   friend long double norm2(const vtype & v) {
    return v._x*v._x+v._y*v._y+v._z*v._z; //x^2+y^2+z^2
  }

  // norm: calculates the magnitude of a vector  
  friend long double norm(const vtype & v) {
    return sqrt(v._x*v._x+v._y*v._y+v._z*v._z); //sqrt(x^2+y^2+z^2)
  }

  // s: calculates the dot product of two vectors
  friend long double s(const vtype & v1, const vtype & v2) {
    return v1._x*v2._x+v1._y*v2._y+v1._z*v2._z; // dot= x1*x2+y1*y2+z1*z2
  }

  // v: calculates the cross product of two vectors
  friend vtype v(const vtype & v1, const vtype & v2) {
   vtype res;
   res._x=v1._y*v2._z-v1._z*v2._y; // x=y1*z2-z1*y2
   res._y=v1._z*v2._x-v1._x*v2._z; // y=z1*x2-x1*z2
   res._z=v1._x*v2._y-v1._y*v2._x; // z=x1*y2-y1*x2
   return res;    
  }

  // rot: calculates the rotation of a vector given an axis and angle of rotation, using the Rodrigues rotation formula: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
  //    r: vector to be rotated
  //    rp: some sort of centroid, which is initially subtracted and then added back after vector rotation
  //    cc,ss: cos,sin of angle of rotation
  //    n: unit vector
  friend vtype rot(const vtype& r,const vtype& rp,const long double cc,const long double ss,const vtype& n)
  {
	  vtype rc=r-rp; // substract rp //TODO: Check this
	  return (1-cc)*s(n,rc)*n + cc*rc + ss*v(n,rc) + rp; // add rp back, after rotation
  }

  // swapxy: swaps the X-Y coordinates of vector
  friend void swapxy(vtype & v) {
   vtype hv(v);
   v._x=hv._y;   v._y=hv._x;   v._z=hv._z;	
  }
  
  // swapxz: swaps the X-Z coordinates of vector
  friend void swapxz(vtype & v) {
   vtype hv(v);
   v._x=hv._z;   v._y=hv._y;   v._z=hv._x;	
  }
  
  // swapyz: swaps the Y-Z coordinates of vector
  friend void swapyz(vtype & v) {
   vtype hv(v);
   v._x=hv._x;   v._y=hv._z;   v._z=hv._y;	
  }



public:
  explicit vtype(long double x=0,long double y=0,long double z=0): _x(x), _y(y), _z(z){};

  long double & x() {return _x;}
  long double x() const {return _x;}
  long double & y() {return _y;}
  long double y() const {return _y;}
  long double & z() {return _z;}
  long double z() const {return _z;}

  // add per-dimension (x,y,z) operators += -= *= for vectors
  const vtype & operator += (const vtype & p){
    _x+=p._x; _y+=p._y; _z+=p._z;
    return *this;
  }
  const vtype & operator -= (const vtype & p){
    _x-=p._x; _y-=p._y; _z-=p._z;
    return *this;
  }
  const vtype & operator *= (long double c){
    _x*=c; _y*=c; _z*=c;
    return *this;
  }
private:
  long double _x,_y,_z;
};

const vtype null(0,0,0); //FIXME: Remove if unused
const vtype ex(1,0,0);   //FIXME: Remove if unused
const vtype ey(0,1,0);   //FIXME: Remove if unused
const vtype ez(0,0,1);   //FIXME: Remove if unused
#endif

/* unary plus and minus should be added */
