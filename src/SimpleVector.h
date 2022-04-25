#ifndef _vtype_
#define _vtype_
#include <iostream>
#include <cmath>
using namespace std;

class vtype {
  friend istream & operator >> (istream & is, vtype & v) {
    is >> v._x >> v._y >> v._z;
    return is;
  }
  friend ostream & operator << (ostream & os, const vtype & v) {
    os << v._x << " " << v._y << " " << v._z;
    return os;
  }
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
  friend vtype operator * (long double c, const vtype & p) {
    vtype res=p;
    res*=c;
    return res;
  }
  friend vtype operator * (const vtype & p, long double c) {
    return c*p;
  }
   friend long double norm2(const vtype & v) {
    return v._x*v._x+v._y*v._y+v._z*v._z;
  }
   friend long double norm(const vtype & v) {
    return sqrt(v._x*v._x+v._y*v._y+v._z*v._z);
  }

  friend long double s(const vtype & v1, const vtype & v2) {
    return v1._x*v2._x + v1._y*v2._y+v1._z*v2._z;
  }

  friend vtype v(const vtype & v1, const vtype & v2) {
   vtype res;
   res._x=v1._y*v2._z-v1._z*v2._y;
   res._y=v1._z*v2._x-v1._x*v2._z;
   res._z=v1._x*v2._y-v1._y*v2._x;
   return res;    
  }

  friend vtype rot(const vtype& r,const vtype& rp,const long double cc,const long double ss,const vtype& n)
  {
	  vtype rc=r-rp;
	  return (1-cc)*s(n,rc)*n+cc*rc+ss*v(n,rc)+rp;
  }

   friend void swapxy(vtype & v) {
   vtype hv(v);
   v._x=hv._y;   v._y=hv._x;   v._z=hv._z;	
  }
   friend void swapxz(vtype & v) {
   vtype hv(v);
   v._x=hv._z;   v._y=hv._y;   v._z=hv._x;	
  }
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

const vtype null(0,0,0);
const vtype ex(1,0,0);
const vtype ey(0,1,0);
const vtype ez(0,0,1);
#endif

/* unary plus and minus should be added */
