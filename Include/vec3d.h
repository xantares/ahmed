/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef VEC3D
#define VEC3D

#include <iostream>
#include <cmath>

class vec3d;

std::ostream& operator << (std::ostream&, const vec3d&);

class vec3d
{
  friend std::ostream& operator << (std::ostream&, const vec3d&);

public:

  /* Konstruktoren */
  vec3d() {};
  vec3d(const double d) { x[0] = x[1] = x[2] = d; };
  vec3d(const double*);
  vec3d(const double, const double, const double);
  vec3d(const vec3d&);

  /* Elementfunktionen */
  const double* getdata() const { return &x[0]; }
  vec3d& operator = (const vec3d&);                        // vec.  = vec.
  vec3d& operator+= (const vec3d&);                        // vec. += vec.
  vec3d& operator-= (const vec3d&);                        // vec. -= vec.
  vec3d& operator*= (const double);                        // vec. *= double
  vec3d& operator/= (const double);                        // vec. /= double
  vec3d  operator + (const vec3d&) const;                  // vec. + vec.
  vec3d  operator - (const vec3d&) const;                  // vec. - vec.
  vec3d  operator * (const double) const;                  // vec. * double
  vec3d  operator / (const double) const;                  // vec. / double
  double operator * (const vec3d&) const;                  // Skalarprdkt.
  vec3d  operator % (const vec3d&) const;                  // Kreuzprdkt.
  double spat(const vec3d&, const vec3d&) const;           // Spatprodukt

  vec3d  operator - () const {
    return vec3d(-x[0], -x[1], -x[2]);
  }

  double kreuz  (const vec3d&, const unsigned) const;      // Komp. Kreuzp.
  double sqrsum () const;                                  // eukl. Norm ^ 2
  double norm2  () const {
    return sqrt(sqrsum());   // eukl. Norm
  }
  double norm1  () const;                                  // 1-Norm
  vec3d& normalize ();                                     // normiert
  double& operator [] (const unsigned i) {
    return x[i];   // Komponente
  }
  const double& operator [] (const unsigned i) const {
    return x[i];
  }

protected:
  double x[3];
};


inline vec3d::vec3d(const double* p)
{
  x[0] = *p++;
  x[1] = *p++;
  x[2] = *p;
}

inline vec3d::vec3d(const double d1, const double d2, const double d3)
{
  x[0] = d1;
  x[1] = d2;
  x[2] = d3;
}

inline vec3d::vec3d(const vec3d& v)
{
  x[0] = v.x[0];
  x[1] = v.x[1];
  x[2] = v.x[2];
}


// ==================================================================
//                         Elementfunktionen
// ==================================================================

inline vec3d& vec3d::operator = (const vec3d& v)
{
  x[0] = v.x[0];
  x[1] = v.x[1];
  x[2] = v.x[2];
  return *this;
}

inline vec3d& vec3d::operator+= (const vec3d& v)
{
  x[0] += v.x[0];
  x[1] += v.x[1];
  x[2] += v.x[2];
  return *this;
}

inline vec3d& vec3d::operator-= (const vec3d& v)
{
  x[0] -= v.x[0];
  x[1] -= v.x[1];
  x[2] -= v.x[2];
  return *this;
}

inline vec3d& vec3d::operator*= (const double w)
{
  x[0] *= w;
  x[1] *= w;
  x[2] *= w;
  return *this;
}

inline vec3d& vec3d::operator/= (const double w)
{
  const double s = 1.0/w;
  x[0]*=s;
  x[1]*=s;
  x[2]*=s;
  return *this;
}

inline vec3d vec3d::operator + (const vec3d& v) const
{
  return vec3d(x[0]+v.x[0], x[1]+v.x[1], x[2]+v.x[2]);
}

inline vec3d vec3d::operator - (const vec3d& v) const
{
  return vec3d(x[0]-v.x[0], x[1]-v.x[1], x[2]-v.x[2]);
}

inline vec3d vec3d::operator * (const double w) const
{
  return vec3d(w*x[0], w*x[1], w*x[2]);
}

inline vec3d vec3d::operator / (const double w) const
{
  const double s = 1.0/w;
  return vec3d(s*x[0], s*x[1], s*x[2]);
}

inline double vec3d::operator * (const vec3d& v) const
{
  return x[0]*v.x[0] + x[1]*v.x[1] + x[2]*v.x[2];
}

inline vec3d vec3d::operator % (const vec3d& v) const
{
  const double d1 = x[1]*v.x[2] - x[2]*v.x[1];
  const double d2 = x[2]*v.x[0] - x[0]*v.x[2];
  const double d3 = x[0]*v.x[1] - x[1]*v.x[0];

  return vec3d(d1, d2, d3);
}

inline double vec3d::kreuz(const vec3d& v, const unsigned k) const
{
  double r;

  switch (k) {
  case 0 : {
    r = x[1]*v.x[2] - x[2]*v.x[1];
    break;
  }
  case 1 : {
    r = x[2]*v.x[0] - x[0]*v.x[2];
    break;
  }
  case 2 : {
    r = x[0]*v.x[1] - x[1]*v.x[0];
    break;
  }
  default : {
    r = 0.0;
  }
  }

  return r;
}

inline double vec3d::spat(const vec3d& v, const vec3d& w) const
{
  double r = x[0]*v.x[1]*w.x[2] + x[2]*v.x[0]*w.x[1] + x[1]*v.x[2]*w.x[0];
  r -= x[2]*v.x[1]*w.x[0] + x[0]*v.x[2]*w.x[1] + x[1]*v.x[0]*w.x[2];

  return r;
}

inline double vec3d::sqrsum() const
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

inline double vec3d::norm1() const
{
  return fabs(x[0]) + fabs(x[1]) + fabs(x[2]);
}

inline vec3d& vec3d::normalize()
{
  const double s = 1.0/norm2();

  x[0]*=s;
  x[1]*=s;
  x[2]*=s;
  return *this;
}


// ostream - Operator
inline std::ostream& operator << (std::ostream& os, const vec3d& v)
{
  os << '(' << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ')';
  return os;
}

#endif    // VEC3D
