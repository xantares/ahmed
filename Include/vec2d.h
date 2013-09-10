/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef VEC2D_H
#define VEC2D_H

#include <cmath>

class vec2d
{

  double x[2];

public:
  vec2d() {};
  vec2d(const double d1, const double d2) {
    x[0] = d1;
    x[1] = d2;
  };

  ~vec2d() {};

  double operator[](const unsigned i) const {
    return x[i];
  }

  double sqrsum() const {
    return x[0]*x[0] + x[1]*x[1];
  };
  double norm2() const {
    return sqrt(sqrsum());
  };

  vec2d rotate(const double a) const {
    return vec2d(cos(a)*x[0]-sin(a)*x[1], sin(a)*x[0]+cos(a)*x[1]);
  }

  vec2d operator+ (const vec2d& v) const {
    return vec2d(x[0]+v.x[0], x[1]+v.x[1]);
  }

  vec2d& operator+= (vec2d& v) {
    x[0] += v.x[0];
    x[1] += v.x[1];
    return *this;
  }

  vec2d operator- (const vec2d& v) const {
    return vec2d(x[0]-v.x[0], x[1]-v.x[1]);
  }

  vec2d& operator-= (const vec2d& v) {
    x[0] -= v.x[0];
    x[1] -= v.x[1];
    return *this;
  }

  vec2d operator* (const double a) const {
    return vec2d(a*x[0], a*x[1]);
  }

  double operator* (const vec2d& v) const {
    return x[0]*v.x[0]+x[1]*v.x[1];
  }

  vec2d& operator*= (const double a) {
    x[0] *= a;
    x[1] *= a;
    return *this;
  }
};

#endif
