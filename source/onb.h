#pragma once

#include "vec3.h"

class ONB {
 public:
  ONB() {}

  inline Vec3 operator()(int i) { return axis[i]; }

  Vec3 u() const { return axis[0]; }
  Vec3 v() const { return axis[1]; }
  Vec3 w() const { return axis[2]; }

  Vec3 local(float a, float b, float c) const {
    return a * u() + b * v() + c * w();
  }
  Vec3 local(const Vec3& a) const {
    return a.x() * u() + a.y() * v() + a.z() * w();
  }
  void build_from_w(const Vec3&);

  Vec3 axis[3];
};

void ONB::build_from_w(const Vec3& power) {
  axis[2] = Vec3::unit_vector(power);
  Vec3 a;
  if (fabs(w().x()) > 0.9) {
    a = Vec3(0, 1, 0);
  } else {
    a = Vec3(1, 0, 0);
  }
  axis[1] = Vec3::unit_vector(Vec3::cross(w(), a));
  axis[0] = Vec3::unit_vector(Vec3::cross(w(), v()));
}
