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
    return a.x * u() + a.y * v() + a.z * w();
  }
  void build_from_w(const Vec3 &normal);
  void branchlessONB(const Vec3 &normal);
  void pbrtONB(const Vec3 &normal);


  Vec3 axis[3];
};

void ONB::build_from_w(const Vec3& normal) {
  axis[2] = Vec3::unit_vector(normal);
  Vec3 a;
  if (fabs(w().x) > 0.9) {
    a = Vec3(0, 1, 0);
  } else {
    a = Vec3(1, 0, 0);
  }
  axis[1] = Vec3::unit_vector(Vec3::cross(w(), a));
  axis[0] = Vec3::cross(w(), v());
}

void ONB::branchlessONB(const Vec3 &n)
{
  float sign = copysignf(1.0f, n.z);
  const float a = -1.0f / (sign + n.z);
  const float b = n.x * n.y * a;
  axis[2] = n;
  axis[0] = Vec3(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
  axis[1] = Vec3(b, sign + n.y * n.y * a, - n.y);
}

void ONB::pbrtONB(const Vec3 &v1)
{
  axis[2] = v1;
  if (std::abs(v1.x) > std::abs(v1.y))
  {
    float invLen = 1.0f / std::sqrtf(v1.x * v1.x + v1.z * v1.z);
    axis[1] = Vec3(-v1.z * invLen, 0, v1.x * invLen);
  }
  else
  {
    float invLen = 1.0f / std::sqrtf(v1.y * v1.y + v1.z * v1.z);
    axis[1] = Vec3(0, v1.z * invLen, -v1.y * invLen);
  }
  axis[0] = Vec3::cross(axis[2], axis[1]);
}