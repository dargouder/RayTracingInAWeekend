#pragma once

#include "vec3.h"

class BxDF
{
public:
  // BxDF interface
  // BxDF public data
  virtual Vec3 f(const Vec3 &wo, const Vec3 &wi) const = 0;
  virtual Vec3 Sample_f(const Vec3 &wo, Vec3 &wi, float &pdf) const;
};

class LambertBxdf
{



  float pdf(const Vec3& wi, const Vec3& wo, const Vec3 &n) const
  {
    return 0;
  }
};

