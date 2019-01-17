#pragma once

#include "material.h"


class DiffuseLight : public Material {
public:
  DiffuseLight(Vec3 a) : emit(a) {}
  virtual bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered, float& pdf) const { return false; }
  virtual Vec3 Le(float u, float v, const Vec3& p) const { return emit; }
  Vec3 emit;

  float Pdf(float theta) const
  {
    return 1.0f;
  }
};