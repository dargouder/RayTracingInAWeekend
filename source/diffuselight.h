#pragma once

#include "material.h"


class DiffuseLight : public Material {
public:
  DiffuseLight(Vec3 &a) : emit(a) {}
  virtual bool Scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const { return false; }
  virtual Vec3 emitted(float u, float v, const Vec3& p) const { return emit; }
  Vec3 emit;
};