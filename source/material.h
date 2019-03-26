#pragma once
#include "hitable.h"
#include "ray.h"
class Material {
 public:
  virtual Vec3 f(const Vec3& wo, const Vec3& wi) const = 0;
  virtual bool sample_f(const Ray& r_in, const HitRecord& rec,
                        Vec3& attenuation, Ray& scattered,
                        float& pdf) const = 0;
  virtual Vec3 Le(float u, float v, const Vec3& p) const = 0;
  virtual float ScatteredPdf(const Ray& r_in, const HitRecord& rec,
                             const Ray& scattered) const = 0;

  virtual float pdf(const Vec3& wo, const Vec3& wi) const = 0;
};
