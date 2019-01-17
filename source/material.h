#pragma once
#include "ray.h"
#include "hitable.h"
class Material {
public:
	virtual bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered, float& pdf) const = 0;
  virtual Vec3 Le(float u, float v, const Vec3& p) const = 0;
  virtual float Pdf(float theta) const = 0;
};