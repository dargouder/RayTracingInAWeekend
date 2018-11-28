#pragma once
#include "ray.h"
#include "hitable.h"
class Material {
public:
	virtual bool Scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const = 0;
  virtual Vec3 emitted(float u, float v, const Vec3& p) const = 0;
};