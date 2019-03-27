#pragma once

#include "ray.h"

class BSDF;
class HitRecord
{
public:
	float t;
	Vec3 p;
  Vec3 wo;
	Vec3 normal;
	std::unique_ptr<BSDF> bsdf;
};

class Hitable {
public:
	virtual bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const = 0;
  virtual Vec3 generateSampleOnSurface() const = 0;
  virtual float Pdf() const = 0;
};
