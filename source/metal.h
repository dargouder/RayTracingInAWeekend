#pragma once
#include "material.h"

class Metal : public Material {
public:
	Vec3 albedo;
	float fuzz;
	Metal(const Vec3& a, float f) : albedo(a) {
		if(f < 1){
			fuzz = f;
		} else {
			fuzz = 1;
		}
	}

	virtual bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		Vec3 reflected = reflect(Vec3::unit_vector(r_in.direction()), rec.normal);
		scattered = Ray(rec.p, reflected + fuzz * RandomInUnitSphere());
		attenuation = albedo / M_PI;

		return (Vec3::dot(scattered.direction(), rec.normal) > 0);
	}

  virtual Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0, 0, 0); }

  float Pdf(float theta) const
  {
    return cos(theta) * fuzz / M_PI;
  }
};