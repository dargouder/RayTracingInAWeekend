#pragma once
#include "material.h"

class Metal : public Material {
public:
	Vec3 m_kd;
	float fuzz;
	Metal(const Vec3& a, float f) : m_kd(a) {
		if(f < 1){
			fuzz = f;
		} else {
			fuzz = 1;
		}
	}

	virtual bool sample_f(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered, float& pdf) const {
		Vec3 reflected = reflect(Vec3::unit_vector(r_in.direction()), rec.normal);
		scattered = Ray(rec.p, reflected + fuzz * RandomInUnitSphere());
		attenuation = m_kd;

		return (Vec3::dot(scattered.direction(), rec.normal) > 0);
	}

  virtual Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0, 0, 0); }

  float ScatteredPdf(const Ray &r_in, const HitRecord &rec, const Ray &scattered) const
  {
    return 1; //cos(theta)  fuzz / M_PI;
  }
};
