#pragma once
#include "material.h"

class Lambertian : public Material {

public:
	Lambertian(const Vec3& a) : m_kd(a) {

	}

	bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered, float& pdf) const {
		Vec3 target = rec.p + rec.normal + RandomInUnitSphere();
		scattered = Ray(rec.p, target - rec.p);
    scattered.direction().make_unit_vector();
    float cos_theta = Vec3::dot(scattered.direction(), rec.normal);
    pdf = Pdf(cos_theta);
		attenuation = m_kd / M_PI;

		return true;
	}

  Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0,0,0); }

  float Pdf(float theta) const
  {
    return theta / M_PI;
  }

	Vec3 m_kd;

};