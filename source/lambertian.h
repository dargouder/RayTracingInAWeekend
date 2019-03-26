#pragma once

#include "material.h"
#include "onb.h"

class Lambertian : public Material {
 public:
  explicit Lambertian(const Vec3 &a) : m_kd(a) {}

  virtual Vec3 f(const Vec3 &wo, const Vec3 &wi) const { return Vec3(0, 0, 0); }

  bool sample_f(const Ray &r_in, const HitRecord &rec, Vec3 &attenuation,
                Ray &scattered, float &pdf) const override {
    ONB onb, shadingONB;
    shadingONB.branchlessONB(Vec3(0, 0, 1));
    onb.branchlessONB(rec.normal);

	Vec3 wo = shadingONB.local(-r_in.direction());

    float u = RAND();
    float v = RAND();
    Vec3 wi = CosineSampleHemispherePBRT(u, v);
	if ( wo.z < 0 )
	{
		wi.z *= -1;
	}

  
    Vec3 target = Vec3::unit_vector(onb.local(wi));

    scattered = Ray(rec.p, target);
    attenuation = m_kd / PI;
    attenuation = m_kd;
    pdf = std::abs(Vec3::dot(rec.normal, scattered.direction())) / PI;
    float untransformedPdf = std::fabsf(wi.z) / PI;
    return true;
  }

  Vec3 Le(float u, float v, const Vec3 &p) const override {
    return Vec3(0, 0, 0);
  }

  float ScatteredPdf(const Ray &r_in, const HitRecord &rec,
                     const Ray &scattered) const override {
    float cosine = Vec3::dot(rec.normal, scattered.direction());
    if (cosine < 0) {
      return 0;
    }
    return cosine / float(PI);
  }

  float pdf(const Vec3 &wo, const Vec3 &wi) const override {
    return Vec3::dot(Vec3(0, 0, 1), wi) / PI;
  }

  Vec3 m_kd;
};
