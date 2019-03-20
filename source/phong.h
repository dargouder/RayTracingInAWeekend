#pragma once
#include "material.h"
#include "onb.h"

#include <cassert>

class Phong : public Material {
 public:
  Phong(const Vec3& a, const Vec3& ks) : m_kd(a), m_ks(ks) {
    power = 50;
    diffuseRatio =
        m_kd.getLuminance() / (m_kd.getLuminance() + m_ks.getLuminance());
    const Vec3 combined = m_kd + m_ks;

    assert(combined.r() <= 1);
    assert(combined.g() <= 1);
    assert(combined.b() <= 1);
  }

  virtual Vec3 f(const Vec3& wo, const Vec3& wi) const {
    return Vec3(0.0f, 0.0f, 0.0f);
  }

  bool sample_f(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
                Ray& scattered, float& pdf) const {
    ONB onb;
    onb.branchlessONB(rec.normal);

    // compute single luminance value to see how much we will reflect
    float diffuse_lum = m_kd.getLuminance();
    float specular_lum = m_ks.getLuminance();

    // calculate the perfect specular reflection
    Vec3 perfectReflection =
        reflect(Vec3::unit_vector(r_in.direction()), rec.normal);
    perfectReflection.make_unit_vector();

    Vec3 scatteredDir;
    Vec3 scattered_target;

    pdf = 0.0f;

    // generate a random number to see whether we will reflect with specular or
    // diffuse
    float reflectionDecision = RAND();
    attenuation = m_kd / M_PI;

    if (reflectionDecision < diffuseRatio || diffuseRatio == 1.0f) {
      // calculating only diffuse component
      float u = RAND();
      float v = RAND();

      scatteredDir = onb.local(CosineSampleHemispherePhong(u, v));
      scatteredDir.make_unit_vector();

      pdf = Vec3::dot(scatteredDir, rec.normal) * INV_PI;
    } else {
      // calculating only specular component
      float u = RAND();
      float v = RAND();

      float sin_alpha = std::sqrt(1 - std::powf(u, 2.0f / (power + 1.0f)));
      float cos_alpha = std::powf(u, 1.0f / (power + 1));
      float phi = 2 * M_PI * v;
      float cos_phi = std::cos(phi);
      float sin_phi = std::sin(phi);

      scatteredDir =
          onb.local(Vec3(sin_alpha * cos_phi, sin_alpha * sin_phi, cos_alpha));
      scatteredDir.make_unit_vector();

      float cos_theta = Vec3::dot(scatteredDir, perfectReflection);

      cos_theta = clamp(cos_theta, M_PI / 2.0f, -M_PI / 2.0f);

      pdf = (power + 2.0f) * INV_2PI * powf(cos_theta, power + 1.0f);
      attenuation += m_ks * (power + 2.0f) * INV_2PI * pow(cos_theta, power);
    }
    pdf = std::max(pdf, 0.00001f);
    scattered = Ray(rec.p, scatteredDir);

    return true;
  }

  Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0, 0, 0); }

  float ScatteredPdf(const Ray& r_in, const HitRecord& rec,
                     const Ray& scattered) const {
    return 0;  // cos(theta) / M_PI;
  }

  float pdf(const Vec3& wo, const Vec3& wi) const {
    // const float diffusePDF = Vec3::dot(wi, normal) / M_PI;
    return 0.0f;
  }

  Vec3 m_kd, m_ks;
  float power;
  float diffuseRatio;
};
