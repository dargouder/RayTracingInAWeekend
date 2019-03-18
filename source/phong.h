#pragma once
#include "material.h"
#include "onb.h"

class Phong : public Material {
 public:
  Phong(const Vec3& a, const Vec3& ks) : m_kd(a), m_ks(ks) { n = 1;
      diffuseRatio = m_kd.getLuminance()/(m_kd.getLuminance() + m_ks.getLuminance());
      const Vec3 combined = m_kd + m_ks;

      assert(combined.r() <= 1);
      assert(combined.g() <= 1);
      assert(combined.b() <= 1);
  }

  bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
          Ray& scattered, float& pdf) const {

    ONB onb;
    onb.build_from_w(rec.normal);

    // compute single luminance value to see how much we will reflect
    float diffuse_lum = m_kd.getLuminance();
    float specular_lum = m_ks.getLuminance();

    assert(diffuseRatio != 0);
    assert(diffuse_lum + specular_lum <= 1);

    // calculate the perfect specular reflection
    Vec3 perfectReflection =
        reflect(Vec3::unit_vector(r_in.direction()), rec.normal);
    perfectReflection.make_unit_vector();

    Vec3 scatteredDir;
    Vec3 scattered_target;

    pdf = 0.0f;

    // generate a random number to see whether we will reflect with specular or diffuse
    float reflectionDecision = RAND();

    if (reflectionDecision < diffuseRatio || diffuseRatio == 1.0f) {
      // calculating only diffuse component
      float u = RAND();
      float v = RAND();

      scatteredDir = onb.local(CosineSampleHemispherePhong(u, v));
      scatteredDir.make_unit_vector();

      pdf = Vec3::dot(scatteredDir, rec.normal) / M_PI;
    } else  {
      // calculating only specular component
      float u = RAND();
      float v = RAND();
      float sin_alpha = std::sqrt(1 - std::pow(u, 2.0f / (n + 1.0f)));
      float cos_alpha = std::pow(u, 1 / (n + 1));
      float phi = 2 * M_PI * v;
      float cos_phi = std::cos(phi);
      float sin_phi = std::sin(phi);

      scatteredDir =
          onb.local(Vec3(sin_alpha * cos_phi, sin_alpha * sin_phi, cos_alpha));
      scatteredDir.make_unit_vector();

      float alpha = std::max(std::min(1.0f, Vec3::dot(scatteredDir, perfectReflection)), 0.0f);
      if (alpha >= 0) {
        pdf = ((n + 2.0f) / (2.0f * M_PI)) * powf(alpha, n);
      }
    }

    // the alpha for the specular reflection
    float alpha;
    scattered = Ray(rec.p, scatteredDir);
    alpha = Vec3::dot(perfectReflection, perfectReflection);

    alpha = alpha >  M_PI / 2.0f ? M_PI / 2.0f : alpha;

    // calculate the fr
    attenuation =
        m_kd / M_PI + m_ks * (n + 2.0f) * (1.0f / (2.0f * M_PI)) * pow(alpha, n);

    return true;
  }

  Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0, 0, 0); }

  float ScatteredPdf(const Ray& r_in, const HitRecord& rec,
                     const Ray& scattered) const {
    return 0;  // cos(theta) / M_PI;
  }

  Vec3 m_kd, m_ks;
  float n;
  float diffuseRatio;
};
