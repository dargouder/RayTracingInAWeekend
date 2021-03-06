#pragma once
#include "material.h"

class Phong : public Material {
 public:
  Phong(const Vec3& a, const Vec3& ks) : m_kd(a), m_ks(ks) { n = 30; }

  bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
          Ray& scattered, float& pdf) const {
    float epsilon = RAND();

    float u = RAND();
    float v = RAND();

    float diffuse_lum = m_kd.getLuminance();
    float specular_lum = m_ks.getLuminance();

    Vec3 perfectReflection =
        reflect(Vec3::unit_vector(r_in.direction()), rec.normal);
    perfectReflection.make_unit_vector();

    Vec3 scattered_dir;
    Vec3 scattered_target;
    float alpha;
    pdf = 0.0f;
    if (epsilon < diffuse_lum) {
      // Diffuse component
      float sin_theta = sqrtf(1 - u);
      float phi = 2 * M_PI * v;
      float cos_phi = std::cos(phi);
      float sin_phi = std::sin(phi);

      scattered_dir =
          rec.normal + Vec3(sin_theta * cos_phi, sin_theta * sin_phi, sqrt(u));
      scattered_dir.make_unit_vector();

      pdf = Vec3::dot(scattered_dir, rec.normal) / M_PI;
    } else if (epsilon <= diffuse_lum + specular_lum) {
      // Specular component
      float sin_alpha = std::sqrt(1 - std::pow(u, 2.0f / (n + 1.0f)));
      float cos_alpha = std::pow(u, 1 / (n + 1));
      float phi = 2 * M_PI * v;
      float cos_phi = std::cos(phi);
      float sin_phi = std::sin(phi);

      scattered_target =
          Vec3(sin_alpha * cos_phi, sin_alpha * sin_phi, cos_alpha);
      scattered_dir = rec.normal + scattered_target;
      scattered_dir.make_unit_vector();

      float alpha = Vec3::dot(scattered_dir, perfectReflection);
      if (alpha > 0) {
        pdf = (n + 1.0f) / (2.0f * M_PI) * powf(alpha, n);
      }
    } else {
      attenuation = Vec3(0.0f, 0.0f, 0.0f);
      return false;
    }

    scattered = Ray(rec.p, scattered_dir);
    alpha = acos(Vec3::dot(perfectReflection, scattered_dir));

    alpha = fmax(alpha, 2 * M_PI);

    attenuation =
        m_kd / M_PI + m_ks * (n + 2.0f) / (2.0f * M_PI) * pow(cos(alpha), n);

    return true;
  }

  Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0, 0, 0); }

  float ScatteredPdf(const Ray &r_in, const HitRecord &rec, const Ray &scattered) const {
      return 0;//cos(theta) / M_PI;
  }

  Vec3 m_kd, m_ks;
  float n;
};
