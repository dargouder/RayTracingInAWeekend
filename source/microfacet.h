#pragma once
#include "material.h"

void BeckmannShadowMasking(Vec3 v, Vec3 m, Vec3 n)
{
    float dot_vm = Vec3::dot(v, m);
    float dot_vn = Vec3::dot(v, n);

    // theta_v is the angle between v and n
    float cos_theta_v = dot_vn / (v.length() * n.length());
    float theta_v = acos(cos_theta_v);
    float a = 1.0f / (alpha_b * tan(theta_v));
}

void BeckmannDistribution(Vec3 n, float u, float v, float alpha_b)
{
    float theta_m = arctan(sqrtf(-(alpha_b*alpha_b * log( 1 - u))));
    float phi_m = 2 * M_PI * v;

    Vec3 m = Vec3( cos(theta_m) * cos(phi_m), sin(theta_m * cos(phi_m), cos(phi)));

    float dot_mn = Vec3::dot(m, n);
    float chi = dot_mn > 0 ? 1 : 0;

    auto D_m = chi / (M_PI * alpha_b * alpha_b * powf(cos(theta_m), 4)) * exp(-powf(tan(theta_m), 2) / (alpha_b * alpha_b));

    // G(i,o,m) = G_1 (i,m) * G_1(o,m)
    // G_1 is the shadow masking function
}

class Microfacet : public Material {
 public:
  Microfacet(const Vec3 kd, const Vec3 ks, const float n) : m_kd(kd), m_ks(ks),
	  m_n(n) {}

  virtual bool Scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
                       Ray& scattered) const {

    return true;
  }

  virtual Vec3 emitted(float u, float v, const Vec3& p) const {
    return Vec3(0, 0, 0);
  }

  Vec3 m_kd, m_ks;
  float m_n;

};