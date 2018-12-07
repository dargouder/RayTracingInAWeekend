#pragma once
#include "material.h"

class Microfacet : public Material {
 public:
  Microfacet(const Vec3 kd, const Vec3 ks, const float alpha_b)
      : m_kd(kd), m_ks(ks), m_alpha_b(alpha_b) {}


  float FresnelTerm(Vec3 i, Vec3 m) const
  {
    float c = abs(Vec3::dot(i, m));

    float m_ior = 1.3f, m_ioi = 1.0f;
    float g_squared = (m_ior*m_ior / m_ioi * m_ioi) - 1 + c * c;
    if (g_squared < 0.0f)
    {
      return 1.0f;
    }
    float g = sqrtf(g_squared);

    float gpc = g+c;
    float gmc = g-c;

    float num = (c * gpc - 1) * ( c * gpc - 1);
    float den = (c * gmc - 1) * ( c * gmc - 1);

    float F = 0.5f * ( gmc*gmc / gpc*gpc ) * ( 1.0f + num / den );

    return F;
  }

  virtual bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
                       Ray& scattered) const {
    float u = RAND();
    float v = RAND();

    float theta_m = atan(sqrtf(-(m_alpha_b * m_alpha_b * log(1 - u))));
    float phi_m = 2 * M_PI * v;
    
    Vec3 m =
      Vec3(cos(theta_m) * cos(phi_m), sin(theta_m) * cos(phi_m), cos(phi_m));

    m.make_unit_vector();

    Vec3 i = r_in.direction();
    i.make_unit_vector();

    Vec3 o_r = 2.0f * Vec3::dot(i, m) * m - i;
    o_r.make_unit_vector();
    
    Vec3 n = rec.normal;
    n.make_unit_vector();

    Vec3 h_r = copysign(1, Vec3::dot(i, n)) * (i + o_r);
    h_r.make_unit_vector();

    float G = BeckmannShadowMasking(i, h_r, n);
    float D = BeckmannDistribution(h_r, n, theta_m);
    float c = acos(abs(Vec3::dot(i, h_r)));
    //float F = FresnelTerm(i, h_r);
    float F = schlick(c, 2.0f);

    float brdf = F * G * D /  (4 * abs(Vec3::dot(i, n)) * abs(Vec3::dot(o_r, n)));

    scattered = Ray(rec.p, o_r);
    attenuation = m_kd;

    return true;
  }

  virtual Vec3 Le(float u, float v, const Vec3& p) const {
    return Vec3(0, 0, 0);
  }

  float BeckmannShadowMasking(Vec3 v, Vec3 m, Vec3 n) const {
    float dot_vm = Vec3::dot(v, m);
    float dot_vn = Vec3::dot(v, n);

    // theta_v is the angle between v and n
    float cos_theta_v = dot_vn / (v.length() * n.length());
    float theta_v = acos(cos_theta_v);
    float a = 1.0f / (m_alpha_b * tan(theta_v));
       

    float result = 1;
    if (a < 1.6)
    {
      result = (3.535f * a + 2.181f * a * a) / (1.0f + 2.276f * a + 2.577f * a * a);
    }

    return result;
  }

  float BeckmannDistribution(Vec3 m, Vec3 n, float theta_m) const {


    float dot_mn = Vec3::dot(m, n);
    float chi = dot_mn > 0 ? 1 : 0;

     float result = (chi / (M_PI * m_alpha_b * m_alpha_b * powf(cos(theta_m), 4))) *
               exp(-powf(tan(theta_m), 2) / (m_alpha_b * m_alpha_b));
     return result;
  }

  float Pdf(float theta) const
  {
    return 0.0f;
  }

  Vec3 m_kd, m_ks;
  float m_alpha_b,  m_ior;
};
