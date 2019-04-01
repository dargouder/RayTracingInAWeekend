#pragma once
#include "material.h"

class Microfacet : public Material {
 public:
  Microfacet(const Vec3 kd, const Vec3 ks, const float alpha_b)
      : m_kd(kd), m_ks(ks), m_alpha_b(alpha_b) {}

  float FresnelTerm(const float& dot_i_m, const float& refr_index_ratio) const {
    float c = abs(dot_i_m);

    float g_squared = (refr_index_ratio * refr_index_ratio) - 1 + c * c;

    if (g_squared < 0.0f) {
      return 1.0f;
    }

    float g = sqrtf(g_squared);

    float gpc = g + c;
    float gmc = g - c;

    // auto D_m = chi / (M_PI * alpha_b * alpha_b * powf(cos(theta_m), 4)) *
    // exp(-powf(tan(theta_m), 2) / (alpha_b * alpha_b));

    float num = c * (g + c) - 1;
    float den = c * (g - c) + 1;

    float F =
        0.5f * (gmc * gmc / gpc * gpc) * (1.0f + (num * num) / (den * den));

    return F;
  }

  void refl_brdf(const Vec3& incident_ray, const Vec3& micro_normal,
                 const Vec3& macro_normal, Vec3& reflected_ray, float& brdf,
                 float& pdf) {
    //Vec3 reflected_ray =
    //    2.0f * Vec3::dot(incident_ray, micro_normal) * micro_normal -
    //    incident_ray;
    reflected_ray.make_unit_vector();

    Vec3 reflection_half_vec =
        copysign(1, Vec3::dot(incident_ray, macro_normal)) *
        (incident_ray + reflected_ray);
    reflection_half_vec.make_unit_vector();
  }

  void tr_calc(const Vec3& incident_ray, const Vec3& micro_normal,
               const Vec3& macro_normal, const float dot_i_m,
               const float dot_i_n, Vec3& transmitted_ray, float& brdf,
               float& pdf) const {
    float refractive_index_ratio = m_eta_in / m_eta_out;
    float sign_i_m = dot_i_m < 0 ? -1.0f : 1.0f;
    float sign_i_n = dot_i_n < 0 ? -1.0f : 1.0f;
    float sqrt_res =
        sqrtf(1 + refractive_index_ratio * (sign_i_m * sign_i_m - 1));
    transmitted_ray = micro_normal * (refractive_index_ratio * dot_i_m -
                                      sign_i_m * sqrt_res) -
                      (incident_ray * refractive_index_ratio);

    Vec3 trans_half_normal = sign_i_n * (incident_ray + transmitted_ray);

    float dot_i_h = Vec3::dot(incident_ray, trans_half_normal);
    float dot_o_h = Vec3::dot(transmitted_ray, trans_half_normal);
    float dot_o_n = Vec3::dot(transmitted_ray, macro_normal);

/*    brdf = ((abs(dot_i_h) * abs(dot_o_h)) / (abs(dot_i_n) * abs(dot_o_n))) *
           (m_eta_out * m_eta_out) *
           (1 - FresnelTerm(dot_i_h, refractive_index_ratio)) *
           BeckmannShadowMasking(incident_ray, transmitted_ray,
                                 trans_half_normal) *
           BeckmannDistribution(macro_normal, trans_half_normal, theta_m)*/;
    float den_jacobian = m_eta_in * dot_i_h + m_eta_out * dot_o_h;
    float trans_jacobian =
        m_eta_out * m_eta_out * abs(dot_o_h) / (den_jacobian * den_jacobian);
  }

  virtual bool sample_f(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
                  Ray& scattered) const {
    float u = RAND();
    float v = RAND();

    float theta_m = atan(sqrtf(-(m_alpha_b * m_alpha_b * log(1 - u))));
    float phi_m = 2 * PI * v;

    Vec3 microsurface_normal =
        Vec3(cos(theta_m) * cos(phi_m), sin(theta_m) * cos(phi_m), cos(phi_m));
    microsurface_normal.make_unit_vector();

    Vec3 incident_ray = r_in.direction();
    incident_ray.make_unit_vector();

    Vec3 macrosurface_normal = rec.normal;
    macrosurface_normal.make_unit_vector();

    float dot_i_n = Vec3::dot(incident_ray, microsurface_normal);
    float dot_i_m = Vec3::dot(incident_ray, microsurface_normal);

    Vec3 transmitted_ray;
    float tr_brdf, tr_pdf;

    tr_calc(incident_ray, microsurface_normal, macrosurface_normal, dot_i_m,
            dot_i_n, transmitted_ray, tr_brdf, tr_pdf);

    //float ref_id = m_eta_in / m_eta_out;
    //// float dot_incident_microsurface =
    //// float c = acos(abs(Vec3::dot(incident_ray, reflection_half_vec)));

    //// Vec3 o_t = ref_id * c;
    //float G = BeckmannShadowMasking(incident_ray, reflection_half_vec,
    //                                macrosurface_normal);
    //float D =
    //    BeckmannDistribution(reflection_half_vec, macrosurface_normal, theta_m);

    //// float F = FresnelTerm(i, h_r);
    //// float F = schlick(c, 2.0f);

    //float diffuse_brdf =
    //    F * G * D /
    //    (4 * abs(Vec3::dot(incident_ray, macrosurface_normal)) *
    //     abs(Vec3::dot(reflected_ray, macrosurface_normal)));
    //float specular_brdf = abs(Vec3::dot(incident_ray, transmitted_ray)) *
    //                      abs(Vec3::dot(reflected_ray, transmitted_ray)) /

    //                      // float specular_brdf;

    //                      scattered = Ray(rec.p, reflected_ray);
    //attenuation = m_kd;

    return true;
  }

  virtual Vec3 Le(float u, float v, const Vec3& p) const {
    return Vec3(0, 0, 0);
  }

  float BeckmannShadowMasking(Vec3 v, Vec3 m, Vec3 power) const {
    float dot_vm = Vec3::dot(v, m);
    float dot_vn = Vec3::dot(v, power);

    // theta_v is the angle between v and n
    float cos_theta_v = dot_vn / (v.length() * power.length());
    float theta_v = acos(cos_theta_v);
    float a = 1.0f / (m_alpha_b * tan(theta_v));

    float result = 1;
    if (a < 1.6) {
      result =
          (3.535f * a + 2.181f * a * a) / (1.0f + 2.276f * a + 2.577f * a * a);
    }

    return result;
  }

  float BeckmannDistribution(Vec3 m, Vec3 power, float theta_m) const {
    float dot_mn = Vec3::dot(m, power);
    float chi = dot_mn > 0 ? 1.0f : 0.0f;

    float result =
        (chi / (PI * m_alpha_b * m_alpha_b * powf(cos(theta_m), 4))) *
        exp(-powf(tan(theta_m), 2) / (m_alpha_b * m_alpha_b));
    return result;
  }

  float ScatteringPdf(float theta) const { return 0.0f; }

  Vec3 m_kd, m_ks;
  float m_alpha_b, m_ior;
  float m_eta_in, m_eta_out;
};

