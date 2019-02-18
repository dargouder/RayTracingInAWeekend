#pragma once
#include "material.h"
#include "onb.h"

class Lambertian : public Material {
 public:
  Lambertian(const Vec3& a) : m_kd(a) {}

  bool fr(const Ray& r_in, const HitRecord& rec, Vec3& attenuation,
          Ray& scattered, float& pdf) const {
    ONB uvw;
    uvw.build_from_w(rec.normal);
    Vec3 direction = uvw.local(CosineSampleHemisphere());
    direction.make_unit_vector();
    //Vec3 target = rec.normal + CosineSampleHemisphere ();// RandomInUnitSphere();
    scattered = Ray(rec.p, direction);
    //scattered.direction().make_unit_vector();
    float cos_theta = Vec3::dot(scattered.direction(), rec.normal);
    pdf = Vec3::dot(uvw.w(), scattered.direction()) / M_PI;
    attenuation = m_kd;

    return true;
  }

  Vec3 Le(float u, float v, const Vec3& p) const { return Vec3(0, 0, 0); }

  float Pdf(float theta) const { return theta / M_PI; }

  Vec3 m_kd;
};
