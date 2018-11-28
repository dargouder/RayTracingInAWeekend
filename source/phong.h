#pragma once
#include "material.h"

class Phong : public Material {
 public:
  Phong(const Vec3 kd, const Vec3 ks, const float n) : m_kd(kd), m_ks(ks), 
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