#pragma once

#include "material.h"
#include "onb.h"

class Lambertian : public Material {
public:
    explicit Lambertian(const Vec3 &a) : m_kd(a) {}


    bool fr(const Ray &r_in, const HitRecord &rec, Vec3 &attenuation,
            Ray &scattered, float &pdf) const override 
    {
      ONB onb;
      onb.build_from_w(rec.normal);

      float u = RAND();
      float v = RAND();
      Vec3 generatedCosinePoint = CosineSampleHemispherePhong(u, v);
      Vec3 transformedPoint = onb.local(generatedCosinePoint);
      Vec3 target = Vec3::unit_vector(transformedPoint);

      scattered = Ray(rec.p, target);
      attenuation = m_kd / M_PI;
      pdf = Vec3::dot(rec.normal, scattered.direction()) / M_PI;
      return true;
    }

    Vec3 Le(float u, float v, const Vec3 &p) const override { return Vec3(0, 0, 0); }

    float ScatteredPdf(const Ray &r_in, const HitRecord &rec, const Ray &scattered) const override 
    {
        float cosine = Vec3::dot(rec.normal, scattered.direction());
        if (cosine < 0) {
           return 0;
        }
        return cosine / float(M_PI);
    }

    Vec3 m_kd;
};
