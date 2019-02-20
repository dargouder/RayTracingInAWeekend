#pragma once

#include "material.h"
#include "onb.h"

class Lambertian : public Material {
public:
    explicit Lambertian(const Vec3 &a) : m_kd(a) {}

    bool fr(const Ray &r_in, const HitRecord &rec, Vec3 &attenuation,
            Ray &scattered, float &pdf) const override {

        Vec3 target = rec.p + rec.normal + RandomInUnitSphere();
        scattered = Ray(rec.p, target - rec.p);
        scattered.direction().make_unit_vector();
        pdf = Vec3::dot(rec.normal, scattered.direction()) / float(M_PI);
        attenuation = m_kd;

        return true;
    }
/*    bool fr(const Ray &r_in, const HitRecord &rec, Vec3 &attenuation,
            Ray &scattered, float &pdf) const override {
        ONB uvw;
        uvw.build_from_w(rec.normal);
        Vec3 direction = uvw.local(CosineSampleHemisphere());
        direction.make_unit_vector();
        scattered = Ray(rec.p, direction );

        pdf = Vec3::dot(uvw.w(), scattered.direction()) / float(M_PI);
        attenuation = m_kd;

        return true;
   } */
    Vec3 Le(float u, float v, const Vec3 &p) const override { return Vec3(0, 0, 0); }

    float ScatteredPdf(const Ray &r_in, const HitRecord &rec, const Ray &scattered) const override {
        float cosine = Vec3::dot(rec.normal, Vec3::unit_vector(scattered.direction()));
        if (cosine < 0) {
           return 0;
        }
        return cosine / float(M_PI);
    }

    Vec3 m_kd;
};
