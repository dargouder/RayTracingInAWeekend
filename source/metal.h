#pragma once
#include "material.h"

class Metal : public Material {
public:
	Vec3 albedo;
	float fuzz;
	Metal(const Vec3& a, float f) : albedo(a) {
		if(f < 1){
			fuzz = f;
		} else {
			fuzz = 1;
		}
	}

	virtual bool Scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		Vec3 reflected = reflect(Vec3::unit_vector(r_in.direction()), rec.normal);
		scattered = Ray(rec.p, reflected + fuzz * RandomInUnitSphere());
		attenuation = albedo;

		return (Vec3::dot(scattered.direction(), rec.normal) > 0);
	}
};