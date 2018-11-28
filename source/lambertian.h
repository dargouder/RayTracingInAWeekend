#pragma once
#include "material.h"

class Lambertian : public Material {

public:
	Lambertian(const Vec3& a) : albedo(a) {

	}

	virtual bool Scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		Vec3 target = rec.p + rec.normal + RandomInUnitSphere();
		scattered = Ray(rec.p, target - rec.p);
		attenuation = albedo;

		return true;
	}

  virtual Vec3 emitted(float u, float v, const Vec3& p) const { return Vec3(0,0,0); }

	Vec3 albedo;

};