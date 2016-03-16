#pragma once

#include "material.h"

class Dielectric : public Material {
public:
	float refl_idx;

	Dielectric(float ri) : refl_idx(ri) {

	}

	virtual bool Scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		Vec3 outward_normal;
		Vec3 reflected = reflect(r_in.direction(), rec.normal);
		float ni_over_nt;
		attenuation = Vec3(1.0,1.0,1.0);
		Vec3 refracted;
		float reflect_prob;
		float cosine;
		if(Vec3::dot(r_in.direction(), rec.normal) > 0){
			outward_normal = -rec.normal;
			ni_over_nt = refl_idx;
			cosine = refl_idx * Vec3::dot(r_in.direction(), rec.normal) / r_in.direction().length();
		} else 
		{
			outward_normal = rec.normal;
			ni_over_nt = 1.0 / refl_idx;
			cosine = -Vec3::dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}

		if(refract(r_in.direction(), outward_normal, ni_over_nt, refracted)){
			reflect_prob = schlick(cosine, refl_idx);
		} else {
			scattered = Ray(rec.p, reflected);
			reflect_prob = 1.0;
		}

		if(drand48() < reflect_prob)
		{
			scattered = Ray(rec.p, reflected);
		} else {
			scattered = Ray(rec.p, refracted);
		}
		return true;
	}
};