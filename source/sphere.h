#pragma once

#include "hitable.h"
#include <memory>

class Sphere : public Hitable {
public:
	Vec3 center;
	float radius;
	std::unique_ptr<Material> material;
	Sphere() {}
	Sphere(Vec3 cen, float r, std::unique_ptr<Material> mat) : center(cen), radius(r), material(std::move(mat)) {}

	bool hit(const Ray& r, float tmin, float tmax, HitRecord& rec) const {
		Vec3 oc = r.origin() - center;
		float a = Vec3::dot(r.direction(), r.direction());
		float b = Vec3::dot(oc, r.direction());
		float c = Vec3::dot(oc,oc) - radius*radius;
		float discriminant = b*b - a*c;

		if (discriminant > 0) {
			float temp = (-b - sqrt(discriminant)) / a;
			if (temp < tmax && temp > tmin) {
				rec.t = temp;
				rec.p = r.point_at_parameter(rec.t);
				rec.normal = (rec.p - center) / radius;
				rec.mat_ptr = material.get();
				return true;
			}
		
			temp = (-b + sqrt(discriminant))/a;
			if(temp < tmax && temp > tmin) {
				rec.t = temp;
				rec.p = r.point_at_parameter(rec.t);
				rec.normal = (rec.p - center) / radius;
				rec.mat_ptr = material.get();
				return true;
			}
		}
		return false;
	}


};