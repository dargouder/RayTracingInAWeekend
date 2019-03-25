#pragma once

#include "vec3.h"

class Ray 
{
public:

	Ray() {}
	Ray(const Vec3& a, const Vec3& b) { A = a; B = b;
    B.make_unit_vector(); }
	Vec3 origin() const {return A;}
	Vec3 direction() const { return B;}
	Vec3 point_at_parameter(float t) const {
		return A + t*B; 
	}

	Vec3 A;
	Vec3 B;
};

inline Vec3 OffsetRayOrigin(const Vec3 &p, const Vec3 &pError, const Vec3 &normal, const Vec3 &w)
{

}