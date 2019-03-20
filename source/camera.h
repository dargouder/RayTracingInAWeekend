#pragma once
#include "transform.h"

class Camera {
	public:
		Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect,
			float aperture, float focus_dist) {
			lens_radius = aperture / 2.0;
			float theta = vfov * M_PI/180.0;
			float half_height = tan(theta/2.0);
			float half_width = aspect * half_height;
			origin = lookfrom;
			w = Vec3::unit_vector(lookfrom - lookat);
			u = Vec3::unit_vector(Vec3::cross(vup, w));
			v = Vec3::cross(w,u);
			lower_left_corner = origin - half_width * focus_dist * u -
								half_height * focus_dist * v - focus_dist * w;
			horizontal = 2 * half_width * focus_dist * u;
			vertical = 2 * half_height * focus_dist * v;

		}

		Ray GetRay(float s, float t){
			Vec3 rd = lens_radius * RandomInUnitDisk();
			Vec3 offset = u * rd.x + v * rd.y;
			return Ray(origin + offset,
				lower_left_corner + s*horizontal +
				t*vertical - origin - offset);
		}

		Vec3 origin;
		Vec3 lower_left_corner;
		Vec3 horizontal;
		Vec3 vertical;
		Vec3 u, v, w;
		float lens_radius;
};
