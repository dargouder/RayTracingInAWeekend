#pragma once
#include "ray.h"
Vec3 RandomInUnitDisk(){
	Vec3 p;
	do {
		p = 2.0 * Vec3(drand48(), drand48(), 0) - Vec3(1,1,0);
	} while (Vec3::dot(p,p) >= 1.0);

	return p;
}
class Camera {
public:

	Camera() {
		lower_left_corner = Vec3(-2.0, -1.0, -1.0);
		horizontal = Vec3(4.0,0.0, 0.0);
		vertical = Vec3(0.0, 2.0, 0.0);
		origin = Vec3(0.0, 0.0, 0.0);
	}


	Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect) {
		Vec3 u,v,w;
		float theta = vfov * M_PI/180;
		float half_height = tan(theta/2);
		float half_width = aspect * half_height;
		origin = lookfrom;
		w = Vec3::unit_vector(lookfrom - lookat);
		u = Vec3::unit_vector(Vec3::cross(vup, w));
		v = Vec3::cross(w,u);

		lower_left_corner = origin - half_width*u - half_height*v-w;
		horizontal = 2*half_width*u;
		vertical  = 2*half_height*v;
		
	}

	/*Ray GetRay(float u, float v){
		return Ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
	}*/

	Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect, 
			float aperture, float focus_dist) {
			lens_radius = aperture / 2;
			float theta = vfov * M_PI/180;
			float half_height = tan(theta/2);
			float half_width = aspect * half_height;
			origin = lookfrom;
			w = Vec3::unit_vector(lookfrom - lookat);
			u = Vec3::unit_vector(Vec3::cross(vup, w));
			v = Vec3::cross(w,u);
			lower_left_corner = origin - half_width /* focus_dist*/ * u - 
								half_height /* focus_dist*/ * v - /*focus_dist **/ w;
			horizontal = 2 * half_width /* focus_dist*/ * u;
			vertical = 2 * half_height /* focus_dist*/ * v;
			
	}

	Ray GetRay(float s, float t){
		Vec3 rd = lens_radius * RandomInUnitDisk();
		Vec3 offset = u * rd.x() + v * rd.y();
		return Ray(origin /* offset*/, lower_left_corner + s*vertical + 
			t* horizontal- origin);
	}

	Vec3 origin;
	Vec3 lower_left_corner;
	Vec3 horizontal;
	Vec3 vertical;
	Vec3 u, v, w;
	float lens_radius;

	/*Vec3 origin;
	Vec3 lower_left_corner;
	Vec3 horizontal;
	Vec3 vertical;*/
};
/*

class Camera {
	public:
		Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect, 
			float aperture, float focus_dist) {
			lens_radius = aperture / 2;
			float theta = vfov * M_PI/180;
			float half_height = tan(theta/2);
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
			Vec3 offset = u * rd.x() + v * rd.y();
			return Ray(origin * offset, lower_left_corner + s*horizontal + 
				t*vertical - origin - offset);
		}

		Vec3 origin;
		Vec3 lower_left_corner;
		Vec3 horizontal;
		Vec3 vertical;
		Vec3 u, v, w;
		float lens_radius;
};
*/