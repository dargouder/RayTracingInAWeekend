#pragma once
#include "transform.h"
//
//class Camera{
//public:
//	int height;
//	int width;
//	float fov;
//	Transform camera_to_world;
//	Transform raster_to_camera;
//	Transform screen_to_raster;
//	Transform raster_to_screen;
//	Transform camera_to_screen;
//
//	float aspect_ratio;
//	float inv_aspect_ratio;
//	float scale;
//
//	Vec3 position;
//
//
//	Vec3 pos, top, dx, dy;
//	~Camera() {
//
//	}
//
//
//
//	Camera(const int p_width, const int p_height, const float p_fov, const Vec3& p_position) : width(p_width), height(p_height) {
//		camera_to_world = camera_to_world.Translate(p_position);
//
//		aspect_ratio = static_cast<float>(p_width) / static_cast<float>(p_height);
//		camera_to_screen = camera_to_screen.Perspective(p_fov, 1e-2f, 1000);
//
//
//		float screen_window_max_x;
//		float screen_window_max_y;
//		float screen_window_min_x;
//		float screen_window_min_y;
//		if (aspect_ratio > 1.f) {
//			screen_window_min_x = -aspect_ratio;
//			screen_window_max_x = aspect_ratio;
//			screen_window_min_y = -1.f;
//			screen_window_max_y = 1.f;
//		}
//		else {
//			screen_window_min_x = -1.f;
//			screen_window_max_x = 1.f;
//			screen_window_min_y = -1.f / aspect_ratio;
//			screen_window_max_y = 1.f / aspect_ratio;
//		}
//		
//		screen_to_raster =
//			screen_to_raster.Scale(p_width, p_height, 1) *
//			screen_to_raster.Scale(1 / (screen_window_max_x - screen_window_min_x),
//				1 / (screen_window_min_y - screen_window_max_y), 1);
//
//		raster_to_screen = Transform::Inverse(screen_to_raster);
//
//		raster_to_camera = Transform::Inverse(camera_to_screen) * raster_to_screen;
//		//camera_to_world = camera_to_world.LookAt(Vec3(0, 2, 0), Vec3(0, 0, -1), Vec3(0, 1, 0));
//		scale = tan(p_fov * (M_PI / 180.0f) * 0.5f);
//		
//		inv_aspect_ratio = 1 / aspect_ratio;
//	}
//
//	Ray GenerateRay(const double& p_x, const double& p_y) const  {
//
//		Vec3 p_film = Vec3(p_x, p_y, 0);
//		Vec3 p_camera = raster_to_camera(p_film);
//		Ray r(Vec3(0, 0, 0), Vec3::unit_vector(Vec3(p_camera)));
//		//float x = (2 * (p_x) / static_cast<float>(width) - 1) * scale;
//		//float y = (1 - 2 * (p_y) / static_cast<float>(height)) * scale * inv_aspect_ratio;
//
//		//Vec3 direction;
//		//r.A = Vec3(0, 0, 0);
//		//r.B = Vec3::unit_vector(Vec3(x, y, -1));
//		r = camera_to_world(r);
//
//
//		return r;
//	}
//
//	void Update() {
//		position = camera_to_world(position);
//	}
//
//};
//
//class Camera {
//public:
//
//	Camera() {
//		lower_left_corner = Vec3(-2.0, -1.0, -1.0);
//		horizontal = Vec3(4.0,0.0, 0.0);
//		vertical = Vec3(0.0, 2.0, 0.0);
//		origin = Vec3(0.0, 0.0, 0.0);
//	}
//
//
//	/*Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect) {
//		Vec3 u,v,w;
//		float theta = vfov * M_PI/180;
//		float half_height = tan(theta/2);
//		float half_width = aspect * half_height;
//		origin = lookfrom;
//		w = Vec3::unit_vector(lookfrom - lookat);
//		u = Vec3::unit_vector(Vec3::cross(vup, w));
//		v = Vec3::cross(w,u);
//
//		lower_left_corner = origin - half_width*u - half_height*v-w;
//		horizontal = 2*half_width*u;
//		vertical  = 2*half_height*v;
//		
//	}*/
//
//	/*Ray GetRay(float u, float v){
//		return Ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
//	}*/
//
//	Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect) {
//
//			float theta = vfov * M_PI/180;
//			
//			float half_height = tan(theta/2);
//			scale = half_height;
//			float half_width = aspect * half_height;
//			origin = lookfrom;
//			w = Vec3::unit_vector(lookfrom - lookat);
//			u = Vec3::unit_vector(Vec3::cross(vup, w));
//			v = Vec3::cross(w,u);
//			lower_left_corner = origin - half_width * u - half_height * v -  w;
//			horizontal = 2.0 * half_width * u;
//			vertical = 2.0 * half_height * v;
//			
//	}
//
//	Ray GetRay(float s, float t){
//		float x = 2 * s * scale;
//		float y = (1 - 2 * t)*scale *  float(1 / (1.5));
//		return Ray(origin ,Vec3(x,y,-1));
//	}
//	float scale;
//	Vec3 origin;
//	Vec3 lower_left_corner;
//	Vec3 horizontal;
//	Vec3 vertical;
//	Vec3 u, v, w;
//	float lens_radius;
//
//	/*Vec3 origin;
//	Vec3 lower_left_corner;
//	Vec3 horizontal;
//	Vec3 vertical;*/
//};


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
			Vec3 offset = u * rd.x() + v * rd.y();
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
