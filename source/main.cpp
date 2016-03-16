#include <iostream>
#include "hitable_list.h"
#include "sphere.h"
#include "float.h"
#include "camera.h"
#include "metal.h"
#include "lambertian.h"
#include "dielectric.h"

Vec3 colour(const Ray& ray, const Hitable& world, int depth)
{
	HitRecord rec;
	if(world.hit(ray, 0.001, FLT_MAX, rec)) 
	{
		Ray scattered;
		Vec3 attenuation;
		if(depth < 50 && rec.mat_ptr->Scatter(ray, rec, attenuation, scattered))
		{
			return attenuation * colour(scattered, world, depth+1);
		}
		 else {
		 	return Vec3(0,0,0);
		 }
	} else {
		Vec3 unit_direction = Vec3::unit_vector(ray.direction());
		float t = 0.5*(unit_direction.y() + 1.0);
		return (1.0-t)*Vec3(1.0, 1.0, 1.0) + t*Vec3(0.5, 0.7, 1.0);
	}
}

void RandomScene(HitableList& list){
	int n = 500;
	list.list.push_back(std::make_unique<Sphere>(Vec3(0,-1000,0), 1000, 
		std::make_unique<Lambertian>(Vec3(0.5,0.5,0.5))));
	
	for(int a = -11; a < 11; a++){
		for(int b = -11; b < 11; b++){
			float choose_mat = drand48();
			Vec3 center( a+ 0.9 * drand48(), 0.2, b + 0.9 * drand48());
			if((center-Vec3(4,0.2,0)).length() > 0.9) {
				if(choose_mat < 0.8) {
					list.list.push_back(std::make_unique<Sphere>(center, 0.2, 
							std::make_unique<Lambertian>(Vec3(drand48() * drand48(),
							 drand48() * drand48(),
							 drand48() * drand48()))));
				} else if(choose_mat < 0.95) {
					list.list.push_back(std::make_unique<Sphere>(center, 0.2,
						 std::make_unique<Metal>(Vec3(0.5 * (1+drand48()), 
						 	0.5 * (1+drand48()),
						 	0.5 * (1+drand48())),
						 0.5 * drand48())));
				} else {
					list.list.push_back(std::make_unique<Sphere>(center, 0.2, 
							std::make_unique<Dielectric>(1.5)));
				}
			}
		}
	}

	list.list.push_back(std::make_unique<Sphere>(Vec3(0, 1, 0), 1.0, std::make_unique<Dielectric>(1.5)));
	list.list.push_back(std::make_unique<Sphere>(Vec3(-4, 1, 0), 1.0, std::make_unique<Lambertian>(Vec3(0.4, 0.2, 0.1))));
	list.list.push_back(std::make_unique<Sphere>(Vec3(4, 1, 0), 1.0, std::make_unique<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));
}

int main() {

	int nx = 600;
	int ny = 400;
	int ns = 20;
	std::cout << "P3\n" << nx << " " << ny << "\n255\n";

	HitableList world;
	Vec3 lookfrom(0,3,1);
	Vec3 lookat(0,0,-1);
	float dist_to_focus = (lookfrom - lookat).length();
	float aperture = 1.0;

	float R = cos(M_PI / 4);
	//world.list.push_back(std::make_unique<Sphere>(Vec3(-R,0,-1), R, std::make_unique<Lambertian>(Vec3(0.0, 0.0, 1))));
	//world.list.push_back(std::make_unique<Sphere>(Vec3(R,0,-1), R, std::make_unique<Lambertian>(Vec3(1, 0.0, 0.0))));

	world.list.push_back(std::make_unique<Sphere>(Vec3(0,0,-1), 0.5, std::make_unique<Lambertian>(Vec3(0.1, 0.2, 0.5))));
	world.list.push_back(std::make_unique<Sphere>(Vec3(0,-100.5,-1), 100, std::make_unique<Lambertian>(Vec3(0.8, 0.8, 0.0))));
	world.list.push_back(std::make_unique<Sphere>(Vec3(-1,0,-1), 0.5, std::make_unique<Metal>(Vec3(0.8, 0.6, 0.2), 1.0)));
	world.list.push_back(std::make_unique<Sphere>(Vec3(-1,0,-1), 0.5, std::make_unique<Dielectric>(1.5)));
	world.list.push_back(std::make_unique<Sphere>(Vec3(-1,0,-1), -0.45, std::make_unique<Dielectric>(1.5)));


	RandomScene(world);
	//world.list.push_back(Vec3(0,0,-1), 0.5, std::make_unique<Lambertian>(Vec3(0.1, 0.2, 0.5)));
	//world.list.push_back(Vec3(0,-100.5,-1), 100, std::make_unique<Lambertian>(Vec3(0.8, 0.8, 0.0)));
	//world.list.push_back(Vec3(1,0,-1), 0.5, std::make_unique<Metal>(Vec3(0.8, 0.6, 0.2)));
	//world.list.push_back(Vec3(-1,0,-1), 0.5, std::make_unique<Dielectric>(1.5));

	//Camera cam(Vec3(-2,2,1), Vec3(0,0,-1), Vec3(0,1,0), 20, float(nx)/float(ny));
	Camera cam(lookfrom, lookat, Vec3(0,1,0), 90, float(nx)/float(ny), aperture, dist_to_focus);

	for(int j = ny-1; j >= 0; j--){
		for(int i = 0; i < nx; i++){
			Vec3 col(0.0, 0.0, 0.0);
			for(int s = 0; s < ns; s++) {
				double u = double(i + drand48()) / double(nx);
				double v = double(j + drand48()) / double(ny);

				Ray r = cam.GetRay(u,v);
			
				Vec3 p = r.point_at_parameter(2.0);
				col += colour(r, world, 0);				
			}

			col /= float(ns);
			col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);

			std::cout << ir << " " << ig << " " << ib << "\n";
		}
	}
}