#include <iostream>
#include "hitable_list.h"
#include "sphere.h"
#include "float.h"
#include "camera.h"
#include "metal.h"
#include "lambertian.h"
#include "dielectric.h"
#include <fstream>

Vec3 colour(const Ray& ray, const Hitable& world, int depth)
{
	HitRecord rec;
	if (world.hit(ray, 0.001, FLT_MAX, rec))
	{
		Ray scattered;
		Vec3 attenuation;
		if (depth < 50 && rec.mat_ptr->Scatter(ray, rec, attenuation, scattered))
		{
			return attenuation * colour(scattered, world, depth + 1);
		}

		return Vec3(0, 0, 0);

	}
	else {
		Vec3 unit_direction = Vec3::unit_vector(ray.direction());
		float t = 0.5*(unit_direction.y() + 1.0);
		return (1.0 - t)*Vec3(1.0, 1.0, 1.0) + t*Vec3(0.5, 0.7, 1.0);
	}
}

void RandomScene(HitableList& list) {
	int n = 500;
	list.list.push_back(std::make_unique<Sphere>(Vec3(0, -1000, 0), 1000,
		std::make_unique<Lambertian>(Vec3(0.5, 0.5, 0.5))));

	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			float choose_mat = RAND();
			Vec3 center(a + 0.9 * RAND(), 0.2, b + 0.9 * RAND());
			if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
				if (choose_mat < 0.8) {
					list.list.push_back(std::make_unique<Sphere>(center, 0.2, std::make_unique<Lambertian>(Vec3(RAND() * RAND(), RAND() * RAND(), RAND() * RAND()))));
				}
				else if (choose_mat < 0.95) {
					list.list.push_back(std::make_unique<Sphere>(center, 0.2, std::make_unique<Metal>(Vec3(0.5 * (1 + RAND()), 0.5 * (1 + RAND()), 0.5 * (1 + RAND())), 0.5 * RAND())));
				}
				else {
					list.list.push_back(std::make_unique<Sphere>(center, 0.2, std::make_unique<Dielectric>(1.5)));
				}
			}
		}
	}

	list.list.push_back(std::make_unique<Sphere>(Vec3(0, 1, 0), 1.0, std::make_unique<Dielectric>(1.5)));
	list.list.push_back(std::make_unique<Sphere>(Vec3(-4, 1, 0), 1.0, std::make_unique<Lambertian>(Vec3(0.4, 0.2, 0.1))));
	list.list.push_back(std::make_unique<Sphere>(Vec3(4, 1, 0), 1.0, std::make_unique<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));
}

int main() {
	std::ofstream os;
	os.open("image6.ppm", std::ios::binary);
	int nx = 1920;
	int ny = 1080;
	int ns = 100;
	//os << "P3\n" << nx << " " << ny << "\n255\n";

	os << "P6" << std::endl;
	os << nx << " " << ny << std::endl;
	os << "255" << std::endl;

	HitableList world;

	world.list.push_back(std::make_unique<Sphere>(Vec3(-1.5, 0.0, -1), 0.5, std::make_unique<Lambertian>(Vec3(0.1, 0.2, 0.5))));
	//world.list.push_back(std::make_unique<Sphere>(Vec3(1, -100.5, -1), 100, std::make_unique<Lambertian>(Vec3(0.8, 0.8, 0.0))));
	world.list.push_back(std::make_unique<Sphere>(Vec3(1.5, 0.0, -1), 0.5, std::make_unique<Metal>(Vec3(0.8, 0.6, 0.2), 1.0)));
	world.list.push_back(std::make_unique<Sphere>(Vec3(0.0, 0.0, -1), 0.5, std::make_unique<Dielectric>(1.5)));
	//world.list.push_back(std::make_unique<Sphere>(Vec3(-1, 0.0, -5), 3.45, std::make_unique<Lambertian>(Vec3(0.6, 0.4, 5.5))));


	RandomScene(world);

	Camera cam(nx, ny, 90, Vec3(-6, 3, -5));

	//for (int j = ny - 1; j >= 0; j--) {
		//for (int i = 0; i < nx; i++) {

	for (int j = 0; j < ny; j++) {
		for (int i = nx - 1; i >= 0; i--) {
			Vec3 col(0.0, 0.0, 0.0);
			for (int s = 0; s < ns; s++) {
				double u = double(i + RAND());
				double v = double(j + RAND());

				Ray r = cam.GenerateRay(u, v);

				Vec3 p = r.point_at_parameter(2.0);
				col += colour(r, world, 0);
			}

			col /= float(ns);
			col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);

			unsigned char red = static_cast<unsigned char>(std::min(1.f, col[0]) * 255);
			unsigned char green = static_cast<unsigned char>(std::min(1.f, col[1]) * 255);
			unsigned char blue = static_cast<unsigned char>(std::min(1.f, col[2]) * 255);

			//os << ir << " " << ig << " " << ib << "\n";
			os << red << green << blue;
		}
	}

	os.close();
	return 0;
}