#include "camera.h"
#include "dielectric.h"
#include "diffuselight.h"
#include "float.h"
#include "hitable_list.h"
#include "lambertian.h"
#include "metal.h"
#include "microfacet.h"
#include "phong.h"
#include "quad.h"
#include "rects.h"
#include "sphere.h"

#include <fstream>
#include <iostream>

Vec3 colourRecursive(const Ray &ray, const Hitable &world, int depth) {
  HitRecord rec;
  if (world.hit(ray, 0.001, FLT_MAX, rec)) {
    Ray scattered;
    Vec3 attenuation;
    float u = 0.0f, v = 0.0f;
    rec.normal.make_unit_vector();
    Vec3 emitted = rec.mat_ptr->Le(0, 0, rec.p);
    float pdf;
    if (depth < 10 && rec.mat_ptr->fr(ray, rec, attenuation, scattered, pdf)) {

      Vec3 f = attenuation / pdf;
      assert(pdf > 0);
      float cos_theta = Vec3::dot(rec.normal, scattered.direction());
      return emitted + f * cos_theta * colourRecursive(scattered, world, depth+1);
    } else {
      return emitted;
    }
  } else {
    Vec3 unit_direction = Vec3::unit_vector(ray.direction());
    float t = 0.5f * (unit_direction.y() + 1.0f);
    Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
    return col;
  }
}

Vec3 colourRecursiveOld(const Ray &ray, const Hitable &world, int depth) {
  HitRecord rec;
  if (world.hit(ray, 0.001, FLT_MAX, rec)) {
    Ray scattered;
    Vec3 attenuation;
    float u = 0.0f, v = 0.0f;
    rec.normal.make_unit_vector();
    Vec3 emitted = rec.mat_ptr->Le(0, 0, rec.p);
    float pdf;
    if (depth < 10 && rec.mat_ptr->fr(ray, rec, attenuation, scattered, pdf)) {
      return emitted +
             attenuation * colourRecursive(scattered, world, depth + 1);
    } else {
      return emitted;
    }
  } else {
    Vec3 unit_direction = Vec3::unit_vector(ray.direction());
    float t = 0.5f * (unit_direction.y() + 1.0f);
    Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
    return col;
  }
}

void CornellBox(HitableList &list, Hitable **light) {
  // floor
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-10.0f, 0.0f, -10.0f), Vec3(-10.0f, 0.0f, 10.0f),
      Vec3(10.0f, 0.0f, 10.0f), Vec3(10.0f, 0.0f, -10.0f),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
 // list.list.push_back(std::make_unique<Sphere>(
 //     Vec3(0.0f, 0.5f, 0.5f), 0.5f,
 //     std::make_unique<Lambertian>(Vec3(0.0f, 0.9f, 0.0f))));
   list.list.push_back(std::make_unique<Sphere>(Vec3(0.0f, 0.5f, 0.5f), 0.5f,
   std::make_unique<Phong>(Vec3(0.0f, 0.2f, 0.0f), Vec3(0.0f, 0.1f, 0.0f))));
}

int main() {
  std::ofstream os;
  os.open("directLighting.ppm", std::ios::binary);
  const int nx = 512;
  const int ny = 512;
  const int ns = 64;

  int *image = new int[nx * ny * 3];

  os << "P3" << std::endl;
  os << nx << " " << ny << std::endl;
  os << "255" << std::endl;

  HitableList world;
  Hitable *light = nullptr;

  CornellBox(world, &light);

  light = world.list[world.list.size() - 1].get();

  Vec3 lookfrom(0.0f, 4.5f, -4.0f), lookat(0.0f, 0.5f, 0.5f);

  float dist_to_focus = (lookfrom - lookat).length();
  float aperture = 0.0f;
  Camera cam(lookfrom, lookat, Vec3(0, 1, 0), 35, float(nx / ny), aperture,
             dist_to_focus);

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; s++) {
        double u = float(i + RAND()) / float(nx);
        double v = float(j + RAND()) / float(ny);
        //       double u = float(i) / float(nx);
        //        double v = float(j) / float(ny);

        Ray r = cam.GetRay(u, v);
        col += colourRecursive(r, world, 0);
      }

      col /= float(ns);
      col[0] = col[0] > 1 ? 1 : col[0];
      col[1] = col[1] > 1 ? 1 : col[1];
      col[2] = col[2] > 1 ? 1 : col[2];
      float exponent = 1.0f / 2.2f;
      col = Vec3(powf(col[0], exponent), powf(col[1], exponent),
                 powf(col[2], exponent));

      int y_pixel = ny - j - 1;
      int index = (i + y_pixel * nx) * 3;
      image[index] = int(255.99 * col[0]);
      image[index + 1] = int(255.99 * col[1]);
      image[index + 2] = int(255.99 * col[2]);
    }
  }

  for (int i = 0; i < nx * ny * 3; i += 3) {
    os << image[i] << " " << image[i + 1] << " " << image[i + 2] << "\n";
  }

  os.close();

  delete[] image;
  return 0;
}
