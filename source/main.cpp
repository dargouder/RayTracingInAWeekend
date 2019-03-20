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

#include "sphere.h"

#include <cassert>
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
    if (depth < 10 &&
        rec.mat_ptr->sample_f(ray, rec, attenuation, scattered, pdf)) {
      Vec3 f = attenuation / pdf;
      //assert(pdf > 0);
      float cos_theta = Vec3::dot(rec.normal, scattered.direction());
      return emitted +
             f * cos_theta * colourRecursive(scattered, world, depth + 1);
    } else {
      return emitted;
    }
  } else {
    Vec3 unit_direction = Vec3::unit_vector(ray.direction());
    float t = 0.5f * (unit_direction.y + 1.0f);
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
     list.list.push_back(std::make_unique<Sphere>(
         Vec3(0.0f, 0.5f, 0.5f), 0.5f,
         std::make_unique<Lambertian>(Vec3(0.0f, 0.9f, 0.0f))));
  //list.list.push_back(std::make_unique<Sphere>(
  //    Vec3(0.0f, 0.5f, 0.5f), 0.5f,
  //    std::make_unique<Phong>(Vec3(0.4f, 0.4f, 0.4f), Vec3(0.05f, 0.05f, 0.05f))));
}

void testPdf() {
  float power = 10;
  Vec3 normal(0, 1, 0);

  ONB onb;
  onb.branchlessONB(normal);

  Vec3 wo(-1.0f, -1.0f, 0.0f);
  wo.make_unit_vector();
  float pdf = 0.0f;
  Vec3 perfectReflection = reflect(wo, normal);
  Vec3 m_kd(0.4, 0.4, 0.4);
  Vec3 m_ks(0.05, 0.05, 05);

  Vec3 scatteredDir;

  float diffuseRatio =
      m_kd.getLuminance() / (m_kd.getLuminance() + m_ks.getLuminance());
  for (int i = 0; i < 100000; ++i) {
    float reflectionDecision = RAND();

    if (reflectionDecision < diffuseRatio || diffuseRatio == 1.0f) {
      // calculating only diffuse component
      float u = RAND();
      float v = RAND();

      scatteredDir = onb.local(CosineSampleHemispherePhong(u, v));
      scatteredDir.make_unit_vector();

      pdf = Vec3::dot(scatteredDir, normal) * INV_PI;
    } else {
      // calculating only specular component
      float u = RAND();
      float v = RAND();

      float sin_alpha = std::sqrt(1 - std::powf(u, 2.0f / (power + 1.0f)));
      float cos_alpha = std::powf(u, 1.0f / (power + 1));
      float phi = 2 * M_PI * v;
      float cos_phi = std::cos(phi);
      float sin_phi = std::sin(phi);

      scatteredDir =
          onb.local(Vec3(sin_alpha * cos_phi, sin_alpha * sin_phi, cos_alpha));
      scatteredDir.make_unit_vector();

      float cos_theta = Vec3::dot(scatteredDir, perfectReflection);

      cos_theta = clamp(cos_theta, M_PI / 2.0f, -M_PI / 2.0f);

      pdf += (power + 2.0f) * INV_2PI * powf(cos_theta, power + 1.0f);
    }
  }

  int i = 0;
}

int main() {
  //testPdf();
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

  Vec3 lookfrom(0.0f, 0.5f, -4.0f), lookat(0.0f, 0.5f, 0.5f);

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
        Ray r = cam.GetRay(u, v);
        col += colourRecursive(r, world, 0);
      }

      col /= float(ns);
      col.x = col.x > 1 ? 1 : col.x;
      col.y = col.x > 1 ? 1 : col.y;
      col.z = col.x > 1 ? 1 : col.z;
      float exponent = 1.0f / 2.2f;
      col = Vec3(powf(col.x, exponent), powf(col.y, exponent),
                 powf(col.z, exponent));

      int y_pixel = ny - j - 1;
      int index = (i + y_pixel * nx) * 3;
      image[index] = int(255.99 * col.x);
      image[index + 1] = int(255.99 * col.y);
      image[index + 2] = int(255.99 * col.z);
    }
  }

  for (int i = 0; i < nx * ny * 3; i += 3) {
    os << image[i] << " " << image[i + 1] << " " << image[i + 2] << "\n";
  }

  os.close();

  delete[] image;
  return 0;
}
