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

Vec3 directLighting(int shadowSamples, const Hitable &world,
                    const Hitable *light, Vec3 attenuation, float pdf,
                    HitRecord rec) {
  Vec3 shadowAccumulation(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < shadowSamples; ++i) {
    // Generate shadow ray
    Vec3 lightSample = light->generateSampleOnSurface();
    Ray shadowRay(rec.p, Vec3::unit_vector(lightSample - rec.p));
    HitRecord shadowRec;

    if (world.hit(shadowRay, 0.001, FLT_MAX, shadowRec)) {
      // if the closest hit was the light
      if (shadowRec.mat_ptr == ((Quad *)light)->material.get()) {
        float lightPdf = light->Pdf();
        Vec3 Le = shadowRec.mat_ptr->Le(0.0f, 0.0f, rec.p);
        float distanceSquared = (rec.p - shadowRec.p).squared_length();
        float cos_theta_x =
            std::max(Vec3::dot(rec.normal, shadowRay.direction()), 0.0f);
        float cos_theta_y =
            std::max(Vec3::dot(shadowRec.normal, -shadowRay.direction()), 0.0f);
        float G = (cos_theta_x * cos_theta_y) / distanceSquared;
        shadowAccumulation += (Le * attenuation * pdf * G) / lightPdf;
      } else {
        shadowAccumulation += Vec3(0.0, 0.0, 0.0f);
      }
    }
  }

  return shadowAccumulation / float(shadowSamples);
}

Vec3 indirectLighting(const Ray &ray, const Hitable &world,
                      const Hitable *light, int depth) {
  return Vec3(0.0f, 0.0f, 0.0f);
}

Vec3 render(const Ray &ray, const Hitable &world, const Hitable *light,
            int depth) {
  HitRecord rec;

  // Shoot primary ray
  if (world.hit(ray, 0.001, FLT_MAX, rec)) {
    Ray scattered;
    Vec3 attenuation;
    float u = 0.0f, v = 0.0f;
    Vec3 emitted = rec.mat_ptr->Le(u, v, rec.p);
    float pdf;
    int shadowSamples = 4;

	
    if (depth < 10 && rec.mat_ptr->fr(ray, rec, attenuation, scattered, pdf)) {
    // We hit a normal surface so calculate the direct lighting.
		Vec3 direct =
          directLighting(shadowSamples, world, light, attenuation, pdf, rec);

      return direct;
    } else {
		// When we hit a light, just return what it emits
      return emitted;
    }
  } else {
    return Vec3(0, 0, 0);
  }
}

void CornellBox(HitableList &list, Hitable *light) {
  // floor
  list.list.push_back(std::make_unique<Quad>(
      Vec3(552.8f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f),
      Vec3(0.0f, 0.0f, 559.2f), Vec3(549.6f, 0.0f, 559.2f),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));

  // ceiling
  list.list.push_back(std::make_unique<Quad>(
      Vec3(556.0f, 548.8f, 0.0f), Vec3(556.0f, 548.8f, 559.2f),
      Vec3(0.0f, 548.8f, 559.2f), Vec3(0.6f, 548.8f, 0.0f),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));

  // back wall
  list.list.push_back(std::make_unique<Quad>(
      Vec3(549.6f, 0.0f, 559.2f), Vec3(0.0f, 0.0f, 559.2f),
      Vec3(0.0f, 548.8f, 559.2f), Vec3(556.6f, 548.8f, 559.2f),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));

  // right wall
  list.list.push_back(std::make_unique<Quad>(
      Vec3(0.0f, 0.0f, 559.2f), Vec3(0.0f, 0.0f, 0.0f),
      Vec3(0.0f, 548.8f, 0.0f), Vec3(0.0f, 548.8f, 559.2f),
      std::make_unique<Lambertian>(Vec3(0.14f, 0.45f, 0.091f))));

  // left wall
  list.list.push_back(std::make_unique<Quad>(
      Vec3(549.6f, 0.0f, 0.0f), Vec3(549.6f, 0.0f, 559.2f),
      Vec3(549.6f, 548.8f, 559.2f), Vec3(549.6f, 548.8f, 0.0f),
      std::make_unique<Lambertian>(Vec3(0.63f, 0.065f, 0.05f))));

  // short block
  list.list.push_back(std::make_unique<Quad>(
      Vec3(130.0, 165.0, 65.0), Vec3(82.0, 165.0, 225.0),
      Vec3(240.0, 165.0, 272.0), Vec3(290.0, 165.0, 114.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(290.0, 0.0, 114.0), Vec3(290.0, 165.0, 114.0),
      Vec3(240.0, 165.0, 272.0), Vec3(240.0, 0.0, 272.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(130.0, 0.0, 65.0), Vec3(130.0, 165.0, 65.0),
      Vec3(290.0, 165.0, 114.0), Vec3(290.0, 0.0, 114.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(82.0, 0.0, 225.0), Vec3(82.0, 165.0, 225.0),
      Vec3(130.0, 165.0, 65.0), Vec3(130.0, 0.0, 65.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(240.0, 0.0, 272.0), Vec3(240.0, 165.0, 272.0),
      Vec3(82.0, 165.0, 225.0), Vec3(82.0, 0.0, 225.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  // long block
  list.list.push_back(std::make_unique<Quad>(
      Vec3(423.0, 330.0, 247.0), Vec3(265.0, 330.0, 296.0),
      Vec3(314.0, 330.0, 456.0), Vec3(472.0, 330.0, 406.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(423.0, 0.0, 247.0), Vec3(423.0, 330.0, 247.0),
      Vec3(472.0, 330.0, 406.0), Vec3(472.0, 0.0, 406.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(472.0, 0.0, 406.0), Vec3(472.0, 330.0, 406.0),
      Vec3(314.0, 330.0, 456.0), Vec3(314.0, 0.0, 456.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  list.list.push_back(std::make_unique<Quad>(
      Vec3(314.0, 0.0, 456.0), Vec3(314.0, 330.0, 456.0),
      Vec3(265.0, 330.0, 296.0), Vec3(265.0, 0.0, 296.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));
  ;
  list.list.push_back(std::make_unique<Quad>(
      Vec3(265.0, 0.0, 296.0), Vec3(265.0, 330.0, 296.0),
      Vec3(423.0, 330.0, 247.0), Vec3(423.0, 0.0, 247.0),
      std::make_unique<Lambertian>(Vec3(0.725f, 0.71f, 0.68f))));

  // light
  list.list.push_back(std::make_unique<Quad>(
      Vec3(343.0f, 548.8f, 227.0f), Vec3(343.0f, 548.8f, 332.0f),
      Vec3(213.0f, 548.8f, 332.2f), Vec3(213.0f, 548.8f, 227.0f),
      std::make_unique<DiffuseLight>(Vec3(17.0f, 12.0f, 4.0f))));

  // light = list.list[list.list.size() - 1].get();
}

int main() {
  std::ofstream os;
  os.open("directLighting.ppm", std::ios::binary);
  const int nx = 1024;
  const int ny = 1024;
  const int ns = 512;

  int *image = new int[nx * ny * 3];

  os << "P3" << std::endl;
  os << nx << " " << ny << std::endl;
  os << "255" << std::endl;

  HitableList world;
  Hitable *light = nullptr;

  CornellBox(world, light);

  light = world.list[world.list.size() - 1].get();

  Vec3 lookfrom(278.0f, 292.0f, -800.0f), lookat(278.0f, 292.0f, -799.0f);

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
        col += render(r, world, light, 0);
      }

      col /= float(ns);
      col[0] = col[0] > 1 ? 1 : col[0];
      col[1] = col[1] > 1 ? 1 : col[1];
      col[2] = col[2] > 1 ? 1 : col[2];
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

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
