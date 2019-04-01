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
#include "reflection.h"
#include "sphere.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <cassert>
#include <fstream>
#include <iostream>

ONB shadingONB;
Vec3 colourNormal(const Ray &ray, const Hitable &world, int depth) {
  HitRecord rec;
  if (world.hit(ray, 0.001, FLT_MAX, rec)) {

    rec.normal.make_unit_vector();

    return (rec.normal + Vec3(1.0f, 1.0f, 1.0f)) * 0.5;

  } else {
    Vec3 unit_direction = Vec3::unit_vector(ray.direction());
    float t = 0.5f * (unit_direction.y + 1.0f);
    Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
    return col;
  }
}

Vec3 colourRecursive(const Ray &ray, const Hitable &world, int depth) {
  HitRecord rec;
  if (world.hit(ray, 0.001, FLT_MAX, rec)) {
    Ray scattered;
    Vec3 attenuation;
    float u = 0.0f, v = 0.0f;
    rec.normal.make_unit_vector();

    if (depth < 10) {
      Vec3 wo = shadingONB.local(-ray.direction());
      wo.make_unit_vector();
      Vec3 wi;
      float pdf;

      Vec3 f = rec.bsdf->bxdf->Sample_f(wo, wi, pdf);
      // wi.make_unit_vector();
      f = f * AbsCosTheta(wi) / pdf;
      // if (!(pdf > 0))
      //{
      //  return Vec3(0.0f, 0.0f, 0.0f);
      //}
      // assert(pdf > 0);
      ONB shapeONB;
      shapeONB.branchlessONB(rec.normal);
      Vec3 origin = rec.p + rec.normal * 0.001f;
      scattered = Ray(origin, shapeONB.local(wi));

      float pdfAfter =
          std::abs(Vec3::dot(scattered.direction(), rec.normal)) / PI;

      // assert(abs(pdf - pdfAfter) < 0.000001f);
      // scattered.B.make_unit_vector();
      return f * colourRecursive(scattered, world, depth + 1);
    } else {
      return Vec3(0, 0, 0);
    }
  } else {
    Vec3 unit_direction = Vec3::unit_vector(ray.direction());
    float t = 0.5f * (unit_direction.y + 1.0f);
    Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
    return col;
  }
}

// Vec3 colourRecursive(const Ray &ray, const Hitable &world, int depth) {
//  HitRecord rec;
//  if (world.hit(ray, 0.001, FLT_MAX, rec)) {
//    Ray scattered;
//    Vec3 attenuation;
//    float u = 0.0f, v = 0.0f;
//    rec.normal.make_unit_vector();
//    Vec3 emitted = rec.mat_ptr->Le(0, 0, rec.p);
//    float pdf;
//    if (depth < 10 &&
//        rec.mat_ptr->sample_f(ray, rec, attenuation, scattered, pdf)) {
//      Vec3 f = attenuation / pdf;
//      // assert(pdf > 0);
//      float cos_theta = Vec3::dot(rec.normal, scattered.direction());
//      return emitted +
//             f * cos_theta * colourRecursive(scattered, world, depth + 1);
//    } else {
//      return emitted;
//    }
//  } else {
//    Vec3 unit_direction = Vec3::unit_vector(ray.direction());
//    float t = 0.5f * (unit_direction.y + 1.0f);
//    Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
//    return col;
//  }
//}

// void CornellBox(HitableList &list, Hitable **light) {
//  // floor
//  list.list.push_back(std::make_unique<Quad>(
//      Vec3(-10.0f, 0.0f, -10.0f), Vec3(-10.0f, 0.0f, 10.0f),
//      Vec3(10.0f, 0.0f, 10.0f), Vec3(10.0f, 0.0f, -10.0f),
//      std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
//  //list.list.push_back(std::make_unique<Sphere>(
//  //    Vec3(0.0f, 0.5f, 0.5f), 0.5f,
//  //    std::make_unique<LambertianReflection>(Vec3(0.0f, 0.9f, 0.0f))));
//
//  //list.list.push_back(std::make_unique<Quad>(
//  //    Vec3(-0.5f, 1.0f, -0.5f), Vec3(-0.5f, 1.0f, 0.5f),
//  //    Vec3(0.5f, 1.0f, 0.5f), Vec3(0.5f, 1.0f, -0.5f),
//  //    std::make_unique<LambertianReflection>(Vec3(1.0f, 1.0f, 0.0f))));
//
//  list.list.push_back(std::make_unique<Quad>(
//      Vec3(-0.5f, 1.0f, -0.5f), Vec3(-0.5f, 1.0f, 0.5f),
//      Vec3(0.5f, 1.0f, 0.5f), Vec3(0.5f, 1.0f, -0.5f),
//      std::make_unique<AreaLight>(Vec3(1.0f, 1.0f, 1.0f))));
//
//  // list.list.push_back(std::make_unique<Sphere>(
//  //    Vec3(0.0f, 0.5f, 0.5f), 0.5f,
//  //    std::make_unique<Phong>(Vec3(0.4f, 0.4f, 0.4f), Vec3(0.05f, 0.05f,
//  //    0.05f))));
//}

void CornellBox(HitableList &list, Hitable *light) {
  //// floor
  // list.list.push_back(std::make_unique<Quad>(
  //    Vec3(552.8f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f),
  //    Vec3(0.0f, 0.0f, 559.2f), Vec3(549.6f, 0.0f, 559.2f),
  //    std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));

  //// ceiling
  // list.list.push_back(std::make_unique<Quad>(
  //    Vec3(556.0f, 548.8f, 0.0f), Vec3(556.0f, 548.8f, 559.2f),
  //    Vec3(0.0f, 548.8f, 559.2f), Vec3(0.0f, 548.8f, 0.0f),
  //    std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));

  //// back wall
  // list.list.push_back(std::make_unique<Quad>(
  //    Vec3(549.6f, 0.0f, 559.2f), Vec3(0.0f, 0.0f, 559.2f),
  //    Vec3(0.0f, 548.8f, 559.2f), Vec3(549.6f, 548.8f, 559.2f),
  //    std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));

  // Sphere
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(0, 0, -2), 1,
      std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));

  //// right wall
  // list.list.push_back(std::make_unique<Quad>(
  //    Vec3(0.0f, 0.0f, 559.2f), Vec3(0.0f, 0.0f, 0.0f),
  //    Vec3(0.0f, 548.8f, 0.0f), Vec3(0.0f, 548.8f, 559.2f),
  //    std::make_unique<LambertianReflection>(Vec3(0.14f, 0.45f, 0.091f))));

  //// left wall
  // list.list.push_back(std::make_unique<Quad>(
  //    Vec3(549.6f, 0.0f, 0.0f), Vec3(549.6f, 0.0f, 559.2f),
  //    Vec3(549.6f, 548.8f, 559.2f), Vec3(549.6f, 548.8f, 0.0f),
  //    std::make_unique<LambertianReflection>(Vec3(0.63f, 0.065f, 0.05f))));

  // short block
  // list.list.push_back(std::make_unique<Quad>(
  //  Vec3(130.0, 165.0, 65.0), Vec3(82.0, 165.0, 225.0),
  //  Vec3(240.0, 165.0, 272.0), Vec3(290.0, 165.0, 114.0),
  //  std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
  // list.list.push_back(std::make_unique<Quad>(
  //  Vec3(290.0, 0.0, 114.0), Vec3(290.0, 165.0, 114.0),
  //  Vec3(240.0, 165.0, 272.0), Vec3(240.0, 0.0, 272.0),
  //  std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
  // list.list.push_back(std::make_unique<Quad>(
  //  Vec3(130.0, 0.0, 65.0), Vec3(130.0, 165.0, 65.0),
  //  Vec3(290.0, 165.0, 114.0), Vec3(290.0, 0.0, 114.0),
  //  std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
  // list.list.push_back(std::make_unique<Quad>(
  //  Vec3(82.0, 0.0, 225.0), Vec3(82.0, 165.0, 225.0),
  //  Vec3(130.0, 165.0, 65.0), Vec3(130.0, 0.0, 65.0),
  //  std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
  // list.list.push_back(std::make_unique<Quad>(
  //  Vec3(240.0, 0.0, 272.0), Vec3(240.0, 165.0, 272.0),
  //  Vec3(82.0, 165.0, 225.0), Vec3(82.0, 0.0, 225.0),
  //  std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
  // long block
  /* list.list.push_back(std::make_unique<Quad>(
     Vec3(423.0, 330.0, 247.0), Vec3(265.0, 330.0, 296.0),
     Vec3(314.0, 330.0, 456.0), Vec3(472.0, 330.0, 406.0),
     std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
   list.list.push_back(std::make_unique<Quad>(
     Vec3(423.0, 0.0, 247.0), Vec3(423.0, 330.0, 247.0),
     Vec3(472.0, 330.0, 406.0), Vec3(472.0, 0.0, 406.0),
     std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
   list.list.push_back(std::make_unique<Quad>(
     Vec3(472.0, 0.0, 406.0), Vec3(472.0, 330.0, 406.0),
     Vec3(314.0, 330.0, 456.0), Vec3(314.0, 0.0, 456.0),
     std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
   list.list.push_back(std::make_unique<Quad>(
     Vec3(314.0, 0.0, 456.0), Vec3(314.0, 330.0, 456.0),
     Vec3(265.0, 330.0, 296.0), Vec3(265.0, 0.0, 296.0),
     std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));
   ;
   list.list.push_back(std::make_unique<Quad>(
     Vec3(265.0, 0.0, 296.0), Vec3(265.0, 330.0, 296.0),
     Vec3(423.0, 330.0, 247.0), Vec3(423.0, 0.0, 247.0),
     std::make_unique<LambertianReflection>(Vec3(0.725f, 0.71f, 0.68f))));*/

  // light
  // list.list.push_back(std::make_unique<Quad>(
  // Vec3(343.0f, 548.8f, 227.0f), Vec3(343.0f, 548.8f, 332.0f),
  // Vec3(213.0f, 548.8f, 332.2f), Vec3(213.0f, 548.8f, 227.0f),
  // std::make_unique<AreaLight>(Vec3(17.0f, 12.0f, 4.0f))));

  // light = list.list[list.list.size() - 1].get();
}

//#define SINGLE_RAY
int main() {
  shadingONB.pbrtONB(Vec3(0, 0, 1));

  int nx = 512;
  int ny = 512;
  int ns = 64;
  std::vector<char> img;

#ifdef SINGLE_RAY
  nx = 1;
  ny = 1;
  ns = 1;
#endif

  int *normal_image = new int[nx * ny * 3];
  int *image = new int[nx * ny * 3];

  HitableList world;
  Hitable *light = nullptr;

  CornellBox(world, light);

  light = world.list[world.list.size() - 1].get();

  Vec3 lookfrom(0.0f, 0.0f, -6.0f), lookat(0.0f, 0.0f, -1.0f);
  //Vec3 lookfrom(278.0f, 292.0f, -800.0f), lookat(278.0f, 292.0f, -799.0f);
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
      Vec3 normal(0.0f, 0.0f, 0.0f);
      {
        // Trace normal
        float u = float(i + 0.5f) / float(nx);
        float v = float(j + 0.5f) / float(ny);

        Ray r = cam.GetRay(u, v);

        normal = colourNormal(r, world, 0);

        int y_pixel = ny - j - 1;
        int index = (i + y_pixel * nx) * 3;
        normal_image[index] = int(255.99 * normal.x);
        normal_image[index + 1] = int(255.99 * normal.y);
        normal_image[index + 2] = int(255.99 * normal.z);
      }
      for (int s = 0; s < ns; s++) {
#ifdef SINGLE_RAY
        double u = float(i) / float(nx);
        double v = float(j) / float(ny);

#else
        double u = float(i + RAND()) / float(nx);
        double v = float(j + RAND()) / float(ny);
#endif
        Ray r = cam.GetRay(u, v);
        // r.B = Vec3(0,0,1);
        col += colourRecursive(r, world, 0);
      }

      col /= float(ns);
      // col.x = col.x > 1 ? 1 : col.x;
      // col.y = col.x > 1 ? 1 : col.y;
      // col.z = col.x > 1 ? 1 : col.z;
      col.x = clamp(col.x, 1.0f, 0.0f);
      col.y = clamp(col.y, 1.0f, 0.0f);
      col.z = clamp(col.z, 1.0f, 0.0f);

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

  for (int i = 0; i < nx * ny * 3; ++i) {
    img.push_back((char)image[i]);
  }

  stbi_write_bmp("test.bmp", nx, ny, 3, img.data());

  img.clear();
  for (int i = 0; i < nx * ny * 3; ++i) {
    img.push_back((char)normal_image[i]);
  }

  stbi_write_bmp("normals.bmp", nx, ny, 3, img.data());

  delete[] normal_image;
  delete[] image;
  return 0;
}
