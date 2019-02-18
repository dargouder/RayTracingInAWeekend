#include "camera.h"
#include "dielectric.h"
#include "diffuselight.h"
#include "float.h"
#include "hitable_list.h"
#include "lambertian.h"
#include "metal.h"
#include "microfacet.h"
#include "quad.h"
#include "sphere.h"
#include "phong.h"

#include <fstream>
#include <iostream>

Vec3 colourRecursive(const Ray &ray, const Hitable &world, int depth) {
  HitRecord rec;
  if (world.hit(ray, 0.001, FLT_MAX, rec)) {
    Ray scattered;
    Vec3 attenuation;
    float u, v;
    Vec3 emitted = rec.mat_ptr->Le(u, v, rec.p);
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
    ;
  }
}

Vec3 colour(const Ray &ray, const Hitable &world, int depth) {
  Ray r = ray;
  Vec3 currentAttenuation = Vec3(1.0f, 1.0f, 1.0f);
  const int maxBounces = 10;
  for (int i = 0; i < maxBounces; ++i) {
    HitRecord rec;
    if (world.hit(r, 0.001, FLT_MAX, rec)) {
      Ray scattered;
      Vec3 attenuation;
      float u, v;
      Vec3 emitted = rec.mat_ptr->Le(u, v, rec.p);
      rec.normal.make_unit_vector();
      float pdf = 0.0f;
      if (rec.mat_ptr->fr(r, rec, attenuation, scattered, pdf)) {
        scattered.direction().make_unit_vector();
        rec.normal.make_unit_vector();
        float theta = Vec3::dot(scattered.direction(), rec.normal);
        currentAttenuation *= (attenuation / pdf) * cos(theta);
        currentAttenuation += emitted;
        r = scattered;
      } else {
        return emitted;  // Vec3(0, 0, 0);
      }

    } else {
      Vec3 unit_direction = Vec3::unit_vector(r.direction());
      float t = 0.5f * (unit_direction.y() + 1.0f);
      Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
      return currentAttenuation * col;
    }
  }
  return currentAttenuation;
}

Vec3 colourNormal(const Ray &ray, const Hitable &world, int depth) {
  Ray r = ray;
  Vec3 currentAttenuation = Vec3(1.0f, 1.0f, 1.0f);

  HitRecord rec;
  if (world.hit(r, 0.001, FLT_MAX, rec)) {
    // return Vec3(1.0, 0.0, 0.0);
    // return rec.normal;
    currentAttenuation =
        0.5f * Vec3(rec.normal.x() + 1.0f, rec.normal.y() + 1.0f,
                    rec.normal.z() + 1.0f);
    return currentAttenuation;
  } else {
    Vec3 unit_direction = Vec3::unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
    return currentAttenuation * col;
  }
}

void RandomScene(HitableList &list) {
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(0, -1000, 0), 1000,
      std::make_unique<Lambertian>(Vec3(0.5, 0.5, 0.5))));

  for (int a = -10; a < 0; a++) {
    for (int b = -10; b < 0; b++) {
      float choose_mat = RAND();
      Vec3 center(a + 0.9 * RAND(), 0.2, b + 0.9 * RAND());
      if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
        if (choose_mat < 0.8) {
          list.list.push_back(std::make_unique<Sphere>(
              center, 0.2,
              std::make_unique<Lambertian>(
                  Vec3(RAND() * RAND(), RAND() * RAND(), RAND() * RAND()))));
        } else if (choose_mat < 0.95) {
          list.list.push_back(std::make_unique<Sphere>(
              center, 0.2,
              std::make_unique<Metal>(
                  Vec3(0.5 * (1 + RAND()), 0.5 * (1 + RAND()),
                       0.5 * (1 + RAND())),
                  0.5 * RAND())));
        } else {
          list.list.push_back(std::make_unique<Sphere>(
              center, 0.2, std::make_unique<Dielectric>(1.5)));
        }
      }
    }
  }

  list.list.push_back(std::make_unique<Sphere>(
      Vec3(0, 1, 0), 1.0, std::make_unique<Dielectric>(1.5)));
  ;
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(4, 1, 0), 1.0, std::make_unique<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));
}
void SimpleScene(HitableList &list) {
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-10.0f, -1.2f, 10.0f), Vec3(10.0f, -1.2f, 10.0f),
      Vec3(10.0f, -1.2f, -10.0f), Vec3(-10.0f, -1.2f, -10.0f),
      std::make_unique<Lambertian>(Vec3(0.4f, 0.4f, 0.4f))));
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(0, 0, 0), 1.0,
		std::make_unique<Lambertian>(Vec3(0.8, 0.1, 0.1))));

  //list.list.push_back(std::make_unique<Sphere>(
  //    Vec3(0, 0, 0), 1.0,
		//std::make_unique<Phong>(Vec3(0.8f, 0.1f, 0.1f), Vec3(0.0f, 0.0f, 0.0f))));

/*    list.list.push_back(std::make_unique<Sphere>(
      Vec3(0, 0, 0), 1.0,
      std::make_unique<Phong>(Vec3(0.2f, 0.2f, 0.2f), Vec3(0.8f, 0.8f, 0.8f))))*/;

//  list.list.push_back(std::make_unique<Sphere>(
//      Vec3(0, 0, 0), 1.0,
//      std::make_unique<Dielectric>(1.5f)));
  //list.list.push_back(std::make_unique<Sphere>(
  //    Vec3(0, 0, 0), 1.0, std::make_unique<Metal>(Vec3(1.0, 1.0, 1.0), 0.0)));

  //  list.list.push_back(std::make_unique<Sphere>(
  //   Vec3(0, 0, 0), 1.0, std::make_unique<Microfacet>(Vec3(0.8, 0.1, 0.1),
  //   Vec3(1.0, 1.0f, 1.0f), 0.52f)));
}
void MISScene(HitableList &list) {
  // Lights
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(10.0f, 10.0f, 4.0f), 0.5,
      std::make_unique<DiffuseLight>(Vec3(800.0, 800.0, 800.0))));
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(-1.25f, 0.0f, 0.0f), 0.1,
      std::make_unique<DiffuseLight>(Vec3(100.0, 100.0, 100.0))));
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(-3.75f, 0.0f, 0.0f), 0.03333f,
      std::make_unique<DiffuseLight>(Vec3(900.0, 900.0, 900.0))));
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(1.25f, 0.0f, 0.0f), 0.3f,
      std::make_unique<DiffuseLight>(Vec3(10.0f, 10.0f, 10.0f))));
  list.list.push_back(std::make_unique<Sphere>(
      Vec3(3.75f, 0.0f, 0.0f), 0.9f,
      std::make_unique<DiffuseLight>(Vec3(1.0f, 1.0f, 1.0f))));

  // Floor
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-10.0f, -4.2f, 10.0f), Vec3(10.0f, -4.2f, 10.0f),
      Vec3(10.0f, -4.2f, -10.0f), Vec3(-10.0f, -4.2f, -10.0f),
      std::make_unique<Lambertian>(Vec3(0.4f, 0.4f, 0.4f))));

  list.list.push_back(std::make_unique<Quad>(
      Vec3(-10.0f, -10.0f, -2.0f), Vec3(10.0f, -10.0f, -2.0f),
      Vec3(10.0f, 10.0f, -2.0f), Vec3(-10.0f, 10.0f, -2.0f),
      std::make_unique<Lambertian>(Vec3(0.4f, 0.4f, 0.4f))));

  // Plate 1
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-4.0f, -2.70651f, 0.25609f), Vec3(4, -2.70651f, 0.25609f),
      Vec3(4.0f, -2.08375f, -0.526323f), Vec3(-4.0f, -2.08375f, -0.526323f),
      std::make_unique<Metal>(Vec3(1.0, 1.0, 1.0), 0.005f)));

  // Plate 2
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-4.0f, -3.28825f, 1.36972f), Vec3(4.0f, -3.28825f, 1.36972f),
      Vec3(4.0f, -2.83856f, 0.476536f), Vec3(-4.0f, -2.83856f, 0.476536f),
      std::make_unique<Metal>(Vec3(1.0, 1.0, 1.0), 0.02f)));

  // Plate 3
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-4.0f, -3.73096f, 2.70046f), Vec3(4.0f, -3.73096f, 2.70046f),
      Vec3(4.0f, -3.43378f, 1.74564f), Vec3(-4.0f, -3.43378f, 1.74564f),
      std::make_unique<Metal>(Vec3(1.0, 1.0, 1.0), 0.05f)));

  // Plate 4
  list.list.push_back(std::make_unique<Quad>(
      Vec3(-4.0f, -3.99615f, 4.0667f), Vec3(4.0f, -3.99615f, 4.0667f),
      Vec3(4.0f, -3.82069f, 3.08221f), Vec3(-4.0f, -3.82069f, 3.08221f),
      std::make_unique<Metal>(Vec3(1.0, 1.0, 1.0), 0.125f)));
}

int main() {
  std::ofstream os;
  os.open("mis.ppm", std::ios::binary);
  int nx = 512;
  int ny = 512;
  int ns = 32;

  os << "P3" << std::endl;
  os << nx << " " << ny << std::endl;
  os << "255" << std::endl;

  HitableList world;

  // MISScene(world);
  SimpleScene(world);
  // RandomScene(world);
  // Vec3 lookfrom(0.0f, 2.0f, 15.0f), lookat(0.0f, -2.0f, 2.5f);
  Vec3 lookfrom(0.0f, 0.0f, -10.0f), lookat(0.0f, 0.0f, 0.0f);
  // Vec3 lookfrom(13.0f, 2.0f, 3.0f), lookat(0.0f, 0.0f, 0.0f);
  float dist_to_focus = (lookfrom - lookat).length();
  float aperture = 0.0f;
  Camera cam(lookfrom, lookat, Vec3(0, 1, 0), 35, float(nx / ny), aperture,
             dist_to_focus);

  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; s++) {
        double u = float(i + RAND()) / float(nx);
        double v = float(j + RAND()) / float(ny);

        Ray r = cam.GetRay(u, v);
        col += colour(r, world, 0);
      }

      col /= float(ns);
      col[0] = col[0] > 1 ? 1 : col[0];
      col[1] = col[1] > 1 ? 1 : col[1];
      col[2] = col[2] > 1 ? 1 : col[2];
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

      int ir = int(255.99 * col[0]);  // > 256 ? 256 : int(255.99 * col[0]);
      int ig = int(255.99 * col[1]);  // > 256 ? 256 : int(255.99 * col[1]);
      int ib = int(255.99 * col[2]);  // > 256 ? 256 : int(255.99 * col[2]);
      os << ir << " " << ig << " " << ib << "\n";
    }
  }

  os.close();
  return 0;
}
