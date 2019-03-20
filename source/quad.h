#pragma once

#include <memory>
#include "hitable.h"

class Quad : public Hitable {
 public:
  Vec3 m_v00, m_v10, m_v11, m_v01;

  std::unique_ptr<Material> material;

  Quad() {}
  Quad(Vec3 v00, Vec3 v10, Vec3 v11, Vec3 v01, std::unique_ptr<Material> mat)
      : m_v00(v00),
        m_v10(v10),
        m_v11(v11),
        m_v01(v01),
        material(std::move(mat)) {}

  bool hit(const Ray& r, float tmin, float tmax, HitRecord& rec) const {
    const float epsilon = 1e-6;
    float u, v;
    // Reject rays using the barycentric coordinates of the intersection point
    // with respect to T
    Vec3 e01 = m_v10 - m_v00;
    Vec3 e03 = m_v01 - m_v00;
    Vec3 P = Vec3::cross(r.direction(), e03);
    float det = Vec3::dot(e01, P);
    if (std::abs(det) < epsilon) {
      return false;
    }

    Vec3 T = r.origin() - m_v00;
    float alpha = Vec3::dot(T, P) / det;
    if (alpha < 0.0f) {
      return false;
    }
    Vec3 Q = Vec3::cross(T, e01);

    float beta = Vec3::dot(r.direction(), Q) / det;
    if (beta < 0.0f) {
      return false;
    }

    // Reject rays using the barycentric coordinates of the interection point
    // with respect to T'
    if ((alpha + beta) > 1.0f) {
      Vec3 E23 = m_v01 - m_v11;
      Vec3 E21 = m_v10 - m_v11;
      Vec3 P1 = Vec3::cross(r.direction(), E21);
      float det1 = Vec3::dot(E23, P1);
      if (std::abs(det1) < epsilon) {
        return false;
      }
      Vec3 T1 = r.origin() - m_v11;
      float alpha1 = Vec3::dot(T1, P1) / det1;

      if (alpha1 < float(0.0)) return false;

      Vec3 Q1 = Vec3::cross(T1, E23);
      float beta1 = Vec3::dot(r.direction(), Q1) / det1;
      if (beta1 < 0.0f) {
        return false;
      }
    }

    // Compute the ray parameter of the intersection point
    float t = Vec3::dot(e03, Q) / det;
    if (t < 0.0f) {
      return false;
    }

    Vec3 E01 = m_v10 - m_v00;
    Vec3 E02 = m_v11 - m_v00;
    Vec3 E03 = m_v01 - m_v00;

    float alpha11, beta11;

    Vec3 N = Vec3::cross(E01, E03);
    if ((std::abs(N.x) >= std::abs(N.y)) &&
        (std::abs(N.x) >= std::abs(N.z))) {
      alpha11 = (E02.y * E03.z - E02.z * E03.y) / N.x;
      beta11 = (E01.y * E02.z - E01.z * E02.y) / N.x;
    } else if ((std::abs(N.y) >= std::abs(N.x)) &&
               (std::abs(N.y) >= std::abs(N.z))) {
      alpha11 = (E02.z * E03.x - E02.x * E03.z) / N.y;
      beta11 = (E01.z * E02.x - E01.x * E02.z) / N.y;
    } else {
      alpha11 = (E02.x * E03.y - E02.y * E03.x) / N.z;
      beta11 = (E01.x * E02.y - E01.y * E02.x) / N.z;
    }

    // Compute the bilinear coordinaes of the intersection point
    if (std::abs(alpha11 - 1.0f) < epsilon) {
      u = alpha;
      if (abs(beta11 - 1.0f) < epsilon) {
        v = beta;
      } else {
        v = beta / (u * (beta11 - 1.0f) + 1.0f);
      }
    } else if (std::abs(beta11 - 1.0f) < epsilon) {
      v = beta;
      u = alpha / ((v * (alpha11 - 1.0f)) + 1.0f);
    } else {
      float A = 1.0f - beta11;
      float B = (alpha * (beta11 - 1.0f)) - (beta * (alpha11 - 1.0f)) - 1.0f;
      float C = alpha;
      float determinant = B * B - 4.0f * A * C;
      float signB = B < 0.0f ? -1.0f : 1.0f;
      float Q = -0.5f * (B + signB * std::sqrt(determinant));
      u = Q / A;
      if ((u < 0.0f) || (u > 1.0f)) {
        u = C / Q;
      }

      v = beta / ((u * (beta11 - 1.0f)) + 1.0f);
    }

    Vec3 poi = (1.0f - u) * (1.0f - v) * m_v00 + u * (1.0f - v) * m_v10 +
               u * v * m_v11 + (1.0f - u) * v * m_v01;

    rec.t = (poi - r.origin()).length();

    if (rec.t > tmax) {
      return false;
    }

    rec.mat_ptr = material.get();
    rec.normal = Vec3::cross(m_v10 - m_v00, m_v01 - m_v00);
    rec.normal.make_unit_vector();
    rec.p = r.point_at_parameter(rec.t);  // + epsilon*rec.normal;
    // rec.p = poi;
    // rec.p = poi;// + rec.normal;
    // rec.p = poi + Vec3(1, 1, 1)*epsilon + rec.normal;
    return true;
  }

  Vec3 generateSampleOnSurface() const {
    float u = RAND();
    float v = RAND();

    // u = 0;
    // v = 0;

    Vec3 poi = (1.0f - u) * (1.0f - v) * m_v00 + u * (1.0f - v) * m_v10 +
               u * v * m_v11 + (1.0f - u) * v * m_v01;

    return poi;
  }

  float Pdf() const {
    Vec3 E01 = m_v10 - m_v00;
    Vec3 E03 = m_v01 - m_v00;

    float length_01 = E01.length();
    float length_03 = E03.length();

    return 1.0f / (length_01 * length_03);
  }
};
