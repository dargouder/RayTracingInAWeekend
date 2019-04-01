#pragma once

#include "vec3.h"

/* Spherical coordinate functions */
inline static Vec3 SphericalDirection(float sinTheta, float cosTheta,
                                      float phi) {
  return Vec3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

inline static Vec3 SphericalDirection(float sinTheta, float cosTheta, float phi,
                                      const Vec3 &x, const Vec3 &y,
                                      const Vec3 &z) {
  return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
         cosTheta * z;
}

inline static float SphericalTheta(const Vec3 &v) {
  return std::acos(clamp(v.z, 1, -1));
}

inline static float SphericalPhi(const Vec3 &v) {
  float p = std::atan2(v.y, v.x);
  return (p < 0) ? p + 2 * PI : PI;
}

/*

Reflection utility functions.
These functions assume that the coordinate system of the vectors passed have z
up. The 2 tangent vectors and the normal vector at the point being shaded are
aligned with the x, y and z axes.

All directions passed to and returned by the BxDFs will respect this coordinate
system.

theta is the angle measured from the given direction to the z axis (normal). phi
os the angle formed with the x axis after projection of the direction onto the
xy plane.

*/

// cos theta = ( n . w ) = ( (0,0,1) . w) = w's z omponent
inline float CosTheta(const Vec3 &w) { return w.z; }

inline float Cos2Theta(const Vec3 &w) { return w.z * w.z; }

inline float AbsCosTheta(const Vec3 &w) { return fabsf(w.z); }

// sin theta^2 = 1 - cos theta^2
// the std::max avoids sqrt of a -ve number
inline float Sin2Theta(const Vec3 &w) {
  return std::max(static_cast<float>(0.0f),
                  static_cast<float>(1 - Cos2Theta(w)));
}

inline float SinTheta(const Vec3 &w) { return std::sqrt(Sin2Theta(w)); }

// tan theta = (sin theta) / (cos theta)
inline float TanTheta(const Vec3 &w) { return SinTheta(w) / CosTheta(w); }

inline float Tan2Theta(const Vec3 &w) { return Sin2Theta(w) / Cos2Theta(w); }

// cos phi = x / sin theta
inline float CosPhi(const Vec3 &w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 1 : clamp(w.x / sinTheta, -1, 1);
}

// sin phi = y / sin theta
inline float SinPhi(const Vec3 &w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 0 : clamp(w.y / sinTheta, -1, 1);
}

inline float Cos2Phi(const Vec3 &w) { return CosPhi(w) * CosPhi(w); }

inline float Sin2Phi(const Vec3 &w) { return SinPhi(w) * SinPhi(w); }

// used to check if the wo and wi lie on the same hemisphere.
inline bool SameHemisphere(const Vec3 &w, const Vec3 &wp) {
  return w.z * wp.z > 0;
}

class BxDF {
 public:
  // BxDF interface

  // returns the value of the distribution function for the given pair of
  // directions.
  virtual Vec3 f(const Vec3 &wo, const Vec3 &wi) const = 0;

  // Computes the direction of the incident light wi, the pdf for that direction
  // and the value of the distribution function for the pair of directions.
  virtual Vec3 Sample_f(const Vec3 &wo, Vec3 &wi, float &pdf) const;

  // Computes the pdf for the given wo and wi.
  virtual float Pdf(const Vec3 &wo, const Vec3 &wi) const;

  virtual bool isLight() { return false; }

   ~BxDF() = default;

};

class LambertianReflection : public BxDF {
 public:
  // Public methods
  explicit LambertianReflection(Vec3 kd) : m_kd(kd) {}

  Vec3 f(const Vec3 &wo, const Vec3 &wi) const { return m_kd / PI; }
  ~LambertianReflection() = default;
 private:
  // private data
  const Vec3 m_kd;


};

class AreaLight : public BxDF {
 public:
  AreaLight(Vec3 a) : emit(a) {}

  Vec3 f(const Vec3 &wo, const Vec3 &wi) const override { return emit; }

  bool isLight() override { return true; }

 private:
  Vec3 emit;
};

class BSDF {
 public:
  BSDF(BxDF *pbxdf) : bxdf(pbxdf) {}


  BxDF *bxdf;

  // shading normal, shading tangent, shading bitangent
  Vec3 ns, ss, ts;
};
