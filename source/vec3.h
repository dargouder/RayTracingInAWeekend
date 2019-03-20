#pragma once

#include <cmath>
#include <cstdlib>
#include <random>

#ifdef _WIN32
const float M_PI = 3.1415926535;
#endif

static const float INV_PI = 1.0f / M_PI;
static const float INV_2PI = 1.0f / (2.0f * M_PI);
static const float PiOver4 = M_PI / 4.0f;
static const float PiOver2 = M_PI / 2.0f;

static std::random_device
    rd;  // Will be used to obtain a seed for the random number engine

inline double RAND() {
  static std::mt19937 gen(
      rd());  // Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> dis(0.0, 1.0);
  return dis(gen);
}

inline float schlick(float cosine, float ref_idx) {
  float r0 = (1 - ref_idx) / (1 + ref_idx);
  r0 = r0 * r0;
  float cos_res = 1 - cosine;
  return r0 + (1 - r0) * cos_res * cos_res * cos_res * cos_res * cos_res;
}

inline float clamp(float power, float max, float min) {
  return std::min(max, std::max(min, power));
}

class Vec3 {
 public:
  float x, y, z;

  Vec3() {}
  Vec3(float e0, float e1, float e2) {
    x = e0;
    y = e1;
    z = e2;
  }

  inline float r() const { return x; }
  inline float g() const { return y; }
  inline float b() const { return z; }

  inline const Vec3& operator+() const { return *this; }

  inline Vec3 operator-() const { return Vec3(-x, -y, -z); }

  inline Vec3 operator+(const Vec3& v2) const {
    return Vec3(x + v2.x, y + v2.y, z + v2.z);
  }

  inline friend Vec3 operator-(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
  }

  inline friend Vec3 operator*(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
  }

  inline friend Vec3 operator/(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
  }

  inline friend Vec3 operator/(const Vec3& v1, float t) {
    return Vec3(v1.x / t, v1.y / t, v1.z / t);
  }

  inline friend Vec3 operator*(const Vec3& v1, float t) {
    return Vec3(v1.x * t, v1.y * t, v1.z * t);
  }

  inline friend Vec3 operator*(float t, const Vec3& v1) {
    return Vec3(v1.x * t, v1.y * t, v1.z * t);
  }

  inline Vec3& operator=(const Vec3& v1) {
    x = v1.x;
    y = v1.y;
    z = v1.z;

    return *this;
  }

  inline Vec3& operator+=(const Vec3& v2) {
    *this = *this + v2;
    return *this;
  }

  inline Vec3& operator-=(const Vec3& v2) {
    *this = *this - v2;
    return *this;
  }
  inline Vec3& operator*=(const Vec3& v2) {
    *this = *this * v2;
    return *this;
  }
  inline Vec3& operator/=(const Vec3& v2) {
    *this = *this / v2;
    return *this;
  }

  inline Vec3& operator*=(const float t) {
    *this = *this * t;
    return *this;
  }
  inline Vec3& operator/=(const float t) {
    *this = *this / t;
    return *this;
  }

  inline float length() const { return sqrt(x * x + y * y + z * z); }

  inline float squared_length() const { return x * x + y * y + z * z; }

  inline static float dot(const Vec3& v1, const Vec3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  inline static float satDot(const Vec3& v1, const Vec3& v2) {
    return clamp(dot(v1, v2), 0.0f, 1.0f);
  }

  inline static Vec3 cross(const Vec3& v1, const Vec3& v2) {
    return Vec3((v1.y * v2.z) - (v1.z * v2.y), -((v1.x * v2.z) - (v1.z * v2.x)),
                (v1.x * v2.y) - (v1.y * v2.x)

    );
  }

  inline void make_unit_vector() { *this = *this / length(); }

  inline static Vec3 unit_vector(Vec3 v) { return v / v.length(); }

  inline float getLuminance() const {
    return x * 0.212671f + y * 0.715160f + z * 0.072169f;
  }
};

inline Vec3 RandomInUnitDisk() {
  Vec3 p;
  do {
    p = 2.0 * Vec3(RAND(), RAND(), 0) - Vec3(1, 1, 0);
  } while (Vec3::dot(p, p) >= 1.0);

  return p;
}

static Vec3 reflect(const Vec3& wo, const Vec3& n) {
  return wo - 2 * Vec3::dot(wo, n) * n;
}

static Vec3 UniformSampleHemisphere(float u1, float u2) {
  const float r = std::sqrt(1.0f - u1 * u1);
  const float phi = 2 * M_PI * u2;

  return Vec3(cos(phi) * r, sin(phi) * r, u1);
}

static Vec3 CosineSampleHemisphere(float u, float v) {
  float z = sqrt(1 - v);
  float phi = 2 * M_PI * u;

  float x = cos(phi) * 2 * sqrt(v);
  float y = sin(phi) * 2 * sqrt(v);

  return Vec3(x, y, z);
}

static Vec3 CosineSampleHemispherePhong(float u, float v) {
  const float alpha = sqrtf(1.0f - u);
  const float beta = 2 * M_PI * v;

  const float x = alpha * cos(beta);
  const float y = alpha * sin(beta);
  const float z = sqrt(u);

  return Vec3(x, y, z);
}

static Vec3 CosineSampleHemisphereDriscoll(float u, float v) {
  const float r = sqrt(u);
  const float theta = 2 * M_PI * v;

  const float x = r * cos(theta);
  const float y = r * sin(theta);
  const float z = sqrt(1 - u);

  return Vec3(x, y, z);
}

static Vec3 ConcentricSampleDisk(float u, float v) {
  // Map uniform random numbers to [-1, 1]^2
  float uOffset = 2.0f * u - 1.0f;
  float vOffset = 2.0f * v - 1.0f;

  // Handle degeneracy at the origin
  if (uOffset == 0 && vOffset == 0) {
    return Vec3(0, 0, 0);
  }

  // Apply concentric mapping to point
  float theta, r;
  if (std::abs(uOffset) > std::abs(vOffset)) {
    r = uOffset;
    theta = PiOver4 * (vOffset / uOffset);
  } else {
    r = vOffset;
    theta = PiOver2 - PiOver4 * (uOffset / vOffset);
  }

  return Vec3(r * std::cos(theta), r *  std::sin(theta), 0.0f);
}

static Vec3 CosineSampleHemispherePBRT(float u, float v) {
  Vec3 d = ConcentricSampleDisk(u, v);
  float z = std::sqrt(std::max(
      (float)0, 1 - d.x * d.x - d.y * d.y));
  return Vec3(d.x, d.y, z);
}

static Vec3 UniformSampleSphere(float u, float v) {
  float theta = 2 * M_PI * u;
  float phi = acos(1 - 2 * v);
  float x = sin(phi) * cos(theta);
  float y = sin(phi) * sin(theta);
  float z = cos(phi);
  return Vec3(x, y, z);
}

static Vec3 RandomInUnitSphere() {
  Vec3 p;
  do {
    p = 2.0 * Vec3(RAND(), RAND(), RAND()) - Vec3(1, 1, 1);
  } while (Vec3::dot(p, p) >= 1.0);
  return p;
}

static bool refract(const Vec3& v, const Vec3& power, float ni_over_nt,
                    Vec3& refracted) {
  Vec3 uv = Vec3::unit_vector(v);
  float dt = Vec3::dot(uv, power);
  float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if (discriminant > 0) {
    refracted = ni_over_nt * (uv - power * dt) - power * sqrt(discriminant);
    return true;
  } else {
    return false;
  }
}

inline float AbsCosTheta(const Vec3& w) { return fabsf(w.y); }
