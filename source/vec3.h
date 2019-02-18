#pragma once

#include <cmath>
#include <cstdlib>

#ifdef _WIN32
const float M_PI = 3.142;
#endif

inline double RAND() {
#ifdef _WIN32

  float rt = 0.0;
  rt = (float)(rand() % 1000);
  rt /= 1000;

  return rt;
#else
  return drand48();
#endif
}

inline float schlick(float cosine, float ref_idx) {
  float r0 = (1 - ref_idx) / (1 + ref_idx);
  r0 = r0 * r0;
  float cos_res = 1 - cosine;
  return r0 + (1 - r0) * cos_res * cos_res * cos_res * cos_res * cos_res;
}

//#ifdef __linux__
//#else
// const float M_PI = 3.14159265358979323846;
//#endif
class Vec3 {
 public:
  float e[3];

  Vec3() {}
  Vec3(float e0, float e1, float e2) {
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
  }

  inline float x() const { return e[0]; }
  inline float y() const { return e[1]; }
  inline float z() const { return e[2]; }

  inline float r() const { return e[0]; }
  inline float g() const { return e[1]; }
  inline float b() const { return e[2]; }

  inline const Vec3& operator+() const { return *this; }

  inline Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }

  inline float operator[](int i) const { return e[i]; }

  inline float& operator[](int i) { return e[i]; }

  inline Vec3 operator+(const Vec3& v2) const {
    return Vec3(this->e[0] + v2.e[0], this->e[1] + v2.e[1],
                this->e[2] + v2.e[2]);
  }

  inline friend Vec3 operator-(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
  }

  inline friend Vec3 operator*(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
  }

  inline friend Vec3 operator/(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
  }

  inline friend Vec3 operator/(const Vec3& v1, float t) {
    return Vec3(v1.e[0] / t, v1.e[1] / t, v1.e[2] / t);
  }

  inline friend Vec3 operator*(const Vec3& v1, float t) {
    return Vec3(v1.e[0] * t, v1.e[1] * t, v1.e[2] * t);
  }

  inline friend Vec3 operator*(float t, const Vec3& v1) {
    return Vec3(v1.e[0] * t, v1.e[1] * t, v1.e[2] * t);
  }

  inline Vec3& operator=(const Vec3& v1) {
    e[0] = v1.e[0];
    e[1] = v1.e[1];
    e[2] = v1.e[2];

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

  inline float length() const {
    return sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  }

  inline float squared_length() const {
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
  }

  inline static float dot(const Vec3& v1, const Vec3& v2) {
    return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
  }

  inline static Vec3 cross(const Vec3& v1, const Vec3& v2) {
    return Vec3((v1.e[1] * v2.e[2]) - (v1.e[2] * v2.e[1]),
                -((v1.e[0] * v2.e[2]) - (v1.e[2] * v2.e[0])),
                (v1.e[0] * v2.e[1]) - (v1.e[1] * v2.e[0])

    );
  }

  inline void make_unit_vector() { *this = *this / length(); }

  inline static Vec3 unit_vector(Vec3 v) { return v / v.length(); }

  inline float getLuminance() const {
    return e[0] * 0.212671f + e[1] * 0.715160f + e[2] * 0.072169f;
  }
};

inline Vec3 RandomInUnitDisk() {
  Vec3 p;
  do {
    p = 2.0 * Vec3(RAND(), RAND(), 0) - Vec3(1, 1, 0);
  } while (Vec3::dot(p, p) >= 1.0);

  return p;
}

static Vec3 reflect(const Vec3& v, const Vec3& n) {
  return v - 2 * Vec3::dot(v, n) * n;
}

static Vec3 CosineSampleHemisphere() {
  float u = RAND();
  float v = RAND();
  float z = sqrt(1 - v);
  float phi = 2 * M_PI * u;
  float x = cos(phi) * 2 * sqrt(v);
  float y = sin(phi) * 2 * sqrt(v);
  return Vec3(x, y, z);
}
static Vec3 RandomInUnitSphere() {
  Vec3 p;
  do {
    p = 2.0 * Vec3(RAND(), RAND(), RAND()) - Vec3(1, 1, 1);
  } while (Vec3::dot(p, p) >= 1.0);
  return p;
}

static bool refract(const Vec3& v, const Vec3& n, float ni_over_nt,
                    Vec3& refracted) {
  Vec3 uv = Vec3::unit_vector(v);
  float dt = Vec3::dot(uv, n);
  float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if (discriminant > 0) {
    refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    return true;
  } else {
    return false;
  }
}
