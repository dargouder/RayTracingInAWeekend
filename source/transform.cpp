#pragma once
#include "transform.h"
#include "vec3.h"

Transform::Transform() {}

Transform::Transform(float p_matrix[4][4]) {
  matrix =
      Matrix4x4(p_matrix[0][0], p_matrix[0][1], p_matrix[0][2], p_matrix[0][3],
                p_matrix[0][0], p_matrix[1][1], p_matrix[1][2], p_matrix[1][3],
                p_matrix[0][0], p_matrix[2][1], p_matrix[2][2], p_matrix[2][3],
                p_matrix[0][0], p_matrix[3][1], p_matrix[3][2], p_matrix[3][3]);

  matrix_inverse = Matrix4x4::Inverse(matrix);
}

Transform::Transform(Matrix4x4 &p_mat)
    : matrix(p_mat), matrix_inverse(Matrix4x4::Inverse(p_mat)) {}

Transform::Transform(const Matrix4x4 &p_matrix,
                     const Matrix4x4 &p_matrix_inverse)
    : matrix(p_matrix), matrix_inverse(p_matrix_inverse) {}

Transform Transform::Inverse(const Transform &m_transform) {
  return Transform(m_transform.matrix_inverse, m_transform.matrix);
}

Transform Transform::Translate(const Vec3 &p_delta) {
  Matrix4x4 m(1, 0, 0, p_delta.x(), 0, 1, 0, p_delta.y(), 0, 0, 1, p_delta.z(),
              0, 0, 0, 1);

  Matrix4x4 minv(1, 0, 0, -p_delta.x(), 0, 1, 0, -p_delta.y(), 0, 0, 1,
                 -p_delta.z(), 0, 0, 0, 1);

  return Transform(m, minv);
}

Transform Transform::Perspective(float fov, float n, float f) {
  // Perform projective divide
  Matrix4x4 persp = Matrix4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, f / (f - n),
                              -f * n / (f - n), 0, 0, 1, 0);

  // Scale to canonical viewing volume
  float invTanAng = 1.f / tanf((fov * M_PI / 180.0) / 2.f);
  return Scale(invTanAng, invTanAng, 1) * Transform(persp);
}

Matrix4x4 Transpose(const Matrix4x4 &p_matrix) {
  return Matrix4x4(
      p_matrix.matrix[0][0], p_matrix.matrix[1][0], p_matrix.matrix[2][0],
      p_matrix.matrix[3][0], p_matrix.matrix[0][1], p_matrix.matrix[1][1],
      p_matrix.matrix[2][1], p_matrix.matrix[3][1], p_matrix.matrix[0][2],
      p_matrix.matrix[1][2], p_matrix.matrix[2][2], p_matrix.matrix[3][2],
      p_matrix.matrix[0][3], p_matrix.matrix[1][3], p_matrix.matrix[2][3],
      p_matrix.matrix[3][3]);
}

Vec3 Transform::operator()(const Vec3 &p_vector) const {
  float x = p_vector.x();
  float y = p_vector.y();
  float z = p_vector.z();

  float xp = matrix.matrix[0][0] * x + matrix.matrix[0][1] * y +
             matrix.matrix[0][2] * z + matrix.matrix[0][3];
  float yp = matrix.matrix[1][0] * x + matrix.matrix[1][1] * y +
             matrix.matrix[1][2] * z + matrix.matrix[1][3];
  float zp = matrix.matrix[2][0] * x + matrix.matrix[2][1] * y +
             matrix.matrix[2][2] * z + matrix.matrix[2][3];
  float wp = matrix.matrix[3][0] * x + matrix.matrix[3][1] * y +
             matrix.matrix[3][2] * z + matrix.matrix[3][3];
  if (wp == 1.) {
    return Vec3(xp, yp, zp);
  } else {
    return Vec3(xp, yp, zp) / wp;
  }
}

Vec3 Transform::DirectionTransform(const Vec3 &p_vector) const {
  float x = p_vector.x();
  float y = p_vector.y();
  float z = p_vector.z();

  float xp = matrix.matrix[0][0] * x + matrix.matrix[0][1] * y +
             matrix.matrix[0][2] * z;
  float yp = matrix.matrix[1][0] * x + matrix.matrix[1][1] * y +
             matrix.matrix[1][2] * z;
  float zp = matrix.matrix[2][0] * x + matrix.matrix[2][1] * y +
             matrix.matrix[2][2] * z;

  return Vec3(xp, yp, zp);
}

Vec3 Transform::OriginTransform(const Vec3 &p_vector) const {
  return (*this)(p_vector);
}

Ray Transform::operator()(const Ray &p_ray) const {
  Ray ret = p_ray;
  ret.A = (*this).OriginTransform(ret.origin());
  ret.B = (*this).DirectionTransform(ret.direction());
  ret.direction().make_unit_vector();
  return ret;
}

Transform Transform::RotateX(float p_angle) {
  float sin_t = sinf(M_PI / 180 * p_angle);
  float cos_t = cosf((M_PI / 180 * p_angle));

  Matrix4x4 m(1, 0, 0, 0, 0, cos_t, -sin_t, 0, 0, sin_t, cos_t, 0, 0, 0, 0, 1);

  return Transform(m, m.Transpose());
}
Transform Transform::RotateY(float p_angle) {
  float sin_t = sinf((M_PI / 180 * p_angle));
  float cos_t = cosf((M_PI / 180 * p_angle));

  Matrix4x4 m(1, 0, 0, 0, 0, cos_t, -sin_t, 0, 0, sin_t, cos_t, 0, 0, 0, 0, 1);

  return Transform(m, m.Transpose());
}
Transform Transform::RotateZ(float p_angle) {

  float sin_t = sinf((M_PI / 180 * p_angle));
  float cos_t = cosf((M_PI / 180 * p_angle));

  Matrix4x4 m(1, 0, 0, 0, 0, cos_t, -sin_t, 0, 0, sin_t, cos_t, 0, 0, 0, 0, 1);

  return Transform(m, m.Transpose());
}

Transform Transform::LookAt(const Vec3 &pos, const Vec3 &look, const Vec3 &up) {
  float m[4][4];
  // Initialize fourth column of viewing matrix
  m[0][3] = pos.x();
  m[1][3] = pos.y();
  m[2][3] = pos.z();
  m[3][3] = 1;

  // Initialize first three columns of viewing matrix
  Vec3 dir, temp_up;
  dir = look - pos;
  dir.make_unit_vector();
  temp_up = up;
  temp_up.make_unit_vector();

  Vec3 left;

  left = Vec3::cross(temp_up, dir);
  left.make_unit_vector();
  Vec3 newUp = Vec3::cross(dir, left);
  m[0][0] = left.x();
  m[1][0] = left.y();
  m[2][0] = left.z();
  m[3][0] = 0.;
  m[0][1] = newUp.x();
  m[1][1] = newUp.y();
  m[2][1] = newUp.z();
  m[3][1] = 0.;
  m[0][2] = dir.x();
  m[1][2] = dir.y();
  m[2][2] = dir.z();
  m[3][2] = 0.;
  Matrix4x4 camToWorld(m);

  return Transform(Matrix4x4::Inverse(camToWorld), camToWorld);
}

Transform Transform::operator*(const Transform &t2) const {
  Matrix4x4 m1 = Matrix4x4::Multiply(matrix, t2.matrix);
  Matrix4x4 m2 = Matrix4x4::Multiply(t2.matrix_inverse, matrix_inverse);
  return Transform(m1, m2);
}

Transform Transform::Scale(float x, float y, float z) {
  Matrix4x4 m(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1);
  Matrix4x4 minv(1.f / x, 0, 0, 0, 0, 1.f / y, 0, 0, 0, 0, 1.f / z, 0, 0, 0, 0,
                 1);
  return Transform(m, minv);
}