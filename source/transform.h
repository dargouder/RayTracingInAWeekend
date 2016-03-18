#pragma once
#include "matrix4x4.h"
#include "ray.h"
class Transform {
public:
	Matrix4x4 matrix, matrix_inverse;

	Transform();

	Transform(float p_matrix[4][4]);

	Transform(Matrix4x4& p_mat);

	static Transform Inverse(const Transform& m_transform);
	Transform(const Matrix4x4& p_matrix, const Matrix4x4& p_matrix_inverse);

	static Transform Translate(const Vec3& p_delta);
	Transform Scale(float p_x, float p_y, float p_z);

	Vec3 operator()(const Vec3& p_vector) const;
	Vec3 DirectionTransform(const Vec3& p_vector) const;
	Vec3 OriginTransform(const Vec3& p_vector) const;
	Ray operator()(const Ray& p_ray) const;

	Transform RotateX(float p_angle);
	Transform RotateY(float p_angle);
	Transform RotateZ(float p_angle);
	Transform LookAt(const Vec3 &pos, const Vec3 &look, const Vec3 &up);

	Transform Perspective(float fov, float n, float f);

	Transform operator*(const Transform &t2) const;
};
