
#pragma once
#include <cstring>
#include <algorithm>

inline float det3(float a, float b, float c,
	float d, float e, float f,
	float g, float h, float i)
{
	return a*e*i + d*h*c + g*b*f - g*e*c - d*b*i - a*h*f;
}
class Matrix4x4 {
public:
	float matrix[4][4];


	Matrix4x4() {
		matrix[0][0] = matrix[1][1] = matrix[2][2] = matrix[3][3] = 1.f;
		matrix[0][1] = matrix[0][2] = matrix[0][3] = matrix[1][0] =
			matrix[1][2] = matrix[1][3] = matrix[2][0] = matrix[2][1] = matrix[2][3] =
			matrix[3][0] = matrix[3][1] = matrix[3][2] = 0.f;
	}
	Matrix4x4(float p_matrix[4][4]) {
		std::memcpy(matrix, p_matrix, 16 * sizeof(float));
	}
	Matrix4x4(float p_t00, float p_t01, float p_t02, float p_t03,
		float p_t10, float p_t11, float p_t12, float p_t13,
		float p_t20, float p_t21, float p_t22, float p_t23,
		float p_t30, float p_t31, float p_t32, float p_t33) {
		matrix[0][0] = p_t00; matrix[0][1] = p_t01; matrix[0][2] = p_t02; matrix[0][3] = p_t03;
		matrix[1][0] = p_t10; matrix[1][1] = p_t11; matrix[1][2] = p_t12; matrix[1][3] = p_t13;
		matrix[2][0] = p_t20; matrix[2][1] = p_t21; matrix[2][2] = p_t22; matrix[2][3] = p_t23;
		matrix[3][0] = p_t30; matrix[3][1] = p_t31; matrix[3][2] = p_t32; matrix[3][3] = p_t33;
	}


	bool operator==(const Matrix4x4 &m2) const {
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				if (matrix[i][j] != m2.matrix[i][j]) return false;
		return true;
	}
	bool operator!=(const Matrix4x4 &m2) const {
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				if (matrix[i][j] != m2.matrix[i][j]) return true;
		return false;
	}


	Matrix4x4 Transpose() {
		return Matrix4x4(matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0],
			matrix[0][1], matrix[1][1], matrix[2][1], matrix[3][1],
			matrix[0][2], matrix[1][2], matrix[2][2], matrix[3][2],
			matrix[0][3], matrix[1][3], matrix[2][3], matrix[3][3]);
	}


	float Determinant()
	{
		float det;
		det = matrix[0][0] * det3(matrix[1][1], matrix[1][2], matrix[1][3],
			matrix[2][1], matrix[2][2], matrix[2][3],
			matrix[3][1], matrix[3][2], matrix[3][3]);
		det -= matrix[0][1] * det3(matrix[1][0], matrix[1][2], matrix[1][3],
			matrix[2][0], matrix[2][2], matrix[2][3],
			matrix[3][0], matrix[3][2], matrix[3][3]);
		det += matrix[0][2] * det3(matrix[1][0], matrix[1][1], matrix[1][3],
			matrix[2][0], matrix[2][1], matrix[2][3],
			matrix[3][0], matrix[3][1], matrix[3][3]);
		det -= matrix[0][3] * det3(matrix[1][0], matrix[1][1], matrix[1][2],
			matrix[2][0], matrix[2][1], matrix[2][2],
			matrix[3][0], matrix[3][1], matrix[3][2]);
		return det;
	}

	static Matrix4x4 Inverse(const Matrix4x4 &m) {
		int indxc[4], indxr[4];
		int ipiv[4] = { 0, 0, 0, 0 };
		float minv[4][4];
		memcpy(minv, m.matrix, 4 * 4 * sizeof(float));
		for (int i = 0; i < 4; i++) {
			int irow = -1, icol = -1;
			float big = 0.;
			// Choose pivot
			for (int j = 0; j < 4; j++) {
				if (ipiv[j] != 1) {
					for (int k = 0; k < 4; k++) {
						if (ipiv[k] == 0) {
							if (fabsf(minv[j][k]) >= big) {
								big = float(fabsf(minv[j][k]));
								irow = j;
								icol = k;
							}
						}
						else if (ipiv[k] > 1)
							return 0;
					}
				}
			}
			++ipiv[icol];
			// Swap rows _irow_ and _icol_ for pivot
			if (irow != icol) {
				for (int k = 0; k < 4; ++k)
					std::swap(minv[irow][k], minv[icol][k]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (minv[icol][icol] == 0.)
				return 0;

			// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
			float pivinv = 1.f / minv[icol][icol];
			minv[icol][icol] = 1.f;
			for (int j = 0; j < 4; j++)
				minv[icol][j] *= pivinv;

			// Subtract this row from others to zero out their columns
			for (int j = 0; j < 4; j++) {
				if (j != icol) {
					float save = minv[j][icol];
					minv[j][icol] = 0;
					for (int k = 0; k < 4; k++)
						minv[j][k] -= minv[icol][k] * save;
				}
			}
		}
		// Swap columns to reflect permutation
		for (int j = 3; j >= 0; j--) {
			if (indxr[j] != indxc[j]) {
				for (int k = 0; k < 4; k++)
					std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
			}
		}
		return Matrix4x4(minv);
	}

	static Matrix4x4 Multiply(const Matrix4x4 &m1, const Matrix4x4 &m2) {
		Matrix4x4 r;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				r.matrix[i][j] = m1.matrix[i][0] * m2.matrix[0][j] +
				m1.matrix[i][1] * m2.matrix[1][j] +
				m1.matrix[i][2] * m2.matrix[2][j] +
				m1.matrix[i][3] * m2.matrix[3][j];
		return r;
	}
};