#ifndef _GCM_MATRIXES_H
#define _GCM_MATRIXES_H  1

#include <iostream>
#include <math.h>
#include <string.h>

#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "libgcm/Math.hpp"
#include "libgcm/Exception.hpp"
#include "libgcm/Logging.hpp"

#define GCM_MATRIX_SIZE 10

namespace gcm {

	class gcm_matrix {
	public:
		gcm_matrix();
		~gcm_matrix();
		gcm_matrix& operator=(const gcm_matrix &A);
		bool operator==(const gcm_matrix &A) const;
		bool operator!=(const gcm_matrix &A) const;
		float& operator()(int i, int j);
		gcm_matrix operator+(const gcm_matrix &A) const;
		gcm_matrix operator-(const gcm_matrix &A) const;
		gcm_matrix operator*(const gcm_matrix &A) const;
		gcm_matrix operator/(const gcm_matrix &A) const;
		gcm_matrix operator*(const real &a) const;
		gcm_matrix operator/(const real &a) const;
		gcm_matrix operator%(const gcm_matrix &A) const;

		float get(unsigned int i, unsigned int j) const;

		float max_abs_value() const;
		void clear();
		void createE();
		void setColumn(double *Clmn, int num);
		gcm_matrix inv() const;

		float p[GCM_MATRIX_SIZE][GCM_MATRIX_SIZE]; // Data
	private:
		USE_LOGGER;
	};
}

namespace std {

	inline std::ostream& operator<<(std::ostream &os, const gcm::gcm_matrix &matrix) {
		for (int r = 0; r < GCM_MATRIX_SIZE; r++) {
			for (int c = 0; c < GCM_MATRIX_SIZE; c++)
				os << matrix.p[r][c] << " ";
			os << std::endl;
		}
		return os;
	};
}

#endif
