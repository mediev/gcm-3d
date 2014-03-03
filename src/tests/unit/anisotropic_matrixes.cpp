#include <time.h>

#include "util/AnisotropicMatrix3D.h"
#include "util/ElasticMatrix3D.h"

#include <gtest/gtest.h>

#define ITERATIONS 1000

using namespace gcm;

class AnisotropicMatrixIsotropicWrapper: public AnisotropicMatrix3D {
public:
	void prepareIsotropic(gcm_real la, gcm_real mu, gcm_real rho, short int stage)
	{
		prepareMatrix( { 
				la+2*mu, la,      la,      0.0,    0.0,    0.0, 
					 la+2*mu, la,      0.0,    0.0,    0.0, 
						  la+2*mu, 0.0,    0.0,    0.0, 
							   mu,     0.0,    0.0, 
								   mu,     0.0, 
									   mu, 
				rho }, stage);
	}
};

class AnisotropicMatrix3DWrapper: public AnisotropicMatrix3D {
public:
		void decomposite(gsl_matrix* a, gsl_matrix* u, gsl_matrix* l, gsl_matrix* u1) {
			decompositeIt(a, u, l, u1);
		}
};

TEST(AnisotropicMatrix3D, FuzzyMultiplication) {
	srand(time(NULL));

	AnisotropicMatrix3DWrapper m;
	
	gsl_matrix* a = gsl_matrix_alloc (9, 9);
	gsl_matrix* u1 = gsl_matrix_alloc (9, 9);
	gsl_matrix* l = gsl_matrix_alloc (9, 9);
	gsl_matrix* u = gsl_matrix_alloc (9, 9);
	gsl_matrix* tmp = gsl_matrix_alloc(9, 9);
	double x;
	int i;
	
	for(i = 0; i < ITERATIONS; i++ )
	{
		for(int j = 0; j < 9; j++)
			for(int k = 0; k <= j; k++) {
				x = pow(10, k)*(double)rand() / RAND_MAX;
				gsl_matrix_set(a, j, k, x);
				gsl_matrix_set(a, k, j, x);
			}
		m.decomposite(a, u1, l, u);
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, u1, l, 0.0, tmp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmp, u, 0.0, u1);
		
		for(int j = 0; j < 9; j++)
			for(int k = 0; k < 9; k++)
				ASSERT_NEAR(gsl_matrix_get(a, j, k), gsl_matrix_get(u1, j, k), 1e-6);
	}
	
	gsl_matrix_free (a);
	gsl_matrix_free (u1);
	gsl_matrix_free (l);
	gsl_matrix_free (u);
	gsl_matrix_free (tmp);
};

TEST(AnisotropicMatrix3D, IsotropicTransition) {
	srand(time(NULL));

	AnisotropicMatrixIsotropicWrapper anisotropicMatrix;
	ElasticMatrix3D isotropicMatrix;

	gcm_real la = 1.0e+9 * (double)rand() / RAND_MAX;
	gcm_real mu = 1.0e+8 * (double)rand() / RAND_MAX;
	gcm_real rho = 1.0e+4 * (double)rand() / RAND_MAX;

	for(int i = 0; i < 3; i++)
	{
		isotropicMatrix.prepareMatrix({ la, mu, rho }, i);
		anisotropicMatrix.prepareIsotropic(la, mu, rho, i);
		
		for(int j = 0; j < 9; j++)
			for(int k = 0; k < 9; k++)
				ASSERT_NEAR(anisotropicMatrix.getA(j, k), isotropicMatrix.getA(j, k), 1e-6);
	}
};