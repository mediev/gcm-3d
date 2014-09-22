#include "libgcm/rheology/decomposers/NumericalRheologyMatrixDecomposer10.hpp"

#include <cmath>

NumericalRheologyMatrixDecomposer10::NumericalRheologyMatrixDecomposer10()
{
    _a = gsl_matrix_alloc(10, 10);
    _u = gsl_matrix_alloc(10, 10);
    _u1 = gsl_matrix_alloc(10, 10);
    
    eval = gsl_vector_complex_alloc(10);
    evec = gsl_matrix_complex_alloc(10, 10);

    perm = gsl_permutation_alloc(10);
    w = gsl_eigen_nonsymmv_alloc(10);
}

NumericalRheologyMatrixDecomposer10::~NumericalRheologyMatrixDecomposer10()
{
    gsl_matrix_free(_a);
    gsl_matrix_free(_u);
    gsl_matrix_free(_u1);

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    
    gsl_permutation_free(perm);
    gsl_eigen_nonsymmv_free(w);
}

void NumericalRheologyMatrixDecomposer10::gsl2gcm(const gsl_matrix* a, gcm_matrix& b) const
{
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            b(i, j) = gsl_matrix_get(a, i, j);
};

void NumericalRheologyMatrixDecomposer10::gcm2gsl(const gcm_matrix& a, gsl_matrix* b) const
{
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            gsl_matrix_set(b, i, j, a.get(i, j));
};
void NumericalRheologyMatrixDecomposer10::decompose(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1, int stage) const
{
    const float IMAG_THRESHOLD = 1e-8;

    gcm2gsl(a, _a);
    gsl_eigen_nonsymmv(_a, eval, evec, w);
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    l.clear();
    for (int i = 0; i < 10; i++)
    {
        gsl_complex z = gsl_vector_complex_get(eval, i);
        if (fabs(GSL_IMAG(z)) > IMAG_THRESHOLD) {
//			cout << a << endl << endl;
//			for(int k = 0; k < 10; k++)
//				cout << GSL_REAL(gsl_vector_complex_get(eval, k)) << "+" << GSL_IMAG(gsl_vector_complex_get(eval, k)) << "i\n";
            THROW_INVALID_INPUT("Eigenvalue is complex!");
		}
        l(i, i) = GSL_REAL(z);
    }
    
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            gsl_complex z = gsl_matrix_complex_get(evec, i, j);
            if (fabs(GSL_IMAG(z)) > IMAG_THRESHOLD) {
//				cout << endl << a << endl << endl;
//				for(int k = 0; k < 10; k++)
//					cout << GSL_REAL(gsl_vector_complex_get(eval, k)) << "+" << GSL_IMAG(gsl_vector_complex_get(eval, k)) << "i\n";
//				for(int k = 0; k < 10; k++) {
//					for(int l = 0; l < 10; l++)
//						cout << GSL_REAL(gsl_matrix_complex_get(evec, k, l)) << "+" << GSL_IMAG(gsl_matrix_complex_get(evec, k, l)) << "i\t\t";
//					cout << endl;
//				}
                THROW_INVALID_INPUT("Eigenvector is complex!");
			}
            u1(i, j) = GSL_REAL(z);
        }
        for (int j = 6; j < 10; j++)
            u1(i, j) =  0.0;
    }

    switch (stage)
    {
        case 0: u1(6, 6) = u1(7, 7) = u1(8, 8) = u1(9, 9) = 1.0; break;
        case 1: u1(3, 6) = u1(5, 7) = u1(8, 8) = u1(9, 9) = 1.0; break;
        case 2: u1(3, 6) = u1(4, 7) = u1(6, 8) = u1(9, 9) = 1.0; break;

    }

    gcm2gsl(u1, _u1);

    int s;
    gsl_linalg_LU_decomp(_u1, perm, &s);
    gsl_linalg_LU_invert(_u1, perm, _u);

    gsl2gcm(_u, u);
}

void NumericalRheologyMatrixDecomposer10::decomposeX(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1) const
{
    decompose(a, u, l, u1, 0);
}

void NumericalRheologyMatrixDecomposer10::decomposeY(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1) const
{
    decompose(a, u, l, u1, 1);
}

void NumericalRheologyMatrixDecomposer10::decomposeZ(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1) const
{
    decompose(a, u, l, u1, 2);
}

