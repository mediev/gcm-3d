#ifndef ANALYTICALRHEOLOGYMATRIXDECOMPOSER10_HPP
#define	ANALYTICALRHEOLOGYMATRIXDECOMPOSER10_HPP

#include "libgcm/rheology/decomposers/AnalyticalRheologyMatrixDecomposer.hpp"
#include "libgcm/rheology/decomposers/NumericalRheologyMatrixDecomposer.hpp"


namespace gcm
{
    class AnalyticalRheologyMatrixDecomposer10: public IDecomposer
    {
        public:
            AnalyticalRheologyMatrixDecomposer10();
            ~AnalyticalRheologyMatrixDecomposer10();
            /**
             * Computes decomposition for matrix in X direction.
             *
             * @param a Matrix to decompose.
             * @param u Matrix to store \f$U\f$
             * @param l Matrix to store \f$\Lambda\f$
             * @param u1 Matrix to store \f$U^{-1}\f$
             */
            void decomposeX(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1) const;
            /**
             * Computes decomposition for matrix in Y direction.
             *
             * @param a Matrix to decompose.
             * @param u Matrix to store \f$U\f$
             * @param l Matrix to store \f$\Lambda\f$
             * @param u1 Matrix to store \f$U^{-1}\f$
             */
            void decomposeY(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1) const;
            /**
             * Computes decomposition for matrix in Z direction.
             *
             * @param a Matrix to decompose.
             * @param u Matrix to store \f$U\f$
             * @param l Matrix to store \f$\Lambda\f$
             * @param u1 Matrix to store \f$U^{-1}\f$
             */
            void decomposeZ(const gcm_matrix& a, gcm_matrix& u, gcm_matrix& l, gcm_matrix& u1) const;
    };
};

#endif	/* ANALYTICALRHEOLOGYMATRIXDECOMPOSER10_HPP */

