#include "libgcm/calc/contact/SlidingContactCalculator.hpp"

#include "libgcm/node/CalcNode.hpp"

SlidingContactCalculator::SlidingContactCalculator()
{
    fbc = new FreeBorderCalculator();
    U_gsl = gsl_matrix_alloc (20, 20);
    om_gsl = gsl_vector_alloc (20);
    x_gsl = gsl_vector_alloc (20);
    p_gsl = gsl_permutation_alloc (20);
};

SlidingContactCalculator::~SlidingContactCalculator()
{
    delete fbc;
    gsl_matrix_free(U_gsl);
    gsl_vector_free(om_gsl);
    gsl_vector_free(x_gsl);
    gsl_permutation_free(p_gsl);
};

void SlidingContactCalculator::doCalc(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node,
                            RheologyMatrixPtr matrix, vector<CalcNode>& previousNodes, bool inner[],
                            RheologyMatrixPtr virt_matrix, vector<CalcNode>& virtPreviousNodes, bool virt_inner[],
                            float outer_normal[], float scale)
{
    assert_eq(previousNodes.size(), 10);
    assert_eq(virtPreviousNodes.size(), 10);

    float local_n[3][3];
    local_n[0][0] = outer_normal[0];
    local_n[0][1] = outer_normal[1];
    local_n[0][2] = outer_normal[2];

    createLocalBasis(local_n[0], local_n[1], local_n[2]);

    //---------------------------------------Check if nodes fall apart
    // TODO - may be '-'?
    // FIXME_ASAP - it's regarding normals...
    float vel_rel[3] = {
        cur_node.vx + virt_node.vx,
        cur_node.vy + virt_node.vy,
        cur_node.vz + virt_node.vz
    };

    float force_cur[3] = {
        cur_node.sxx*outer_normal[0] + cur_node.sxy*outer_normal[1] + cur_node.sxz*outer_normal[2],
        cur_node.sxy*outer_normal[0] + cur_node.syy*outer_normal[1] + cur_node.syz*outer_normal[2],
        cur_node.sxz*outer_normal[0] + cur_node.syz*outer_normal[1] + cur_node.szz*outer_normal[2]
    };

    float vel_abs = -scalarProduct(vel_rel, outer_normal);
    float force_cur_abs = scalarProduct(force_cur,outer_normal);

    // TODO - remove magic number
    float eps = 0.0005;
    bool free_border = false;

    if (vel_abs < -eps)             //first check relative speed
        free_border = true;
    else if (vel_abs < eps)            //if relative speed is small, we check force
        if (force_cur_abs > -eps)
            free_border = true;

    if (free_border)
    {
        fbc->doCalc(cur_node, new_node, matrix, previousNodes, inner, outer_normal, scale);
        return;
    }
    //--------------------------------------------------------------------


    // Here we will store (omega = Matrix_OMEGA * u)
    float omega[10];
    float virt_omega[10];

    int posInEq18 = 0;
    int curNN = 0;

    // For all omegas of real node
    for(int i = 0; i < 10; i++)
    {
        // If omega is 'inner'
        if(inner[i])
        {
            // omega on new time layer is equal to omega on previous time layer along characteristic
            omega[i] = 0;
            for( int j = 0; j < 10; j++ ) {
                omega[i] += matrix->getU(i,j) * previousNodes[i].values[j];
            }

            // then we must set the corresponding values of the 18x18 matrix
            gsl_vector_set( om_gsl, 7 * curNN + posInEq18, omega[i] );

            for( int j = 0; j < 10; j++ ) {
                gsl_matrix_set( U_gsl, 7 * curNN + posInEq18, j, matrix->getU( i, j ) );
            }
            for( int j = 10; j < 20; j++ ) {
                gsl_matrix_set( U_gsl, 7 * curNN + posInEq18, j, 0 );
            }
            posInEq18++;
        }
    }

    posInEq18 = 0;
    curNN = 1;
    // For all omegas of virtual node
    for(int i = 0; i < 10; i++)
    {
        // If omega is 'inner'
        if(virt_inner[i])
        {
            // omega on new time layer is equal to omega on previous time layer along characteristic
            virt_omega[i] = 0;
            for( int j = 0; j < 10; j++ ) {
                virt_omega[i] += virt_matrix->getU(i,j) * virtPreviousNodes[i].values[j];
            }

            // then we must set the corresponding values of the 18x18 matrix
            gsl_vector_set( om_gsl, 7 * curNN + posInEq18, virt_omega[i] );

            for( int j = 0; j < 10; j++ ) {
                gsl_matrix_set( U_gsl, 7 * curNN + posInEq18, j, 0 );
            }
            for( int j = 10; j < 20; j++ ) {
                gsl_matrix_set( U_gsl, 7 * curNN + posInEq18, j, virt_matrix->getU( i, j - 10 ) );
            }
            posInEq18++;
        }
    }

    // Clear the rest 6 rows of the matrix
    for( int strN = 14; strN < 20                       ; strN++ ) {
        for( int colN = 0; colN < 20; colN++ ) {
            gsl_matrix_set( U_gsl, strN, colN, 0 );
        }
    }

    for( int strN = 14; strN < 20; strN++ ) {
        gsl_vector_set( om_gsl, strN, 0 );
    }

    // Normal velocities are equal
    gsl_matrix_set( U_gsl, 14, 0, local_n[0][0]);
    gsl_matrix_set( U_gsl, 14, 1, local_n[0][1]);
    gsl_matrix_set( U_gsl, 14, 2, local_n[0][2]);
    gsl_matrix_set( U_gsl, 14, 10,  - local_n[0][0]);
    gsl_matrix_set( U_gsl, 14, 11, - local_n[0][1]);
    gsl_matrix_set( U_gsl, 14, 12, - local_n[0][2]);

    // We use outer normal to find total stress vector (sigma * n) - sum of normal and shear - and tell it is equal
    // TODO - is it ok?
    // TODO - never-ending questions - is everything ok with (x-y-z) and (ksi-eta-dzeta) basises?

    // TODO FIXME - it works now because exactly the first axis is the only one where contact is possible
    // and it coincides with outer normal

    // Normal stresses are equal
    gsl_matrix_set(U_gsl, 15, 3, local_n[0][0] * local_n[0][0]);
    gsl_matrix_set(U_gsl, 15, 4, 2 * local_n[0][1] * local_n[0][0]);
    gsl_matrix_set(U_gsl, 15, 5, 2 * local_n[0][2] * local_n[0][0]);
    gsl_matrix_set(U_gsl, 15, 6, local_n[0][1] * local_n[0][1]);
    gsl_matrix_set(U_gsl, 15, 7, 2 * local_n[0][2] * local_n[0][1]);
    gsl_matrix_set(U_gsl, 15, 8, local_n[0][2] * local_n[0][2]);

    gsl_matrix_set(U_gsl, 15, 13, - local_n[0][0] * local_n[0][0]);
    gsl_matrix_set(U_gsl, 15, 14, - 2 * local_n[0][1] * local_n[0][0]);
    gsl_matrix_set(U_gsl, 15, 15, - 2 * local_n[0][2] * local_n[0][0]);
    gsl_matrix_set(U_gsl, 15, 16, - local_n[0][1] * local_n[0][1]);
    gsl_matrix_set(U_gsl, 15, 17, - 2 * local_n[0][2] * local_n[0][1]);
    gsl_matrix_set(U_gsl, 15, 18, - local_n[0][2] * local_n[0][2]);

    // Tangential stresses are zero

    gsl_matrix_set(U_gsl, 16, 3, - (local_n[0][0] * local_n[1][0]) );
    gsl_matrix_set(U_gsl, 16, 4, - (local_n[0][1] * local_n[1][0] + local_n[0][0] * local_n[1][1]) );
    gsl_matrix_set(U_gsl, 16, 5, - (local_n[0][2] * local_n[1][0] + local_n[0][0] * local_n[1][2]) );
    gsl_matrix_set(U_gsl, 16, 6, - (local_n[0][1] * local_n[1][1]) );
    gsl_matrix_set(U_gsl, 16, 7, - (local_n[0][2] * local_n[1][1] + local_n[0][1] * local_n[1][2]) );
    gsl_matrix_set(U_gsl, 16, 8, - (local_n[0][2] * local_n[1][2]) );

    gsl_matrix_set(U_gsl, 17, 3, - (local_n[0][0] * local_n[2][0]) );
    gsl_matrix_set(U_gsl, 17, 4, - (local_n[0][1] * local_n[2][0] + local_n[0][0] * local_n[2][1]) );
    gsl_matrix_set(U_gsl, 17, 5, - (local_n[0][2] * local_n[2][0] + local_n[0][0] * local_n[2][2]) );
    gsl_matrix_set(U_gsl, 17, 6, - (local_n[0][1] * local_n[2][1]) );
    gsl_matrix_set(U_gsl, 17, 7, - (local_n[0][2] * local_n[2][1] + local_n[0][1] * local_n[2][2]) );
    gsl_matrix_set(U_gsl, 17, 8, - (local_n[0][2] * local_n[2][2]) );


    gsl_matrix_set(U_gsl, 18, 13, - (local_n[0][0] * local_n[1][0]) );
    gsl_matrix_set(U_gsl, 18, 14, - (local_n[0][1] * local_n[1][0] + local_n[0][0] * local_n[1][1]) );
    gsl_matrix_set(U_gsl, 18, 15, - (local_n[0][2] * local_n[1][0] + local_n[0][0] * local_n[1][2]) );
    gsl_matrix_set(U_gsl, 18, 16, - (local_n[0][1] * local_n[1][1]) );
    gsl_matrix_set(U_gsl, 18, 17, - (local_n[0][2] * local_n[1][1] + local_n[0][1] * local_n[1][2]) );
    gsl_matrix_set(U_gsl, 18, 18, - (local_n[0][2] * local_n[1][2]) );

    gsl_matrix_set(U_gsl, 19, 13, - (local_n[0][0] * local_n[2][0]) );
    gsl_matrix_set(U_gsl, 19, 14, - (local_n[0][1] * local_n[2][0] + local_n[0][0] * local_n[2][1]) );
    gsl_matrix_set(U_gsl, 19, 15, - (local_n[0][2] * local_n[2][0] + local_n[0][0] * local_n[2][2]) );
    gsl_matrix_set(U_gsl, 19, 16, - (local_n[0][1] * local_n[2][1]) );
    gsl_matrix_set(U_gsl, 19, 17, - (local_n[0][2] * local_n[2][1] + local_n[0][1] * local_n[2][2]) );
    gsl_matrix_set(U_gsl, 19, 18, - (local_n[0][2] * local_n[2][2]) );


    // Tmp value for GSL solver
    int s;
    gsl_linalg_LU_decomp (U_gsl, p_gsl, &s);
    gsl_linalg_LU_solve (U_gsl, p_gsl, om_gsl, x_gsl);

    // Just get first 9 values (real node) and dump the rest 9 (virt node)
    for(int j = 0; j < 10; j++)
        new_node.values[j] = gsl_vector_get(x_gsl, j);

};
