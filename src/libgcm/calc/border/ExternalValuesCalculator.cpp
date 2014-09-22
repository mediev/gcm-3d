#include "libgcm/calc/border/ExternalValuesCalculator.hpp"

#include "libgcm/node/CalcNode.hpp"

ExternalValuesCalculator::ExternalValuesCalculator()
{
    U_gsl = gsl_matrix_alloc (10, 10);
    om_gsl = gsl_vector_alloc (10);
    x_gsl = gsl_vector_alloc (10);
    p_gsl = gsl_permutation_alloc (10);
};

ExternalValuesCalculator::~ExternalValuesCalculator()
{
    gsl_matrix_free(U_gsl);
    gsl_vector_free(om_gsl);
    gsl_vector_free(x_gsl);
    gsl_permutation_free(p_gsl);
};

void ExternalValuesCalculator::setParameters(const xml::Node& params)
{
    
};

void ExternalValuesCalculator::set_parameters(int vars[], float vals[])
{
    for(int i = 0; i < 3; i++) {
        vars_index[i] = vars[i];
        vars_values[i] = vals[i];
    }
};

void ExternalValuesCalculator::doCalc(CalcNode& cur_node, CalcNode& new_node, RheologyMatrixPtr& matrix,
                            vector<CalcNode>& previousNodes, bool inner[],
                            float outer_normal[], float scale)
{
	cout << "external values" << endl;
    assert_eq(previousNodes.size(), 10);

    // Tmp value for GSL solver
    int s;

    int outer_count = 3;

    // Here we will store (omega = Matrix_OMEGA * u)
    float omega[10];

    for(int i = 0; i < 10; i++)
    {
        // If omega is 'inner' one
        if(inner[i])
        {
            // Calculate omega value
            omega[i] = 0;
            for(int j = 0; j < 10; j++)
            {
                omega[i] += matrix->getU(i,j) * previousNodes[i].values[j];
            }
            // Load appropriate values into GSL containers
            gsl_vector_set(om_gsl, i, omega[i]);
            for(int j = 0; j < 10; j++)
                gsl_matrix_set(U_gsl, i, j, matrix->getU(i,j));
        }
        // If omega is 'outer' one
        else
        {
            // omega (as right-hand part of OLE) is zero - it is free border, no external stress
            gsl_vector_set(om_gsl, i, 0);
            // corresponding string in matrix is zero ...
            for(int j = 0; j < 10; j++)
                gsl_matrix_set(U_gsl, i, j, 0);

            // ... except normal and tangential stress
            // We use outer normal to find total stress vector (sigma * n) - sum of normal and shear - and tell it is zero
            // TODO - never-ending questions - is everything ok with (x-y-z) and (ksi-eta-dzeta) basises?

            if ( outer_count == 3 ) {
                gsl_matrix_set(U_gsl, i, vars_index[0], 1);
                gsl_vector_set(om_gsl, i, vars_values[0]);
                outer_count--;
            } else if ( outer_count == 2 ) {
                gsl_matrix_set(U_gsl, i, vars_index[1], 1);
                gsl_vector_set(om_gsl, i, vars_values[1]);
                outer_count--;
            } else if ( outer_count == 1 ) {
                gsl_matrix_set(U_gsl, i, vars_index[2], 1);
                gsl_vector_set(om_gsl, i, vars_values[2]);
                outer_count--;
            }
        }
    }

    // Solve linear equations using GSL tools
    gsl_linalg_LU_decomp (U_gsl, p_gsl, &s);
    gsl_linalg_LU_solve (U_gsl, p_gsl, om_gsl, x_gsl);

    for(int j = 0; j < 10; j++)
        new_node.values[j] = gsl_vector_get(x_gsl, j);

};
