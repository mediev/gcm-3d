#ifndef _GCM_FIXED_BORDER_CALCULATOR_H
#define _GCM_FIXED_BORDER_CALCULATOR_H  1

#include "libgcm/calc/border/BorderCalculator.hpp"
#include <gsl/gsl_linalg.h>

namespace gcm
{
    class FixedBorderCalculator : public BorderCalculator
    {
    public:
        FixedBorderCalculator();
        ~FixedBorderCalculator();
        void doCalc(CalcNode& cur_node, CalcNode& new_node, RheologyMatrixPtr matrix,
                                std::vector<CalcNode>& previousNodes, bool inner[],
                                float outer_normal[], float scale);
        inline std::string getType() {
            return "FixedBorderCalculator";
        }

        void setParameters(const xml::Node& params);

    protected:

    private:
        // Used for border calculation
        gsl_matrix *U_gsl;
        gsl_vector *om_gsl;
        gsl_vector *x_gsl;
        gsl_permutation *p_gsl;
    };
}

#endif
