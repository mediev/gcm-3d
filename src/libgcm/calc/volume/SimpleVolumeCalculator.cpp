#include "libgcm/calc/volume/SimpleVolumeCalculator.hpp"

#include "libgcm/node/CalcNode.hpp"

SimpleVolumeCalculator::SimpleVolumeCalculator() {
    INIT_LOGGER("gcm.SimpleVolumeCalculator");
};

void SimpleVolumeCalculator::doCalc(CalcNode& new_node, RheologyMatrixPtr matrix,
                                                        vector<CalcNode>& previousNodes)
{
    assert_eq(previousNodes.size(), 10);

    LOG_TRACE("Start calc");

    // Here we will store (omega = Matrix_OMEGA * u)
    float omega[10];

    // Calculate omega value
    // TODO - should we use U and U^-1 from old or new node??? Most probably, the 1st from old and the 2nd from new.
    for(int i = 0; i < 10; i++)
    {
        // omega on new time layer is equal to omega on previous time layer along characteristic
        omega[i] = 0;
        for(int j = 0; j < 10; j++)
        {
            omega[i] += matrix->getU(i,j) * previousNodes[i].values[j];
        }
    }
    // Calculate new values
    for(int i = 0; i < 10; i++)
    {
        new_node.values[i] = 0;
        for(int j = 0; j < 10; j++)
        {
            new_node.values[i] += matrix->getU1(i,j) * omega[j]; }
    }
    LOG_TRACE("Calc done");
};
