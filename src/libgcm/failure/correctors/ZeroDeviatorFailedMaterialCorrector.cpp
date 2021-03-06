#include "libgcm/failure/correctors/ZeroDeviatorFailedMaterialCorrector.hpp"

using namespace gcm;

ZeroDeviatorFailedMaterialCorrector::ZeroDeviatorFailedMaterialCorrector() {
    INIT_LOGGER( "gcm.ZeroDeviatorFailedMaterialCorrector" );
}

void ZeroDeviatorFailedMaterialCorrector::applyCorrection(ICalcNode& node, const float tau) {
    if( node.isDestroyed() )
    {
        real p = (node.sxx + node.syy + node.szz) / 3;
        node.sxx = node.syy = node.szz = p;
        node.sxy = node.sxz = node.syz = 0;
    }
}
