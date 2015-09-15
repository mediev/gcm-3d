#ifndef GCM_BASIC_CUBIC_MESH_H_
#define GCM_BASIC_CUBIC_MESH_H_

#include "libgcm/mesh/Mesh.hpp"
#include "libgcm/Math.hpp"
#include "libgcm/Logging.hpp"
#include "libgcm/Exception.hpp"

#include "libgcm/interpolator/LineFirstOrderInterpolator.hpp"
#include "libgcm/interpolator/LineSecondOrderInterpolator.hpp"


#include <utility>
#include <algorithm>

namespace gcm
{
    class CalcNode;

    class BasicCubicMesh: public Mesh {

    // Hopefully we will never need these 'friends' back
/*    friend class VTKSnapshotWriter;
    friend class DataBus;
    friend class CollisionDetector;
    friend class BruteforceCollisionDetector;
*/

    private:
        LineFirstOrderInterpolator* interpolator1;
		LineSecondOrderInterpolator* interpolator2;

    protected:
        void logMeshStats();
        void calcMinH();
        void preProcessGeometry();

        int findNeighbourPoint(CalcNode& node, float dx, float dy, float dz, bool debug, float* coords, bool* innerPoint);

        void findNearestsNodes(const vector3r& coords, int N, std::vector< std::pair<int,float> >& result);

        float meshH;

        // Number of cubes along axis
        uint numX, numY, numZ;

        USE_LOGGER;

    public:
        BasicCubicMesh();
        ~BasicCubicMesh();

        float getRecommendedTimeStep();
        float getMinH();
		float getAvgH() override;
        void doNextPartStep(float tau, int stage);
        void checkTopology(float tau);

        void findBorderNodeNormal(const CalcNode& node, float* x, float* y, float* z, bool debug);

        bool interpolateNode(CalcNode& origin, float dx, float dy, float dz, bool debug,
                                CalcNode& targetNode, bool& isInnerPoint);

        bool interpolateNode(CalcNode& node);

        bool interpolateBorderNode(real x, real y, real z, 
                                real dx, real dy, real dz, CalcNode& node);

        bool interpolateBorderNode(const vector3r& x, const vector3r& dx, CalcNode& node);

        void setNumX(uint _numX);
        void setNumY(uint _numY);
        void setNumZ(uint _numZ);

        uint getNumX() const;
        uint getNumY() const;
        uint getNumZ() const;

        virtual const SnapshotWriter& getSnaphotter() const override;
        virtual const SnapshotWriter& getDumper() const override;
    };
}
#endif
