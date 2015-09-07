#include "libgcm/mesh/cube/BasicCubicMeshGenerator.hpp"

#include "libgcm/node/CalcNode.hpp"

using namespace gcm;

BasicCubicMeshGenerator::BasicCubicMeshGenerator() {
    INIT_LOGGER("gcm.BasicCubicMeshGenerator");
}

BasicCubicMeshGenerator::~BasicCubicMeshGenerator() {
}

void BasicCubicMeshGenerator::loadMesh(BasicCubicMesh* mesh, 
	GCMDispatcher* dispatcher, float h, int numX, int numY, int numZ)
{
    for( int k = 0; k <= numZ; k++ )
        for( int j = 0; j <= numY; j++ )
            for( int i = 0; i <= numX; i++ )
            {
				int n = i*(numY+1)*(numZ+1) + j*(numZ+1) + k;
                float x = i*h;
                float y = j*h;
                float z = k*h;
                CalcNode* node = new CalcNode();//(n, x, y, z);
                node->number = n;
                node->coords[0] = x;
                node->coords[1] = y;
                node->coords[2] = z;
                node->setPlacement(true);
                mesh->addNode( *node );
            }
    mesh->preProcess();
}

void BasicCubicMeshGenerator::preLoadMesh(AABB* scene, int& sliceDirection, 
	int& numberOfNodes, float h, int numX, int numY, int numZ)
{
    sliceDirection = 0;
    numberOfNodes = (numX + 1) * (numY + 1) * (numZ + 1);
    scene->minX = scene->minY = scene->minZ = 0;
    scene->maxX = numX * h;
	scene->maxY = numY * h;
	scene->maxZ = numZ * h;
}
