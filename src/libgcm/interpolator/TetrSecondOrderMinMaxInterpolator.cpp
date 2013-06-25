/* 
 * File:   TetrSecondOrderMinMaxInterpolator.cpp
 * Author: anganar
 * 
 * Created on May 3, 2013, 12:13 AM
 */

#include "TetrSecondOrderMinMaxInterpolator.h"

gcm::TetrSecondOrderMinMaxInterpolator::TetrSecondOrderMinMaxInterpolator() {
	type = "TetrSecondOrderMinMaxInterpolator";
	INIT_LOGGER( "gcm.TetrSecondOrderMinMaxInterpolator" );
}

gcm::TetrSecondOrderMinMaxInterpolator::~TetrSecondOrderMinMaxInterpolator() {
}

void gcm::TetrSecondOrderMinMaxInterpolator::interpolate( ElasticNode* node, 
				ElasticNode* node0, ElasticNode* node1, ElasticNode* node2, ElasticNode* node3,
				ElasticNode* addNode0, ElasticNode* addNode1, ElasticNode* addNode2, 
				ElasticNode* addNode3, ElasticNode* addNode4, ElasticNode* addNode5 )
{
	LOG_TRACE("Start interpolation");
	assert( node != NULL );
	assert( node0 != NULL && node1 != NULL && node2 != NULL && node3 !=NULL );
	assert( addNode0 != NULL && addNode1 != NULL && addNode2 != NULL && addNode3 !=NULL && addNode4 != NULL && addNode5 !=NULL );
	
	float factor[4];
	
	float Vol = tetrVolume(
		(node1->coords[0])-(node0->coords[0]),
		(node1->coords[1])-(node0->coords[1]),
		(node1->coords[2])-(node0->coords[2]),
		(node2->coords[0])-(node0->coords[0]),
		(node2->coords[1])-(node0->coords[1]),
		(node2->coords[2])-(node0->coords[2]),
		(node3->coords[0])-(node0->coords[0]),
		(node3->coords[1])-(node0->coords[1]),
		(node3->coords[2])-(node0->coords[2])
	);

	factor[0] = fabs( tetrVolume(
		(node1->coords[0])-(node->coords[0]),
		(node1->coords[1])-(node->coords[1]),
		(node1->coords[2])-(node->coords[2]),
		(node2->coords[0])-(node->coords[0]),
		(node2->coords[1])-(node->coords[1]),
		(node2->coords[2])-(node->coords[2]),
		(node3->coords[0])-(node->coords[0]),
		(node3->coords[1])-(node->coords[1]),
		(node3->coords[2])-(node->coords[2])
	) / Vol);

	factor[1] = fabs( tetrVolume(
		(node0->coords[0])-(node->coords[0]),
		(node0->coords[1])-(node->coords[1]),
		(node0->coords[2])-(node->coords[2]),
		(node2->coords[0])-(node->coords[0]),
		(node2->coords[1])-(node->coords[1]),
		(node2->coords[2])-(node->coords[2]),
		(node3->coords[0])-(node->coords[0]),
		(node3->coords[1])-(node->coords[1]),
		(node3->coords[2])-(node->coords[2])
	) / Vol);

	factor[2] = fabs( tetrVolume(
		(node1->coords[0])-(node->coords[0]),
		(node1->coords[1])-(node->coords[1]),
		(node1->coords[2])-(node->coords[2]),
		(node0->coords[0])-(node->coords[0]),
		(node0->coords[1])-(node->coords[1]),
		(node0->coords[2])-(node->coords[2]),
		(node3->coords[0])-(node->coords[0]),
		(node3->coords[1])-(node->coords[1]),
		(node3->coords[2])-(node->coords[2])
	) / Vol);

	factor[3] = fabs( tetrVolume(
		(node1->coords[0])-(node->coords[0]),
		(node1->coords[1])-(node->coords[1]),
		(node1->coords[2])-(node->coords[2]),
		(node2->coords[0])-(node->coords[0]),
		(node2->coords[1])-(node->coords[1]),
		(node2->coords[2])-(node->coords[2]),
		(node0->coords[0])-(node->coords[0]),
		(node0->coords[1])-(node->coords[1]),
		(node0->coords[2])-(node->coords[2])
	) / Vol);

	// If we see potential instability
	if(factor[0] + factor[1] + factor[2] + factor[3] > 1.0)
	{
		// If it is small - treat instability as minor and just 'smooth' it
		// TODO - think about it more carefully
		//if( point_in_tetr(node->local_num, node->coords[0], node->coords[1], node->coords[2], tetr) )
		if(factor[0] + factor[1] + factor[2] + factor[3] < 1.1)
		{
			float sum = factor[0] + factor[1] + factor[2] + factor[3];
			for(int i = 0; i < 4; i++)
				factor[i] = factor[i] / sum;
		}
		// If point is not in tetr - throw exception
		else
		{
/*			*logger << "\tTetrVol = " < Vol;
			*logger << "\tfactor[0]=" << factor[0] << " factor[1]=" << factor[1] << " factor[2]=" << factor[2] 	<< " factor[3]=" << factor[3] << " Sum: " < factor[0] + factor[1] + factor[2] + factor[3];

			*logger << "\tnode.x[0]=" << node->coords[0] << " node.x[1]=" << node->coords[1] 
				<< " node.x[2]=" < node->coords[2];
			if( node->isFirstOrder() )
				*logger < "First order node";
			else if( node->isSecondOrder() )
				*logger < "Second order node";

			*logger << "\tv0.x[0]=" << nodes[tetr->vert[0]].coords[0] << " v0.x[1]=" << nodes[tetr->vert[0]].coords[1] << " v0.x[2]=" < nodes[tetr->vert[0]].coords[2];
						
			*logger << "\tv1.x[0]=" << nodes[tetr->vert[1]].coords[0] << " v1.x[1]=" << nodes[tetr->vert[1]].coords[1] << " v1.x[2]=" < nodes[tetr->vert[1]].coords[2];

			*logger << "\tv2.x[0]=" << nodes[tetr->vert[2]].coords[0] << " v2.x[1]=" << nodes[tetr->vert[2]].coords[1] << " v2.x[2]=" < nodes[tetr->vert[2]].coords[2];

			*logger << "\tv3.x[0]=" << nodes[tetr->vert[3]].coords[0] << " v3.x[1]=" << nodes[tetr->vert[3]].coords[1] << " v3.x[2]=" < nodes[tetr->vert[3]].coords[2];*/
			THROW_BAD_MESH("Sum of factors is greater than 1.0");
		}
	}

	baseNodes[0] = node0;
	baseNodes[1] = node1;
	baseNodes[2] = node2;
	baseNodes[3] = node3;
	
	addNodes[0] = addNode0;
	addNodes[1] = addNode1;
	addNodes[2] = addNode2;
	addNodes[3] = addNode3;
	addNodes[4] = addNode4;
	addNodes[5] = addNode5;
	
	for (int i = 0; i < 9; i++)
	{
		float min = baseNodes[0]->values[i];
		float max = baseNodes[0]->values[i];
		for(int z = 1; z < 4; z++) {
			if( baseNodes[z]->values[i] < min )
				min = baseNodes[z]->values[i];
			if( baseNodes[z]->values[i] > max )
				max = baseNodes[z]->values[i];
		}
		for(int z = 0; z < 6; z++) {
			if( addNodes[z]->values[i] < min )
				min = addNodes[z]->values[i];
			if( addNodes[z]->values[i] > max )
				max = addNodes[z]->values[i];
		}

		node->values[i] = ( baseNodes[0]->values[i] * factor[0] * ( 2 * factor[0] - 1 )
					+ baseNodes[1]->values[i] * factor[1] * ( 2 * factor[1] - 1 )
					+ baseNodes[2]->values[i] * factor[2] * ( 2 * factor[2] - 1 )
					+ baseNodes[3]->values[i] * factor[3] * ( 2 * factor[3] - 1 )
					+ addNodes[0]->values[i] * 4 * factor[0] * factor[1]
					+ addNodes[1]->values[i] * 4 * factor[0] * factor[2]
					+ addNodes[2]->values[i] * 4 * factor[0] * factor[3]
					+ addNodes[3]->values[i] * 4 * factor[1] * factor[2]
					+ addNodes[4]->values[i] * 4 * factor[1] * factor[3]
					+ addNodes[5]->values[i] * 4 * factor[2] * factor[3]
					);

		if( node->values[i] < min )
			node->values[i] = min;
		if( node->values[i] > max )
			node->values[i] = max;
	}

	for (int i = 0; i < 3; i++)
	{
		float min = baseNodes[0]->elasticRheologyProperties[i];
		float max = baseNodes[0]->elasticRheologyProperties[i];
		for(int z = 1; z < 4; z++) {
			if( baseNodes[z]->elasticRheologyProperties[i] < min )
				min = baseNodes[z]->elasticRheologyProperties[i];
			if( baseNodes[z]->elasticRheologyProperties[i] > max )
				max = baseNodes[z]->elasticRheologyProperties[i];
		}
		for(int z = 0; z < 6; z++) {
			if( addNodes[z]->elasticRheologyProperties[i] < min )
				min = addNodes[z]->elasticRheologyProperties[i];
			if( addNodes[z]->elasticRheologyProperties[i] > max )
				max = addNodes[z]->elasticRheologyProperties[i];
		}

		node->elasticRheologyProperties[i] = ( baseNodes[0]->elasticRheologyProperties[i] * factor[0] * ( 2 * factor[0] - 1 )
					+ baseNodes[1]->elasticRheologyProperties[i] * factor[1] * ( 2 * factor[1] - 1 )
					+ baseNodes[2]->elasticRheologyProperties[i] * factor[2] * ( 2 * factor[2] - 1 )
					+ baseNodes[3]->elasticRheologyProperties[i] * factor[3] * ( 2 * factor[3] - 1 )
					+ addNodes[0]->elasticRheologyProperties[i] * 4 * factor[0] * factor[1]
					+ addNodes[1]->elasticRheologyProperties[i] * 4 * factor[0] * factor[2]
					+ addNodes[2]->elasticRheologyProperties[i] * 4 * factor[0] * factor[3]
					+ addNodes[3]->elasticRheologyProperties[i] * 4 * factor[1] * factor[2]
					+ addNodes[4]->elasticRheologyProperties[i] * 4 * factor[1] * factor[3]
					+ addNodes[5]->elasticRheologyProperties[i] * 4 * factor[2] * factor[3]
					);

		if( node->elasticRheologyProperties[i] < min )
			node->elasticRheologyProperties[i] = min;
		if( node->elasticRheologyProperties[i] > max )
			node->elasticRheologyProperties[i] = max;
	}
	
	LOG_TRACE("Interpolation done");
}