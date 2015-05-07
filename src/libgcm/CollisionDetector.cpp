#include "libgcm/CollisionDetector.hpp"
#include "libgcm/mesh/tetr/TetrMeshSecondOrder.hpp"
#include "libgcm/node/CalcNode.hpp"

using namespace gcm;
using std::vector;

CollisionDetector::CollisionDetector() {
    INIT_LOGGER("gcm.CollisionDetector");
    static_operation = false;
}

CollisionDetector::~CollisionDetector() {
}

void CollisionDetector::set_static(bool state)
{
    static_operation = state;
}

bool CollisionDetector::is_static()
{
    return static_operation;
}

void CollisionDetector::set_threshold(float value)
{
    threshold = value;
    LOG_DEBUG("Current threshold value: " << threshold);
}

float CollisionDetector::get_threshold()
{
    return threshold;
}

bool CollisionDetector::find_intersection(AABB &outline1, AABB &outline2, AABB &intersection)
{
    // check for intersection
    for(int j = 0; j < 3; j++) {
        intersection.min_coords[j] = fmaxf(outline1.min_coords[j] - threshold, outline2.min_coords[j] - threshold);
        intersection.max_coords[j] = fminf(outline1.max_coords[j] + threshold, outline2.max_coords[j] + threshold);
        if(intersection.min_coords[j] > intersection.max_coords[j])
            return false;
    }

    return true;
}

void CollisionDetector::find_nodes_in_intersection(Mesh* mesh, AABB& intersection, vector<CalcNode>& result)
{
    for(int i = 0; i < mesh->getNodesNumber(); i++)
    {
        CalcNode& node = mesh->getNodeByLocalIndex(i);
        // FIXME
        // only local nodes?
        if ( (node.isLocal ()) && (node.isBorder ()) )
        {
            if(intersection.isInAABB(node))
                result.push_back(node);
        }
    }
}

void CollisionDetector::find_border_nodes_in_intersection(Mesh* mesh, const AABB& intersection, vector<CalcNode>& result)
{
	TetrMeshSecondOrder* soMesh = dynamic_cast<TetrMeshSecondOrder*>(mesh);
	int elemNum = soMesh->borderElements.size();

	for(int i = 0; i < elemNum; i++)
    {
        if(soMesh->borderElements[i].size() != 0) {
			CalcNode& node = soMesh->getNodeByLocalIndex(i);
			if( (node.isLocal ()) && (intersection.isInAABB(node)) )
				result.push_back(soMesh->getNodeByLocalIndex(i));
		}
    }
}

void CollisionDetector::find_nodes_in_intersection_withKD(const Mesh* const mesh, struct kdtree** const kdborder, const AABB& intersection, vector<CalcNode>& result)
{
	double pos;
	const double AABBsize[3] = {intersection.maxX-intersection.minX, intersection.maxY-intersection.minY, intersection.maxZ-intersection.minZ};
	int indx;
	if(AABBsize[0] < AABBsize[1]) {
		if(AABBsize[0] < AABBsize[2])
			indx = 0;
		else
			indx = 2;
	} else if(AABBsize[1] < AABBsize[2])
		indx = 1;
	else
		indx = 2;
	
	struct kdres* set;
	const double intCenter = (intersection.min_coords[indx] + intersection.max_coords[indx]) / 2.0;
	
	set = kd_nearest_range(kdborder[indx], &intCenter, (intersection.max_coords[indx]-intersection.min_coords[indx])/2.0);

	CalcNode* pNode;
	while (!kd_res_end(set)) {	
		pNode = (CalcNode*) kd_res_item(set, &pos);
			
		if( (intersection.isInAABB(pNode)) && (pNode->isLocal()) )
			result.push_back(*pNode);

		kd_res_next(set);
	}
	kd_res_free(set);
}

void CollisionDetector::find_nodes_in_intersection(Mesh* mesh, AABB& intersection, vector<int>& result)
{
    for(int i = 0; i < mesh->getNodesNumber(); i++)
    {
        CalcNode& node = mesh->getNodeByLocalIndex(i);
        // FIXME
        // only local nodes?
        if ( (node.isLocal ()) && (node.isBorder ()) )
        {
            if(intersection.isInAABB(node))
                result.push_back(i);
        }
    }
}

/*void CollisionDetector::renumber_surface(vector<Triangle> &faces, vector<CalcNode> &nodes)
{
    if (!faces.size() || !nodes.size())
        return;
    int max_node_num = -1;
    for(int k = 0; k < nodes.size(); k++)
        if( nodes[k].local_num > max_node_num )
            max_node_num = nodes[k].local_num;
    if(max_node_num < 0)
        return;

    int *renum = new int[max_node_num + 1];
    memset(renum, 0, (max_node_num + 1) * sizeof(int));

    for(int k = 0; k < nodes.size(); k++)
        renum[ nodes[k].local_num ] = k + 1;    // +1 to avoid misinterpreting with zeroed memory

    for(int i = 0; i < faces.size(); i++) {
        for(int j = 0; j < 3; j++) {
            faces[i].vert[j] = renum[ faces[i].vert[j] ] - 1;
            if( faces[i].vert[j] < 0 )
                throw GCMException( GCMException::COLLISION_EXCEPTION, "Can not create correct numbering for surface");
        }
    }

    delete[] renum;
}

void CollisionDetector::renumber_volume(vector<Tetrahedron_1st_order> &tetrs, vector<CalcNode> &nodes)
{
    if (!tetrs.size() || !nodes.size())
        return;
    int max_node_num = -1;
    for(int k = 0; k < nodes.size(); k++)
        if( nodes[k].local_num > max_node_num )
            max_node_num = nodes[k].local_num;
    if(max_node_num < 0)
        return;

    int *renum = new int[max_node_num + 1];
    memset(renum, 0, (max_node_num + 1) * sizeof(int));

    for(int k = 0; k < nodes.size(); k++)
        renum[ nodes[k].local_num ] = k + 1;    // +1 to avoid misinterpreting with zeroed memory

    for(int i = 0; i < tetrs.size(); i++)
    {
        tetrs[i].absolute_num = tetrs[i].local_num;
        tetrs[i].local_num = i;
        for(int j = 0; j < 4; j++) {
            tetrs[i].vert[j] = renum[ tetrs[i].vert[j] ] - 1;
            if( tetrs[i].vert[j] < 0 )
                throw GCMException( GCMException::COLLISION_EXCEPTION, "Can't create correct numbering for volume");
        }
    }

    for(int i = 0; i < nodes.size(); i++)
    {
        nodes[i].absolute_num = nodes[i].local_num;
        nodes[i].local_num = i;
        nodes[i].setPlacement (Local);
    }

    delete[] renum;
}*/
