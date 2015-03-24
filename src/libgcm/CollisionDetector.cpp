#include "libgcm/CollisionDetector.hpp"

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

void CollisionDetector::find_nodes_in_intersection_withKD(const Mesh* const mesh, struct kdtree* const kd, const AABB& intersection, vector<CalcNode>& result)
{
	
	// Find sizes of AABB intersection
	const double AABBsize[3] = {intersection.maxX-intersection.minX, intersection.maxY-intersection.minY, intersection.maxZ-intersection.minZ};
	
	// .. and find the minimum size
	int indx[3];
	if(AABBsize[0] < AABBsize[1]) {
		indx[1] = 1;
		
		if(AABBsize[0] < AABBsize[2]) {
			indx[0] = 0;
			indx[2] = 2;
		} else {
			indx[0] = 2;
			indx[2] = 0;
		}		
	} else if(AABBsize[1] < AABBsize[2]) {
		indx[0] = 1;
		indx[1] = 0;
		indx[2] = 2;
	} else {
		indx[0] = 2;
		indx[1] = 1;
		indx[2] = 0;
	}
	
	// Number of spheres along two axis with non-minimum AABB intersection size
	const double sizesRatio[2] = {ceil(AABBsize[indx[1]]/AABBsize[indx[0]]), ceil(AABBsize[indx[2]]/AABBsize[indx[0]])};
	const double rad = AABBsize[indx[0]] * sqrt(3.0) / 2.0;
	
	vector<CalcNode>::iterator it;
	struct kdres* set;
	double pt[3];
	double pos[3];
	CalcNode* pNode;
			
	// Filling vector of vertices which are situated in AABB intersection 
	pt[indx[0]] = intersection.min_coords[indx[0]] + AABBsize[indx[0]] / 2.0;
	for(int i = 0; i < sizesRatio[0]; i++)
	{
		pt[indx[1]] = intersection.min_coords[indx[1]] + (i + 0.5) * AABBsize[indx[0]];
		for(int j = 0; j < sizesRatio[1]; j++)
		{
			pt[indx[2]] = intersection.min_coords[indx[2]] + (j + 0.5) * AABBsize[indx[0]];	
			set = kd_nearest_range(kd, pt, rad);
			
			while (!kd_res_end(set)) {	
				pNode = (CalcNode*) kd_res_item(set, pos);
				
				if( (intersection.isInAABB(pNode)) && (pNode->isLocal()) && (pNode->isBorder()) )
				{
					// Checking that we haven't found it yet
					it = find (result.begin(), result.end(), *pNode);
					if(it == result.end())
						result.push_back(*pNode);
				}
				kd_res_next(set);
			}
			kd_res_free(set);
		}
	
	}
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
