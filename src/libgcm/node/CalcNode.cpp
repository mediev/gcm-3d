#include "CalcNode.h"

gcm::CalcNode::CalcNode()
{
	la = 0;
	mu = 0;
	rho = 0;
	memset(values, 0, 9*sizeof(float));
	publicFlags = 0;
	privateFlags = 0;
	errorFlags = 0;
	// Border condition '1' is 'default' one, since '0' is reserved for 'failsafe' one
	borderCondId = 1;
	addOwner( GCM );
	elements = new vector<int>;
	border_elements = new vector<int>;
}

gcm::CalcNode::CalcNode(int _num) {
	
}

gcm::CalcNode::CalcNode(int _num, float _x, float _y, float _z) {
	
}

gcm::CalcNode::CalcNode(const CalcNode& src) {
	number = src.number;
	memcpy( coords, src.coords, 3*sizeof(float) );
	memcpy( values, src.values, 9*sizeof(float) );
	memcpy( elasticRheologyProperties, src.elasticRheologyProperties, 3*sizeof(float) );
	publicFlags = src.publicFlags;
	privateFlags = src.privateFlags;
	errorFlags = src.errorFlags;
	borderCondId = src.borderCondId;
	elements = new vector<int>;
	border_elements = new vector<int>;
	for( unsigned int i = 0; i < src.elements->size(); i++ )
		elements->push_back( src.elements->at(i) );
	for( unsigned int i = 0; i < src.border_elements->size(); i++ )
		border_elements->push_back( src.border_elements->at(i) );
}

CalcNode& gcm::CalcNode::operator=(const CalcNode &src)
{
	number = src.number;
	memcpy( coords, src.coords, 3*sizeof(float) );
	memcpy( values, src.values, 9*sizeof(float) );
	memcpy( elasticRheologyProperties, src.elasticRheologyProperties, 3*sizeof(float) );
	publicFlags = src.publicFlags;
	privateFlags = src.privateFlags;
	errorFlags = src.errorFlags;
	borderCondId = src.borderCondId;
	elements = new vector<int>;
	border_elements = new vector<int>;
	for( unsigned int i = 0; i < src.elements->size(); i++ )
		elements->push_back( src.elements->at(i) );
	for( unsigned int i = 0; i < src.border_elements->size(); i++ )
		border_elements->push_back( src.border_elements->at(i) );
	return *this;
}

gcm::CalcNode::~CalcNode()
{
	elements->clear();
	delete elements;
	border_elements->clear();
	delete border_elements;
}

void gcm::CalcNode::clearErrorFlags()
{
	errorFlags = 0;
}