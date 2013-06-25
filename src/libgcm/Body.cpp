#include "Body.h"

string gcm::Body::getId() {
	return id;
}

gcm::Body::Body(string id) {
	this->id = id;
	INIT_LOGGER("gcm.Body");
	LOG_INFO("Body created");
}

gcm::Body::~Body() {
	// clear memory
	foreach(m, meshes)
		delete *m;	
	LOG_INFO("Body destroyed");
}

void gcm::Body::attachMesh(Mesh* mesh) {
	meshes.push_back(mesh);
}

void gcm::Body::setRheology(string rheology) {
	this->rheology = rheology;
}

string gcm::Body::getRheology() {
	return rheology;
}

Mesh* gcm::Body::getMeshes() {
	return meshes.size() ?  meshes[0] : NULL;
}

Mesh* gcm::Body::getMesh(string id) {
	foreach(mesh, meshes)
		if ((*mesh)->getId() == id)
			return (*mesh);
	return NULL;
}

void gcm::Body::setEngine(IEngine* engine) {
	this->engine = engine;
}

IEngine* gcm::Body::getEngine() {
	return engine;
}

void gcm::Body::setInitialState(Area* area, float values[9]) {
	for( unsigned int i = 0; i < meshes.size(); i++ )
		meshes[i]->setInitialState(area, values);
}