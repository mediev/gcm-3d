#include <iostream>
#include "libgcm/config.hpp"

#include "tests/perf/util.hpp"

#include <gmsh/Gmsh.h>

#include <vector>

#if CONFIG_ENABLE_LOGGING
#include <log4cxx/propertyconfigurator.h>
#endif

#include "libgcm/util/KDTree.hpp"
#include "libgcm/mesh/tetr/TetrMeshSecondOrder.hpp"
#include "libgcm/BruteforceCollisionDetector.hpp"

#include "launcher/util/FileFolderLookupService.hpp"
#include "launcher/launcher.hpp"

using namespace gcm;
using namespace std;

int main() {
    MPI::Init();
    GmshInitialize();

    auto& ffls = launcher::FileFolderLookupService::getInstance();
    
    #if CONFIG_ENABLE_LOGGING
    ffls.addPath("src/tests");
    log4cxx::PropertyConfigurator::configure(ffls.lookupFile("log4cxx.properties"));
    #endif
    
    ffls.addPath(".");
    
    cout << "Generating meshes..." << endl;
    cout.flush();

    launcher::Launcher launcher;
    launcher.loadSceneFromFile("tasks/tests/collision-detectors-perf-test.xml");

    auto& engine = Engine::getInstance();

    auto m1 = dynamic_cast<TetrMeshSecondOrder*>(engine.getBodyById("cube1")->getMesh("cube1"));
    auto m2 = dynamic_cast<TetrMeshSecondOrder*>(engine.getBodyById("cube2")->getMesh("cube2"));

    assert_true(m1);
    assert_true(m2);

    assert_gt(m1->getTetrsNumber(), 0);
    assert_gt(m2->getTetrsNumber(), 0);

    cout << "done" << endl;
    
    AABB intersection;

    BruteforceCollisionDetector cd_old;
    cd_old.set_threshold(0.05);
    BruteforceCollisionDetector cd_border;
	cd_border.set_threshold(0.05);
	BruteforceCollisionDetector cd_kd;
	cd_kd.set_threshold(0.05);
	
    assert_false(cd_old.is_static());
    assert_false(cd_border.is_static());
    assert_false(cd_kd.is_static());

    auto o1 = m1->getOutline();
    auto o2 = m2->getOutline();

    cd_old.find_intersection(o1, o2, intersection);
    
    vector<CalcNode> nodes1;
    vector<CalcNode> nodes2;
    vector<CalcNode> nodes3;
    struct kdtree** const kdBorder = m1->getKDborder();
    
    print_test_title("find_nodes_in_intersection");
	
	nodes3.clear();
	cd_old.find_nodes_in_intersection(m1, intersection, nodes3);
	
    auto t = measure_time2(
        [&](){ nodes1.clear(); },
        [&](){ cd_border.find_border_nodes_in_intersection(m1, intersection, nodes1); },
        [&](){ nodes2.clear(); },
        [&](){ cd_kd.find_nodes_in_intersection_withKD(m1, kdBorder, intersection, nodes2); }
    );
      
    print_test_results("BruteforceCollisionDetector_with_border", t.first, "BruteforceCollisionDetector_with_kd", t.second);
    cout << "CD 'with border' found - " << nodes1.size() << " nodes\nCD 'with kd' found - " << nodes2.size() << " nodes" << endl;
    assert_eq(nodes1.size(), nodes2.size());
    assert_eq(nodes2.size(), nodes3.size());
    
    vector<CalcNode>::iterator it1;
    vector<CalcNode>::iterator it2;
    for(int i = 0; i < nodes3.size(); i++) {
		
		it1 = find(nodes1.begin(), nodes1.end(), nodes3[i]);
		if(it1 == nodes1.end()) {
			cout << "Results differ !!!" << endl;
			return -1;
		} else
			nodes1.erase(it1);
			
		it2 = find(nodes2.begin(), nodes2.end(), nodes3[i]);
		if(it2 == nodes2.end()) {
			cout << "Results differ !!!" << endl;
			return -1;
		} else
			nodes2.erase(it2);
	}
	
	assert_eq(nodes1.size(), nodes2.size());
	assert_eq(nodes2.size(), 0);

    return 0;
}
