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

    BruteforceCollisionDetector bcd;
    BruteforceCollisionDetector cd;

    assert_false(bcd.is_static());
    assert_false(cd.is_static());

    auto o1 = m1->getOutline();
    auto o2 = m2->getOutline();

    bcd.find_intersection(o1, o2, intersection);

    vector<CalcNode> nodes1;
    vector<CalcNode> nodes2;
    
    print_test_title("find_nodes_in_intersection");

    auto t = measure_time2(
        [&](){ nodes1.clear(); },
        [&](){ bcd.find_nodes_in_intersection(m1, intersection, nodes1); },
        [&](){ nodes2.clear(); },
        [&](){ cd.find_nodes_in_intersection_withKD(m1, m1->getKDtree(), intersection, nodes2); }
    );

    assert_eq(nodes1.size(), nodes2.size());
    assert_gt(nodes1.size(), 0);
        
    print_test_results("BruteforceCollisionDetector", t.first, "CHAGE NAME", t.second);

    return 0;
}
