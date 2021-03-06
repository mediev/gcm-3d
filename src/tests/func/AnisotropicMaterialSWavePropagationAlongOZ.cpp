#include <vector>
#include <gtest/gtest.h>


#include "tests/func/util.hpp"

#include "libgcm/Exception.hpp"
#include "libgcm/rheology/Material.hpp"
#include "libgcm/node/CalcNode.hpp"

using namespace gcm;

/*
 * Analytics for the test described in tasks/tests/s-wave-test.xml
 * Sets analytical values for CalcNode object provided
 * It does not take into account waves from border, so it works only for wave axis and not too long
 * For high accuracy node should have coordinates (0; 0; z), where -5 < z < 5
 */
void setSWaveAlongZAnalytics(CalcNode& node, float t, Engine& engine)
{
    // Parameters from task file
    float LEFT_MARK_START = 1.0;
    float RIGHT_MARK_START = 3.0;
    float WAVE_AMPLITUDE_SCALE = 0.01;

    const MaterialPtr& mat = engine.getMaterial("trans-isotropic");
    auto p = mat->getRheologyProperties();
    float rho = mat->getRho();

    if (node.coords.z < -5 || node.coords.z > 5)
        THROW_INVALID_INPUT("Z is out of acceptable range");
    if (t < 0 || t > 0.02)
        THROW_INVALID_INPUT("T is out of acceptable range");

    float sWaveXvelocity = sqrt(p.c55 / rho);
    float sWaveYvelocity = sqrt(p.c44 / rho);
    
    float leftMarkX = LEFT_MARK_START - t * sWaveXvelocity;
    float rightMarkX = RIGHT_MARK_START - t * sWaveXvelocity;
	float leftMarkY = LEFT_MARK_START - t * sWaveYvelocity;
    float rightMarkY = RIGHT_MARK_START - t * sWaveYvelocity;

    node.vx = node.vy = node.vz = 0;
    node.sxx = node.sxy = node.sxz = node.syy = node.syz = node.szz = 0;

    if (node.coords.z >= leftMarkX && node.coords.z <= rightMarkX) {
        node.vx = sWaveXvelocity * WAVE_AMPLITUDE_SCALE;
        node.sxz = p.c55 * WAVE_AMPLITUDE_SCALE;
    }
    
    if (node.coords.z >= leftMarkY && node.coords.z <= rightMarkY) {
        node.vy = sWaveYvelocity * WAVE_AMPLITUDE_SCALE;
        node.syz = p.c44 * WAVE_AMPLITUDE_SCALE;
    }
}

TEST(Waves, AnisotropicMaterialSWavePropagationAlongOZ)
{
    // Major test parameters

    // Number of time steps
    int STEPS = 10;

    // Check values and draw graphs along this line
    SnapshotLine line;
    line.startPoint = {0.0, 0.0, -5.0};
    line.endPoint = {0.0, 0.0, 5.0};
    line.numberOfPoints = 100;

    // Thresholds
    float ALLOWED_VALUE_DEVIATION_PERCENT = 0.15;
    int ALLOWED_NUMBER_OF_BAD_NODES = 16; // 2 fronts x 2 nodes per front x 4 without any reason

    runTaskAsTest("tasks/demos/anisotropy-veryfication/s-wave-along-z.xml", setSWaveAlongZAnalytics, STEPS, line, { "vx", "vy", "sxz", "syz" }, ALLOWED_VALUE_DEVIATION_PERCENT, ALLOWED_NUMBER_OF_BAD_NODES);
}
