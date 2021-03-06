#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <vector>

#include "libgcm/mesh/tetr/TetrMesh.hpp"
#include "libgcm/elem/TetrFirstOrder.hpp"
#include "libgcm/Engine.hpp"

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

typedef struct {
    std::vector<float> startPoint;
    std::vector<float> endPoint;
    int numberOfPoints;
} SnapshotLine;

typedef struct {
    float min;
    float max;
} ValueLimit;

int mkpath(std::string s);

bfs::path getTestDataDirName();

bfs::path getDataFileName(int stepNum);

gcm::Engine& loadTaskScenario(std::string taskFile);

void dumpPoint(gcm::CalcNode& analytical, gcm::CalcNode& numerical, SnapshotLine line, int stepNum);

bool shouldDraw(std::string value, std::vector<std::string> valuesToDraw);

void drawValues(std::vector<std::string> valuesToDraw, int stepNum, std::vector<ValueLimit> *valueLimits);

void runTaskAsTest(std::string taskFile, void(*setAnalytical)(gcm::CalcNode&, float, gcm::Engine&),
                        int stepsNum, SnapshotLine line, std::vector<std::string> valuesToDraw,
                        float ALLOWED_VALUE_DEVIATION_PERCENT, int ALLOWED_NUMBER_OF_BAD_NODES );

#endif /* UTILS_H_ */
