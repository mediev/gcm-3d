
#include "launcher/util/FileFolderLookupService.hpp"

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;
using std::string;
using std::ifstream;

void launcher::FileFolderLookupService::addPath(string path) {
    paths.push_back(path);
}

string launcher::FileFolderLookupService::lookupFile(string fname) {
    for(auto path: paths) {
        bfs::path p(path);
        p /= fname;
        if (bfs::is_regular_file(p))
            return p.string();
    }
    THROW_INVALID_ARG("File not found: " + fname);
}

string launcher::FileFolderLookupService::lookupFolder(string fname)
{
    for(auto path: paths) {
        bfs::path p(path);
        p /= fname;
        if (bfs::is_directory(p))
            return p.string();
    }
    THROW_INVALID_ARG("Folder not found: " + fname);
}
