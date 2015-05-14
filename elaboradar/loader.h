#ifndef SRC_LOADER_H
#define SRC_LOADER_H

#include <elaboradar/volume.h>
#include <string>
#include <map>
#include <vector>

namespace elaboradar {
namespace volume {

struct Loader
{
    std::map<std::string, Scans<double>*> to_load;

};


}
}

#endif 

