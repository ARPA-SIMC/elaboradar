#include "logging.h"
#include <cstdio>
#include <stdexcept>

using namespace std;

Logging::Logging()
{
    // Initialize logging
    if (log4c_init() != 0)
        throw runtime_error("failed to set up logging");
}

Logging::~Logging()
{
    if (log4c_fini() != 0)
      fprintf(stderr, "failed to shut down logging.\n");
}
