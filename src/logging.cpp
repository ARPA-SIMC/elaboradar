#include "logging.h"
#include <cstdio>
#include <stdexcept>

using namespace std;

namespace {
static bool initialized = false;
}

/*
 * log4c is currently causing "double free or corruption" if log4c_init is
 * called twice (even with a log4c_fini inbetween).
 *
 * Working around this by calling init just once, and never calling fini
 */

Logging::Logging()
{
    if (!initialized)
    {
        // Initialize logging
        if (log4c_init() != 0)
            throw runtime_error("failed to set up logging");
        initialized = true;
    }
}

Logging::~Logging()
{
    // if (log4c_fini() != 0)
    //   fprintf(stderr, "failed to shut down logging.\n");
}
