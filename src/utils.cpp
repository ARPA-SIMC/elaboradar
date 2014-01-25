#include "utils.h"
#include "logging.h"
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <cstring>
#include <cerrno>

using namespace std;

FILEFromEnv::FILEFromEnv()
    : fd(0)
{
}

FILEFromEnv::FILEFromEnv(const char* envname, const char* mode)
    : fd(0)
{
    open_from_env(envname, mode);
}

bool FILEFromEnv::open_from_env(const char* envname, const char* mode)
{
    LOG_CATEGORY("radar.utils");

    if (fd)
    {
        fclose(fd);
        fd = 0;
        fname.clear();
    }

    if (const char* val = getenv(envname))
    {
        fd = fopen(val, mode);
        if (!fd)
        {
            LOG_ERROR("Cannot open $%s (%s): %s", envname, val, strerror(errno));
        } else {
            fname = val;
        }
    } else {
        LOG_INFO("Cannot open $%s is not set: skipping file.", envname);
    }

    return fd != 0;
}

FILEFromEnv::~FILEFromEnv()
{
    if (fd)
    {
        fclose(fd);
    }
}

const char* getenv_default(const char* envname, const char* default_value)
{
    const char* res = getenv(envname);
    if (res) return res;
    return default_value;
}

FILE* fopen_checked(const char* fname, const char* mode, const char* description)
{
    FILE* res = fopen(fname, mode);
    if (!res)
    {
        string errmsg("cannot open ");
        if (description)
        {
            errmsg += description;
            errmsg += " ";
        }
        errmsg += fname;
        errmsg += " (";
        errmsg += mode;
        errmsg += "): ";
        errmsg += strerror(errno);
        throw runtime_error(errmsg);
    }
    return res;
}
