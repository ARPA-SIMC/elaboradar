#include "utils.h"
#include "logging.h"
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <cstring>
#include <cerrno>

using namespace std;

namespace cumbac {

File::File()
    : logging_category(log4c_category_get("radar.utils"))
{
}

File::File(File&& f)
    : logging_category(f.logging_category),
      fname(std::move(f.fname)),
      fdesc(std::move(f.fdesc)),
      fd(f.fd)
{
    f.fd = 0;
}

File::File(log4c_category_t* logging_category)
    : logging_category(logging_category)
{
}

File::~File()
{
    if (fd) fclose(fd);
}

bool File::open_from_env(const char* varname, const char* mode, const char* desc)
{
    if (fd)
    {
        fclose(fd);
        fd = nullptr;
        fname.clear();
        fdesc.clear();
    }

    const char* envfname = getenv(varname);
    if (!envfname)
    {
        LOG_ERROR("$%s is not set", varname);
        return false;
    }

    fd = fopen(envfname, mode);
    if (!fd)
    {
        LOG_ERROR("Cannot open $%s=%s: %s", varname, envfname, strerror(errno));
        return false;
    }

    fname = envfname;
    if (desc) fdesc = desc;

    return true;
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
