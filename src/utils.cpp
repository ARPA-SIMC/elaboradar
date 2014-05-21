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
    const char* envfname = getenv(varname);
    if (!envfname)
    {
        LOG_ERROR("$%s is not set", varname);
        return false;
    }

    return open(envfname, mode, desc);
}

bool File::open(const std::string& pathname, const char* mode, const char* desc)
{
    if (fd)
    {
        fclose(fd);
        fd = nullptr;
        fname.clear();
        fdesc.clear();
    }

    fd = fopen(pathname.c_str(), mode);
    if (!fd)
    {
        if (desc)
            LOG_ERROR("Cannot open %s (%s): %s", pathname.c_str(), desc, strerror(errno));
        else
            LOG_ERROR("Cannot open %s: %s", pathname.c_str(), strerror(errno));
        return false;
    }

    fname = pathname;
    if (desc) fdesc = desc;

    return true;
}

void File::read_lines(std::function<void (char*, size_t)> line_cb)
{
    char *line = NULL;
    size_t len = 0;

    while (true)
    {
        errno = 0;
        ssize_t read = getline(&line, &len, fd);
        if (read == -1)
        {
            if (errno == 0)
            {
                break;
            } else {
                string errmsg("cannot read ");
                errmsg += fname;
                errmsg += ": ";
                errmsg += strerror(errno);
                if (line) free(line);
                throw runtime_error(errmsg);
            }
        }
        try {
            line_cb(line, read);
        } catch (...) {
            if (line) free(line);
            throw;
        }
    }

    if (line) free(line);
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

void str_split(char* str, const char* sep, std::function<void (const char* tok)> val_cb)
{
    char* saveptr;
    while (true)
    {
        char* tok = strtok_r(str, sep, &saveptr);
        if (tok == NULL) break;
        val_cb(tok);
        str = NULL;
    }
}

}
