#ifndef ARCHIVIATORE_UTILS_H
#define ARCHIVIATORE_UTILS_H

#include <cstdio>
#include <string>
#include "logging.h"
#include "matrix.h"

namespace cumbac {

/**
 * Open a file taking its name from a given env variable.
 *
 * If the variable is not set, assume a 'false' value.
 */
class File
{
protected:
    log4c_category_t* logging_category;
    std::string fname;
    std::string fdesc;
    FILE* fd = nullptr;

public:
    File();
    File(log4c_category_t* logging_category);
    File(const File&) = delete;
    File(File&&);
    ~File();

    File& operator=(const File*) = delete;

    /**
     * Opens a file taking its name from the environment variable envname.
     */
    bool open_from_env(const char* varname, const char* mode, const char* desc=nullptr);

    const char* name() const { return fname.c_str(); }
    const char* desc() const { return fdesc.c_str(); }

    /// Allows FILEFromEnv objects to be used as FILE pointers
    operator FILE*() { return fd; }

    /**
     * Allows FILEFromEnv to be used in an if (...) to check if the file is
     * open
     */
    operator bool() const { return fd; }

};

/**
 * A wrapper of getenv, that returns 'default_value' if the given environment
 * name is not defined.
 */
const char* getenv_default(const char* envname, const char* default_value);

/**
 * A wrapper of fopen that throws an exception if it cannot open the file.
 *
 * If description is provided, it is used in the error message to say what is
 * the file that we were trying to open.
 *
 * It always return a valid, open FILE pointer.
 */
FILE* fopen_checked(const char* fname, const char* mode, const char* description=0);

#endif
