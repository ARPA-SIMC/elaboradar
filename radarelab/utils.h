#ifndef RADARELAB_UTILS_H
#define RADARELAB_UTILS_H

#include <functional>
#include <cstdio>
#include <string>
#include "logging.h"
#include "matrix.h"

namespace radarelab {

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
     *
     * Returns false if anything fails.
     */
    bool open_from_env(const char* varname, const char* mode, const char* desc=nullptr);

    /**
     * Opens a file by its pathname.
     *
     * Returns false if anything fails.
     */
    bool open(const std::string& fname, const char* mode, const char* desc=nullptr);

    const char* name() const { return fname.c_str(); }
    const char* desc() const { return fdesc.c_str(); }

    /// Allows FILEFromEnv objects to be used as FILE pointers
    operator FILE*() { return fd; }

    /**
     * Allows FILEFromEnv to be used in an if (...) to check if the file is
     * open
     */
    operator bool() const { return fd; }

    /**
     * Performs a fread on the file, throwing an exception if anything goes
     * wrong.
     *
     * If the read failed because the end of file was reached, it returns
     * false.
     */
    bool fread(void* buf, size_t size);

    void fseek(size_t seek_par, int origin);

    /**
     * Read the file line by line, calling line_cb on each line read
     */
    void read_lines(std::function<void (char*, size_t)> line_cb);
};

/**
 * Split a string in tokens, skipping the given separator characters
 */
void str_split(char* str, const char* sep, std::function<void (const char* tok)> val_cb);

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

}

#endif
