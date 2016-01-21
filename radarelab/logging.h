/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_LOGGING_H
#define RADARELAB_LOGGING_H

extern "C" {
#include <log4c.h>
}

/*
 * Convenient macros for logging. They assume that a 'logging_category'
 * variable of type `const log4c_category_t*` is accessible in the current
 * scope.
 */

#define LOG_DEBUG(...) log4c_category_log(logging_category, LOG4C_PRIORITY_DEBUG, __VA_ARGS__)
#define LOG_INFO(...) log4c_category_log(logging_category, LOG4C_PRIORITY_INFO, __VA_ARGS__)
#define LOG_WARN(...) log4c_category_log(logging_category, LOG4C_PRIORITY_WARN, __VA_ARGS__)
#define LOG_ERROR(...) log4c_category_log(logging_category, LOG4C_PRIORITY_ERROR, __VA_ARGS__)

/**
 * Define and initialize a logging_category variable for the given category
 * name
 */
#define LOG_CATEGORY(name) log4c_category_t* logging_category = log4c_category_get(name)

class Logging
{
public:
    /**
     * Set up log4c.
     *
     * Throws runtime_error if it fails.
     */
    Logging();

    /**
     * Shut down log4c.
     *
     * Prints a message on stderr if it fails.
     */
    ~Logging();
};

#endif
