#ifndef ARCHIVIATORE_ASSETS_H
#define ARCHIVIATORE_ASSETS_H

#include <string>
#include <ctime>
#include <logging.h>

/**
 * Finds resources, like data files, used by the program.
 */
class Assets
{
protected:
    log4c_category_t* logging_category;
    enum {
        SITE_GAT,
        SITE_SPC,
    } conf_site;
    int conf_year;
    int conf_month;
    int conf_day;

public:
    Assets();

    /**
     * Configure asset lookup with the given details.
     *
     * sito: "GAT" or "SPC" according to the radar site.
     * time: the volume acquisition time
     */
    void configure(const char* sito, time_t acq_time);
};

#endif
