#include "assets.h"
#include <cstring>
#include <stdexcept>

using namespace std;

Assets::Assets()
    : logging_category(log4c_category_get("radar.assets"))
{
}

void Assets::configure(const char* sito, time_t acq_time)
{
    if (strcmp(sito, "GAT") == 0)
        conf_site = SITE_GAT;
    else if (strcmp(sito, "SPC") == 0)
        conf_site = SITE_SPC;
    else
    {
        string errmsg(sito);
        throw domain_error(errmsg + " is not a valid radar site name");
    }

    struct tm* tempo = gmtime(&acq_time);
    conf_year = tempo->tm_year + 1900;
    conf_month = tempo->tm_mon + 1;
    conf_day = tempo->tm_mday;
}
