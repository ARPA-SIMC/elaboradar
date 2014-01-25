#include "assets.h"
#include "utils.h"
#include <cstring>
#include <cstdlib>
#include <cerrno>
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

FILE* Assets::open_file_dem()
{
    const char* fname;
    switch (conf_site)
    {
        case SITE_SPC: fname = getenv_default("FILE_DEM_SPC", "../../PP+BLOC/dati/dem_SanPi.txt"); break;
        case SITE_GAT: fname = getenv_default("FILE_DEM_GAT", "../../PP+BLOC/dati/dem_Gatta.txt"); break;
    }
    LOG_INFO("Opening dem file %s", fname);
    return fopen_checked(fname, "rt", "file dem");
}

FILE* Assets::open_file_first_level()
{
    const char* fname = getenv("FIRST_LEVEL_FILE");
    if (!fname)
    {
        switch (conf_site)
        {
            case SITE_SPC:
                if (1 <= conf_month && conf_month <= 3)
                    fname = "../dati/FIRST_LEVEL_SPC_2006_INV";
                else if (4 <= conf_month && conf_month <= 9)
                    fname = "../dati/FIRST_LEVEL_SPC_2006_PRI-EST";
                else
                    fname = "../dati/FIRST_LEVEL_SPC_2006_AUT";
                break;
            case SITE_GAT:
                if (1 <= conf_month && conf_month <= 3)
                    fname = "../dati/FIRST_LEVEL_GAT_2006_INV";
                else if (4 <= conf_month && conf_month <= 9)
                    fname = "../dati/FIRST_LEVEL_GAT_2006_PRI-EST";
                else
                    fname = "../dati/FIRST_LEVEL_GAT_2006_AUT";
                break;
        }
    }
    LOG_INFO("Opening mappa statica %s", fname);
    return fopen_checked(fname, "rb", "mappa statica");
}

int Assets::read_file_first_level_dim()
{
    const char* fname = getenv("FIRST_LEVEL_DIM_FILE");
    if (!fname) throw runtime_error("FIRST_LEVEL_DIM_FILE is not set");
    FILE* fd = fopen_checked(fname, "rt", "dimensioni mappa statica");

    // We can directly read its contents: it's just one number in ASCII
    int dim;
    int res = fscanf(fd, "%i", &dim);
    if (res == EOF)
    {
        string errmsg("Error reading ");
        errmsg += fname;
        errmsg += ": ";
        errmsg += strerror(errno);
        fclose(fd);
        throw runtime_error(errmsg);
    }

    if (res != 1)
    {
        string errmsg("Error reading ");
        errmsg += fname;
        errmsg += ": the file does not seem to contain a number";
        errmsg += strerror(errno);
        fclose(fd);
        throw runtime_error(errmsg);
    }

    fclose(fd);

    return dim;
}
