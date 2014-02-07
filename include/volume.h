#ifndef ARCHIVIATORE_VOLUME_CLASS_H
#define ARCHIVIATORE_VOLUME_CLASS_H

#include <vector>
#include <ctime>
#include <cstdio>

// TODO: prima o poi arriviamo a far senza di questi define
#define NEL 15                // n0 elevazioni massimo
#define NUM_AZ_X_PPI 400

// TODO: per compatibilità con la libsp20, anche questo toglierlo in futuro
extern int elev_array[NEL];


namespace cumbac {

struct LoadLogEntry
{
    double theta;
    double alpha;

    LoadLogEntry(double theta, double alpha)
        : theta(theta), alpha(alpha)
    {
    }
};

struct LoadLog
{
    std::vector<LoadLogEntry> entries;

    void log(double theta, double alpha)
    {
        entries.push_back(LoadLogEntry(theta, alpha));
    }

    void print(FILE* out);
};

struct Ray
{
    std::vector<unsigned char> ray;
    short alfa_true, teta_true;
    short teta, alfa;

    Ray();
};

struct PolarScan : public std::vector<Ray>
{
    PolarScan();
};

struct VolumeStats
{
    unsigned count_zeros[NEL];
    unsigned count_ones[NEL];
    unsigned count_others[NEL];
    unsigned sum_others[NEL];

    void print(FILE* out);
};

class Volume
{
public:
    time_t acq_date;
    double size_cell;
    bool declutter_rsp; // ?

    //dato di base volume polare, struttura definita in libSP20
    //PolarScan vol_pol[NEL];
    Ray vol_pol[NEL][NUM_AZ_X_PPI];

    // Log of what has been loaded on each beam
    LoadLog load_log[NEL][NUM_AZ_X_PPI];

    //numero raggi per elevazione
    int nbeam_elev[NEL];

    Volume();

    void read_sp20(const char* nome_file);
    void read_odim(const char* nome_file);

    void compute_stats(VolumeStats& stats) const;

protected:
    void fill_beam(double theta, double alpha, unsigned size, const unsigned char* data);
    void merge_beam(int el_num, int az_num, double theta, double alpha, unsigned size, const unsigned char* dati);
};


}

#endif
