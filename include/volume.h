#ifndef ARCHIVIATORE_VOLUME_CLASS_H
#define ARCHIVIATORE_VOLUME_CLASS_H

#ifdef __cplusplus
extern "C" {
#endif
// libreria radar
#include <func_SP20read.h>
#ifdef __cplusplus
}
#endif

// TODO: prima o poi arriviamo a far senza di questi define
#define NEL 15                // n0 elevazioni massimo
#define NUM_AZ_X_PPI 400

// TODO: per compatibilità con la libsp20, anche questo toglierlo in futuro
extern int elev_array[NEL];


namespace cumbac {

class Volume
{
public:
    time_t acq_date;
    double size_cell;
    bool declutter_rsp; // ?

    //dato di base volume polare, struttura definita in libSP20
    struct VOL_POL vol_pol[NEL][NUM_AZ_X_PPI];

    //numero raggi per elevazione
    int nbeam_elev[NEL];

    Volume();

    void fill_beam(double theta, double alpha, unsigned size, const unsigned char* data);
    void merge_beam(VOL_POL* raggio, double theta, double alpha, int az_num, int el_num, unsigned size, const unsigned char* dati);
};


}

#endif
