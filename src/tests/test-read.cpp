#include <wibble/tests.h>
#include "cum_bac.h"

using namespace wibble::tests;

namespace tut {

struct read_sp20_shar {
};
TESTGRP(read_sp20);

template<> template<>
void to::test<1>()
{
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    CUM_BAC* cb = new CUM_BAC;
    bool res = cb->read_sp20_volume(fname, "GAT", 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // Ensure that nbeam_elev has been filled with the right values
    wassert(actual(cb->nbeam_elev[0]) == 400);
    wassert(actual(cb->nbeam_elev[1]) == 400);
    wassert(actual(cb->nbeam_elev[2]) == 400);
    wassert(actual(cb->nbeam_elev[3]) == 400);
    wassert(actual(cb->nbeam_elev[4]) == 400);
    // TODO: check some of the vol_pol contents
    delete cb;
}

}
