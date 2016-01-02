#include <wibble/tests.h>
#include <cmath>
#include <elaboradar/volume.h>
#include <elaboradar/sp20.h>
#include <elaboradar/logging.h>

using namespace wibble::tests;
using namespace elaboradar;

namespace tut {

struct volume_shar {
    Logging logging;
};
TESTGRP(volume);

template<> template<>
void to::test<1>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    volume::SP20Loader sp20;
    volume::Scans<double> vsp20;
    sp20.vol_z = &vsp20;
    sp20.load("../testdata/DBP2_070120141530_GATTATICO");

    //wassert(actual(liSP20.scan(0).get_elevation_rad(0)) == liSP20.scan(0).get_elevation(0) * M_PI / 180);
}

}
