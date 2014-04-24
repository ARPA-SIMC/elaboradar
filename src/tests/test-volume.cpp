#include <wibble/tests.h>
#include <cmath>
#include "volume.h"
#include "site.h"
#include "logging.h"

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct volume_shar {
    Logging logging;
};
TESTGRP(volume);

template<> template<>
void to::test<1>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    Volume<double> vsp20;
    vsp20.read_sp20("testdata/DBP2_070120141530_GATTATICO", VolumeLoadOptions(Site::get("GAT"), false, false));

    wassert(actual(vsp20.scan(0).get_elevation_rad(0)) == vsp20.scan(0).get_elevation(0) * M_PI / 180);
}

}
