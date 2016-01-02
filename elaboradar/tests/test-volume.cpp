#include "elaboradar/utils/tests.h"
#include "elaboradar/volume.h"
#include "elaboradar/sp20.h"
#include "elaboradar/logging.h"
#include <cmath>

using namespace elaboradar::utils::tests;
using namespace elaboradar;

namespace {

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("volume");

void Tests::register_tests() {

add_method("load_sp20", []() {
    // Test loading of a radar volume via SP20
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    volume::SP20Loader sp20;
    volume::Scans<double> vsp20;
    sp20.vol_z = &vsp20;
    sp20.load("../testdata/DBP2_070120141530_GATTATICO");

    //wassert(actual(liSP20.scan(0).get_elevation_rad(0)) == liSP20.scan(0).get_elevation(0) * M_PI / 180);
});

}
}
