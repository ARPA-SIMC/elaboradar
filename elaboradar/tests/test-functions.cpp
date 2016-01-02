#include "elaboradar/utils/tests.h"
#include "elaboradar/algo/utils.h"
#include "elaboradar/volume.h"
#include "elaboradar/logging.h"

using namespace elaboradar::utils::tests;
using namespace elaboradar;

namespace {

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("functions");

void Tests::register_tests() {

add_method("beam_blocking_correction", []() {
    // Test BeamBlockingCorrection
    wassert(actual((unsigned)DBtoBYTE(algo::beam_blocking_correction(BYTEtoDB(128),50))) == 138u);
});

}
}
