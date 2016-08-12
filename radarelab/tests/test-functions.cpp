#include "radarelab/utils/tests.h"
#include "radarelab/algo/dbz.h"
#include "radarelab/volume.h"
#include "radarelab/logging.h"

using namespace radarelab::utils::tests;
using namespace radarelab;

namespace {

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("functions");

void Tests::register_tests() {

add_method("beam_blocking_correction", []() {
    // Test BeamBlockingCorrection
    wassert(actual((unsigned)DBtoBYTE(algo::DBZ::beam_blocking_correction(BYTEtoDB(128),50))) == 138u);
});

}
}
