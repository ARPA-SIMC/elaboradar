#include <wibble/tests.h>
#include <elaboradar/algo/utils.h>
#include <elaboradar/volume.h>
#include <elaboradar/logging.h>

using namespace wibble::tests;
using namespace elaboradar;

namespace tut {

struct functions_shar {
    Logging logging;
};
TESTGRP(functions);

template<> template<>
void to::test<1>()
{
#if 0
    // Test NormalizzoData
    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "SPC");
    wassert(actual(cb->NormalizzoData(0)) == 0); // TODO: check value
    wassert(actual(cb->NormalizzoData(1)) == 0); // TODO: check value
#endif
}

template<> template<>
void to::test<2>()
{
    // Test BeamBlockingCorrection
    wassert(actual((unsigned)DBtoBYTE(algo::beam_blocking_correction(BYTEtoDB(128),50))) == 138);
}
}
