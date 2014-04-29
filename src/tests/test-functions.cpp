#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct functions_shar {
    Logging logging;
};
TESTGRP(functions);

template<> template<>
void to::test<1>()
{
    // Test NormalizzoData
    CUM_BAC* cb = new CUM_BAC("SPC");
    wassert(actual(cb->NormalizzoData(0)) == 0); // TODO: check value
    wassert(actual(cb->NormalizzoData(1)) == 0); // TODO: check value
}

template<> template<>
void to::test<2>()
{
    // Test BeamBlockingCorrection
    CUM_BAC* cb = new CUM_BAC("SPC");
    wassert(actual((unsigned)DBtoBYTE(cb->BeamBlockingCorrection(128,50))) == 138);
}
}
