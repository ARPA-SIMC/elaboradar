#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#ifdef __cplusplus
extern "C" {
#endif
#include <func_Z_R.h>
#ifdef __cplusplus
}
#endif

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
    wassert(actual(NormalizzoData(0)) == 0); // TODO: check value
    wassert(actual(NormalizzoData(1)) == 0); // TODO: check value
}

template<> template<>
void to::test<2>()
{
    // Test BeamBlockingCorrection
    CUM_BAC* cb = new CUM_BAC("SPC");
    wassert(actual(DBtoBYTE(cb->BeamBlockingCorrection(128,50))) == 138);
}
}
