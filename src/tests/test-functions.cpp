#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"

using namespace wibble::tests;

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

}
