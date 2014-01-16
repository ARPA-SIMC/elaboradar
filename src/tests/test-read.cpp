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
    wassert(actual(res).istrue());
    delete cb;
}

}
