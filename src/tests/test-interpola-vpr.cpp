#include <wibble/tests.h>
#include "interpola_vpr.h"
#include "cum_bac.h"
#include "logging.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "setwork.h"

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct intvpr_shar {
    Logging logging;
};
TESTGRP(intvpr);

namespace {
// hvprmax = 900, livmin = 0
// wassert(actual(iv.B) == 0.768);
// wassert(actual(iv.E) == 1.019);
// wassert(actual(iv.G) == 0.246);
// wassert(actual(iv.C) == 0.63);
// wassert(actual(iv.F) == -0.12);
const float vpr0[] = { 0.337, 0.519, 0.578, 0.788, 1.196, 1.192, 0.678, 0.478,
    0.465, 0.461, 0.399, 0.377, 0.330, 0.313, 0.283, 0.261, 0.253, 0.247,
    0.241, 0.235 };
}
// 1300 hvprmax 300 livmin
const float vpr1[] = { -9999.000000, 2.136112, 1.879434, 1.757904, 2.025282, 3.777671, 4.346189, 2.004490, 0.988849, 0.743071, 0.659562, 0.571394, 0.537457, 0.565114, 0.593373, 0.575278, 0.539934, 0.511525, 0.501043, 0.489665, 0.480079, 0.474079, 0.468079, 0.462079, 0.456079, 0.450079, 0.444079, 0.438079, 0.432079, 0.426079, 0.420079, 0.414079, 0.408079, 0.402079, 0.396079, 0.390079, 0.384079, 0.378079, 0.372079, 0.366079, 0.360079, 0.354079, 0.348079, 0.342079, 0.336079, 0.330079, 0.324079, 0.318079, 0.312079, 0.306079, 0.300079, 0.294079, 0.288079, 0.282079, 0.276079, 0.270079, 0.264079, 0.258079, 0.252079, 0.246079, 0.240079, 0.234079, 0.228079, 0.222079, 0.216079, 0.210079, 0.204079, 0.198079, 0.192079, 0.186079 };

template<> template<>
void to::test<1>()
{
    InterpolaVPR_NR iv;
    int res = iv.interpola_VPR(vpr1, 1300, 300);
    wassert(actual(res) == 0);
    wassert(actual(round(iv.B * 10)) ==  33);
    wassert(actual(round(iv.E * 10)) ==  12);
    wassert(actual(round(iv.G * 10)) ==   2);
    wassert(actual(round(iv.C * 10)) ==  23);
    wassert(actual(round(iv.F * 10)) ==  -8);
    //wassert(actual(iv.chisqfin) == 0);
    //wassert(actual(iv.rmsefin) == 0);
    //double vpr_int[NMAXLAYER];
}

template<> template<>
void to::test<2>()
{
    InterpolaVPR_GSL iv;
    int res = iv.interpola_VPR(vpr1, 1300, 300);
    fprintf(stderr, "%f %f %f %f %f\n", iv.B, iv.E, iv.G, iv.C, iv.F);
    wassert(actual(res) == 0);
    wassert(actual(round(iv.B * 10)) ==  33);
    wassert(actual(round(iv.E * 10)) ==  12);
    wassert(actual(round(iv.G * 10)) ==   2);
    wassert(actual(round(iv.C * 10)) ==  23);
    wassert(actual(round(iv.F * 10)) ==  -8);
}

}

