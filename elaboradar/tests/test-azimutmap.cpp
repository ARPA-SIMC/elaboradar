#include <wibble/tests.h>
#include <cmath>
#include <elaboradar/azimuthmap.h>
#include <elaboradar/logging.h>

using namespace std;
using namespace wibble::tests;
using namespace elaboradar;

namespace tut {

struct azimuthmap_shar {
    Logging logging;
};
TESTGRP(azimuthmap);

template<> template<>
void to::test<1>()
{
    typedef azimuthmap::Position Pos;
    UniformAzimuthMap am(360);

    wassert(actual(am.closest(0.0)) == Pos(0, 0));
    wassert(actual(am.closest(0.4)) == Pos(0, 0));
    wassert(actual(am.closest(0.5)) == Pos(1, 1));
    wassert(actual(am.closest(0.6)) == Pos(1, 1));
    wassert(actual(am.closest(359.0)) == Pos(359, 359));
    wassert(actual(am.closest(359.4)) == Pos(359, 359));
    wassert(actual(am.closest(359.5)) == Pos(0, 0));
    wassert(actual(am.closest(360.0)) == Pos(0, 0));
}

template<> template<>
void to::test<2>()
{
    typedef azimuthmap::Position Pos;
    UniformAzimuthMap am(360);

    vector<Pos> res = am.intersecting(0, 2);
    wassert(actual(res.size()) == 3);
    wassert(actual(res[0]) == Pos(359, 359));
    wassert(actual(res[1]) == Pos(0, 0));
    wassert(actual(res[2]) == Pos(1, 1));

    res = am.intersecting(180, 2.1);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == Pos(178, 178));
    wassert(actual(res[1]) == Pos(179, 179));
    wassert(actual(res[2]) == Pos(180, 180));
    wassert(actual(res[3]) == Pos(181, 181));
    wassert(actual(res[4]) == Pos(182, 182));

    res = am.intersecting(360, 2.1);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == Pos(358, 358));
    wassert(actual(res[1]) == Pos(359, 359));
    wassert(actual(res[2]) == Pos(0, 0));
    wassert(actual(res[3]) == Pos(1, 1));
    wassert(actual(res[4]) == Pos(2, 2));
}

template<> template<>
void to::test<3>()
{
    typedef azimuthmap::Position Pos;
    NonuniformAzimuthMap am;
    for (unsigned i = 0; i < 360; ++i)
        am.add(i, i);

    wassert(actual(am.closest(0.0)) == Pos(0, 0));
    wassert(actual(am.closest(0.4)) == Pos(0, 0));
    wassert(actual(am.closest(0.5)) == Pos(1, 1));
    wassert(actual(am.closest(0.6)) == Pos(1, 1));
    wassert(actual(am.closest(359.0)) == Pos(359, 359));
    wassert(actual(am.closest(359.4)) == Pos(359, 359));
    wassert(actual(am.closest(359.5)) == Pos(0, 0));
    wassert(actual(am.closest(360.0)) == Pos(0, 0));
}

template<> template<>
void to::test<4>()
{
    typedef azimuthmap::Position Pos;
    NonuniformAzimuthMap am;
    for (unsigned i = 0; i < 360; ++i)
        am.add(i, i);

    vector<Pos> res = am.intersecting(0, 2);
    wassert(actual(res.size()) == 3);
    wassert(actual(res[0]) == Pos(359, 359));
    wassert(actual(res[1]) == Pos(0, 0));
    wassert(actual(res[2]) == Pos(1, 1));

    res = am.intersecting(180, 2.1);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == Pos(178, 178));
    wassert(actual(res[1]) == Pos(179, 179));
    wassert(actual(res[2]) == Pos(180, 180));
    wassert(actual(res[3]) == Pos(181, 181));
    wassert(actual(res[4]) == Pos(182, 182));

    res = am.intersecting(360, 2.1);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == Pos(358, 358));
    wassert(actual(res[1]) == Pos(359, 359));
    wassert(actual(res[2]) == Pos(0, 0));
    wassert(actual(res[3]) == Pos(1, 1));
    wassert(actual(res[4]) == Pos(2, 2));
}

}

