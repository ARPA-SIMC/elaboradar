#include <wibble/tests.h>
#include <cmath>
#include <elaboradar/algo/azimuth_resample.h>
#include <elaboradar/logging.h>

using namespace std;
using namespace wibble::tests;
using namespace elaboradar;
using namespace elaboradar::algo::azimuthresample;

namespace elaboradar {
namespace algo {
namespace azimuthresample {

ostream& operator<<(ostream& out, const pair<double, unsigned>& p)
{
    return out << "(" << p.first << "," << p.second << ")";
}

}
}
}

namespace tut {

struct azimuthmap_shar {
    Logging logging;
};
TESTGRP(azimuthmap);

template<> template<>
void to::test<1>()
{
    Eigen::VectorXd azimuths(360);
    for (unsigned i = 0; i < 360; ++i)
        azimuths(i) = (double)i;
    AzimuthIndex am(azimuths);

    wassert(actual(am.closest(0.0)) == make_pair(0.0, 0u));
    wassert(actual(am.closest(0.4)) == make_pair(0.0, 0u));
    wassert(actual(am.closest(0.5)) == make_pair(1.0, 1u));
    wassert(actual(am.closest(0.6)) == make_pair(1.0, 1u));
    wassert(actual(am.closest(359.0)) == make_pair(359.0, 359u));
    wassert(actual(am.closest(359.4)) == make_pair(359.0, 359u));
    wassert(actual(am.closest(359.5)) == make_pair(0.0, 0u));
    wassert(actual(am.closest(360.0)) == make_pair(0.0, 0u));
}

template<> template<>
void to::test<4>()
{
    Eigen::VectorXd azimuths(360);
    for (unsigned i = 0; i < 360; ++i)
        azimuths(i) = (double)i;
    AzimuthIndex am(azimuths);

    // From -177.5 to 182.5
    auto res = am.intersecting(0, 4);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == make_pair(358.0, 358u));
    wassert(actual(res[1]) == make_pair(359.0, 359u));
    wassert(actual(res[2]) == make_pair(0.0, 0u));
    wassert(actual(res[3]) == make_pair(1.0, 1u));
    wassert(actual(res[4]) == make_pair(2.0, 2u));

    // From -177.1 to 182.9
    res = am.intersecting(180, 4.8);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == make_pair(178.0, 178u));
    wassert(actual(res[1]) == make_pair(179.0, 179u));
    wassert(actual(res[2]) == make_pair(180.0, 180u));
    wassert(actual(res[3]) == make_pair(181.0, 181u));
    wassert(actual(res[4]) == make_pair(182.0, 182u));

    // From -176.6 to 183.1
    res = am.intersecting(360, 5.2);
    wassert(actual(res.size()) == 7);
    wassert(actual(res[0]) == make_pair(357.0, 357u));
    wassert(actual(res[1]) == make_pair(358.0, 358u));
    wassert(actual(res[2]) == make_pair(359.0, 359u));
    wassert(actual(res[3]) == make_pair(0.0, 0u));
    wassert(actual(res[4]) == make_pair(1.0, 1u));
    wassert(actual(res[5]) == make_pair(2.0, 2u));
    wassert(actual(res[6]) == make_pair(3.0, 3u));
}

}

