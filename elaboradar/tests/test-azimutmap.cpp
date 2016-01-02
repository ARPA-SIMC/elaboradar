#include "elaboradar/utils/tests.h"
#include "elaboradar/algo/azimuth_resample.h"
#include "elaboradar/logging.h"
#include <cmath>

using namespace std;
using namespace elaboradar::utils::tests;
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

namespace {

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("azimuthmap");

void Tests::register_tests() {

add_method("closest", []() {
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
});

add_method("intersecting", []() {
    Eigen::VectorXd azimuths(360);
    for (unsigned i = 0; i < 360; ++i)
        azimuths(i) = (double)i;
    AzimuthIndex am(azimuths);

    // From -177.5 to 182.5
    auto res = am.intersecting(0, 4, 1);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == make_pair(358.0, 358u));
    wassert(actual(res[1]) == make_pair(359.0, 359u));
    wassert(actual(res[2]) == make_pair(0.0, 0u));
    wassert(actual(res[3]) == make_pair(1.0, 1u));
    wassert(actual(res[4]) == make_pair(2.0, 2u));

    // From -177.1 to 182.9
    res = am.intersecting(180, 4.8, 1);
    wassert(actual(res.size()) == 5);
    wassert(actual(res[0]) == make_pair(178.0, 178u));
    wassert(actual(res[1]) == make_pair(179.0, 179u));
    wassert(actual(res[2]) == make_pair(180.0, 180u));
    wassert(actual(res[3]) == make_pair(181.0, 181u));
    wassert(actual(res[4]) == make_pair(182.0, 182u));

    // From -176.6 to 183.1
    res = am.intersecting(360, 5.2, 1);
    wassert(actual(res.size()) == 7);
    wassert(actual(res[0]) == make_pair(357.0, 357u));
    wassert(actual(res[1]) == make_pair(358.0, 358u));
    wassert(actual(res[2]) == make_pair(359.0, 359u));
    wassert(actual(res[3]) == make_pair(0.0, 0u));
    wassert(actual(res[4]) == make_pair(1.0, 1u));
    wassert(actual(res[5]) == make_pair(2.0, 2u));
    wassert(actual(res[6]) == make_pair(3.0, 3u));
});

// Test closest resample from a bigger volume to a smaller
add_method("reduce", []() {
    // Create a volume with one elevation in which each beam contains as
    // samples its azimuth in degrees
    Volume<double> vol(360);
    auto& scan = vol.make_scan(0, 10, 0, 250);
    for (unsigned i = 0; i < 360; ++i)
    {
        scan.azimuths_real(i) = i;
        for (unsigned j = 0; j < 10; ++j)
            scan.set(i, j, i);
    }

    Volume<double> dst(180);
    algo::azimuthresample::Closest<double> resample;
    resample.resample_volume(vol, dst, 1);

    wassert(actual(dst.scan(0).get(179, 0)) == 358);
    wassert(actual(dst.scan(0).get(0, 0)) == 0);
    wassert(actual(dst.scan(0).get(1, 0)) == 2);
});

// Test closest resample from a smaller volume to a bigger
add_method("enlarge", []() {
    // Create a volume with one elevation in which each beam contains as
    // samples its azimuth in degrees
    Volume<double> vol(180);
    auto& scan = vol.make_scan(0, 10, 0, 250);
    for (unsigned i = 0; i < 180; ++i)
    {
        scan.azimuths_real(i) = i * 2;
        for (unsigned j = 0; j < 10; ++j)
            scan.set(i, j, i);
    }

    Volume<double> dst(360);
    algo::azimuthresample::Closest<double> resample;
    resample.resample_volume(vol, dst, 2);

    wassert(actual(dst.scan(0).get(356, 0)) == 178);
    wassert(actual(dst.scan(0).get(357, 0)) == 179);
    wassert(actual(dst.scan(0).get(358, 0)) == 179);
    wassert(actual(dst.scan(0).get(359, 0)) == 0);
    wassert(actual(dst.scan(0).get(0, 0)) == 0);
    wassert(actual(dst.scan(0).get(1, 0)) == 1);
    wassert(actual(dst.scan(0).get(2, 0)) == 1);
    wassert(actual(dst.scan(0).get(3, 0)) == 2);
    wassert(actual(dst.scan(0).get(4, 0)) == 2);
});

// Test closest resample from a smaller volume to a bigger, with gaps
add_method("enlarge_with_gaps", []() {
    // Create a volume with one elevation in which each beam contains as
    // samples its azimuth in degrees
    Volume<double> vol(90);
    auto& scan = vol.make_scan(0, 10, 0, 250);
    for (unsigned i = 0; i < 90; ++i)
    {
        scan.azimuths_real(i) = i * 4;
        for (unsigned j = 0; j < 10; ++j)
            scan.set(i, j, i);
    }

    Volume<double> dst(360);
    algo::azimuthresample::Closest<double> resample;
    resample.resample_volume(vol, dst, 2);

    double missing = BYTEtoDB(1);

    wassert(actual(dst.scan(0).get(356, 0)) == 89);
    wassert(actual(dst.scan(0).get(357, 0)) == 89);
    wassert(actual(dst.scan(0).get(358, 0)) == missing);
    wassert(actual(dst.scan(0).get(359, 0)) == 0);
    wassert(actual(dst.scan(0).get(0, 0)) == 0);
    wassert(actual(dst.scan(0).get(1, 0)) == 0);
    wassert(actual(dst.scan(0).get(2, 0)) == missing);
    wassert(actual(dst.scan(0).get(3, 0)) == 1);
    wassert(actual(dst.scan(0).get(4, 0)) == 1);
});

// Test max-of-closest resample from a bigger volume to a smaller
add_method("max-of-closest-reduce", []() {
    // Create a volume with one elevation in which each beam contains as
    // samples its azimuth in degrees
    Volume<double> vol(360);
    auto& scan = vol.make_scan(0, 10, 0, 250);
    for (unsigned i = 0; i < 360; ++i)
    {
        scan.azimuths_real(i) = i;
        for (unsigned j = 0; j < 10; ++j)
            scan.set(i, j, i);
    }

    Volume<double> dst(180);
    algo::azimuthresample::MaxOfClosest<double> resample;
    resample.resample_volume(vol, dst, 1);

    wassert(actual(dst.scan(0).get(179, 0)) == 359);
    wassert(actual(dst.scan(0).get(0, 0)) == 359);
    wassert(actual(dst.scan(0).get(1, 0)) == 3);
});

// Test max-of-closest resample from a smaller volume to a bigger
add_method("max-of-closest-enlarge", []() {
    // Create a volume with one elevation in which each beam contains as
    // samples its azimuth in degrees
    Volume<double> vol(180);
    auto& scan = vol.make_scan(0, 10, 0, 250);
    for (unsigned i = 0; i < 180; ++i)
    {
        scan.azimuths_real(i) = i * 2;
        for (unsigned j = 0; j < 10; ++j)
            scan.set(i, j, i);
    }

    Volume<double> dst(360);
    algo::azimuthresample::MaxOfClosest<double> resample;
    resample.resample_volume(vol, dst, 2);

    wassert(actual(dst.scan(0).get(356, 0)) == 178);
    wassert(actual(dst.scan(0).get(357, 0)) == 179);
    wassert(actual(dst.scan(0).get(358, 0)) == 179);
    wassert(actual(dst.scan(0).get(359, 0)) == 179);
    wassert(actual(dst.scan(0).get(0, 0)) == 0);
    wassert(actual(dst.scan(0).get(1, 0)) == 1);
    wassert(actual(dst.scan(0).get(2, 0)) == 1);
    wassert(actual(dst.scan(0).get(3, 0)) == 2);
    wassert(actual(dst.scan(0).get(4, 0)) == 2);
});

// Test closest resample from a smaller volume to a bigger, with gaps
add_method("closest-enlarge", []() {
    // Create a volume with one elevation in which each beam contains as
    // samples its azimuth in degrees
    Volume<double> vol(90);
    auto& scan = vol.make_scan(0, 10, 0, 250);
    for (unsigned i = 0; i < 90; ++i)
    {
        scan.azimuths_real(i) = i * 4;
        for (unsigned j = 0; j < 10; ++j)
            scan.set(i, j, i);
    }

    Volume<double> dst(360);
    algo::azimuthresample::MaxOfClosest<double> resample;
    resample.resample_volume(vol, dst, 2);

    double missing = BYTEtoDB(1);

    wassert(actual(dst.scan(0).get(356, 0)) == 89);
    wassert(actual(dst.scan(0).get(357, 0)) == 89);
    wassert(actual(dst.scan(0).get(358, 0)) == missing);
    wassert(actual(dst.scan(0).get(359, 0)) == 0);
    wassert(actual(dst.scan(0).get(0, 0)) == 0);
    wassert(actual(dst.scan(0).get(1, 0)) == 0);
    wassert(actual(dst.scan(0).get(2, 0)) == missing);
    wassert(actual(dst.scan(0).get(3, 0)) == 1);
    wassert(actual(dst.scan(0).get(4, 0)) == 1);
});

}

}
