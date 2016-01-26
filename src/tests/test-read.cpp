#include "radarelab/utils/tests.h"
#include <radarelab/sp20.h>
#include <radarelab/odim.h>
#include <radarelab/logging.h>
#include <radarelab/algo/azimuth_resample.h>
#include "cum_bac.h"
#include "config.h"
#include "site.h"
#include <radarlib/radar.hpp>
#include <stdio.h>
#include <vector>

using namespace radarelab::utils::tests;
using namespace radarelab;
using namespace elaboradar;
using namespace std;

namespace {

void test_0120141530gat_SP20(const Volume<double>& v)
{
    // Ensure that the beam sizes are what we expect
    wassert(actual(v.scan(0).beam_size) == 494);
    wassert(actual(v.scan(1).beam_size) == 494);
    wassert(actual(v.scan(2).beam_size) == 494);
    wassert(actual(v.scan(3).beam_size) == 494);
    wassert(actual(v.scan(4).beam_size) == 494);
    wassert(actual(v.scan(5).beam_size) == 494);

    // Ensure that the beam azimuth are what we expect
    /*
    wassert(actual(v.scan(0)[0].alfa) == 0);
    wassert(actual(v.scan(0)[1].alfa) == 10);
    wassert(actual(v.scan(1)[1].alfa) == 10);
    wassert(actual(v.scan(2)[1].alfa) == 10);
    wassert(actual(v.scan(3)[1].alfa) == 10);
    wassert(actual(v.scan(4)[1].alfa) == 10);
    wassert(actual(v.scan(5)[1].alfa) == 10);
    */

    // Check other header fields
    wassert(actual(v.load_info->acq_date) == 1389108600);
    wassert(actual(v.scan(0).cell_size) == 250);

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats;
    v.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 145426);
    wassert(actual(stats.count_zeros[1]) == 184085);
    wassert(actual(stats.count_zeros[2]) == 193145);
    wassert(actual(stats.count_zeros[3]) == 196213);
    wassert(actual(stats.count_zeros[4]) == 196085);
    wassert(actual(stats.count_zeros[5]) == 196080);
    wassert(actual(stats.count_ones[0]) == 104);
    wassert(actual(stats.count_ones[1]) == 204);
    wassert(actual(stats.count_ones[2]) == 104);
    wassert(actual(stats.count_ones[3]) ==  58);
    wassert(actual(stats.count_ones[4]) ==  61);
    wassert(actual(stats.count_ones[5]) ==  38);
    wassert(actual(stats.count_others[0]) == 52070);
    wassert(actual(stats.count_others[1]) == 13311);
    wassert(actual(stats.count_others[2]) ==  4351);
    wassert(actual(stats.count_others[3]) ==  1329);
    wassert(actual(stats.count_others[4]) ==  1454);
    wassert(actual(stats.count_others[5]) ==  1482);
    wassert(actual(stats.sum_others[0]) == 4759195);
    wassert(actual(stats.sum_others[1]) ==  914997);
    wassert(actual(stats.sum_others[2]) ==  257459);
    wassert(actual(stats.sum_others[3]) ==   46349);
    wassert(actual(stats.sum_others[4]) ==   78749);
    wassert(actual(stats.sum_others[5]) ==   90563);
}

void test_0120141530gat_ODIM(const Volume<double>& v)
{
    // Ensure that the beam sizes are what we expect
    wassert(actual(v.scan(0).beam_size) == 494);
    wassert(actual(v.scan(1).beam_size) == 494);
    wassert(actual(v.scan(2).beam_size) == 494);
    wassert(actual(v.scan(3).beam_size) == 494);
    wassert(actual(v.scan(4).beam_size) == 494);
    wassert(actual(v.scan(5).beam_size) == 494);

    // Ensure that the beam azimuth are what we expect
    /*
    wassert(actual(v.scan(0)[0].alfa) == 0);
    wassert(actual(v.scan(0)[1].alfa) == 10);
    wassert(actual(v.scan(1)[1].alfa) == 10);
    wassert(actual(v.scan(2)[1].alfa) == 10);
    wassert(actual(v.scan(3)[1].alfa) == 10);
    wassert(actual(v.scan(4)[1].alfa) == 10);
    wassert(actual(v.scan(5)[1].alfa) == 10);
    */

    // Check other header fields
    wassert(actual(v.load_info->acq_date) == 1389108600);
    wassert(actual(v.scan(0).cell_size) == 250);

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats;
    v.compute_stats(stats);
    // stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 145420);
    wassert(actual(stats.count_zeros[1]) == 184078);
    wassert(actual(stats.count_zeros[2]) == 193143);
    wassert(actual(stats.count_zeros[3]) == 196220);
    wassert(actual(stats.count_zeros[4]) == 196090);
    wassert(actual(stats.count_zeros[5]) == 196080);
    wassert(actual(stats.count_ones[0]) == 106);
    wassert(actual(stats.count_ones[1]) == 203);
    wassert(actual(stats.count_ones[2]) == 105);
    wassert(actual(stats.count_ones[3]) ==  56);
    wassert(actual(stats.count_ones[4]) ==  58);
    wassert(actual(stats.count_ones[5]) ==  43);
    wassert(actual(stats.count_others[0]) == 52074);
    wassert(actual(stats.count_others[1]) == 13319);
    wassert(actual(stats.count_others[2]) ==  4352);
    wassert(actual(stats.count_others[3]) ==  1324);
    wassert(actual(stats.count_others[4]) ==  1452);
    wassert(actual(stats.count_others[5]) ==  1477);
    wassert(actual(stats.sum_others[0]) == 4740629);
    wassert(actual(stats.sum_others[1]) ==  911724);
    wassert(actual(stats.sum_others[2]) ==  256528);
    wassert(actual(stats.sum_others[3]) ==   46119);
    wassert(actual(stats.sum_others[4]) ==   78528);
    wassert(actual(stats.sum_others[5]) ==   90287);
}

struct Difference
{
    unsigned idx;
    double val;
    Difference() : idx(0), val(0) {}
    Difference(unsigned idx, double val)
        : idx(idx), val(val) {}
};

void test_volumes_equal(const Volume<double>& vsp20, const Volume<double>& vodim)
{
    using namespace std;

    wassert(actual(vsp20.size()) == vodim.size());

    unsigned failed_beams = 0;
    for (unsigned ie = 0; ie < vsp20.size(); ++ie)
    {
        RADARELAB_UTILS_TEST_INFO(testinfo);
        testinfo() << "elevation " << ie;

        for (unsigned ia = 0; ia < vsp20.scan(ie).beam_count; ++ia)
        {
            testinfo() << "elevation " << ie << " angle " << ia;

            //wassert(actual(vsp20.scan(ie)[ia].teta_true) == vodim.scan(ie)[ia].teta_true);
            //wassert(actual(vsp20.scan(ie)[ia].teta) == vodim.scan(ie)[ia].teta);
            //wassert(actual(vsp20.scan(ie)[ia].alfa) == vodim.scan(ie)[ia].alfa);
            wassert(actual(vsp20.scan(ie).beam_size) == vodim.scan(ie).beam_size);

            vector<Difference> vals_sp20;
            vector<Difference> vals_odim;
            for (unsigned ib = 0; ib < vsp20.scan(ie).beam_size; ++ib)
            {
                if (vsp20.scan(ie).get(ia, ib) != vodim.scan(ie).get(ia, ib))
                {
                    vals_sp20.push_back(Difference(ib, vsp20.scan(ie).get(ia, ib)));
                    vals_odim.push_back(Difference(ib, vodim.scan(ie).get(ia, ib)));
                }
            }
            if (!vals_sp20.empty())
            {
#if 0
                printf("sp20 vp[%u][%u]:", ie, ia);
                for (vector<Difference>::const_iterator i = vals_sp20.begin(); i != vals_sp20.end(); ++i)
                    printf(" %u:%d", i->idx, (int)i->val);
                printf("\n");
                printf("odim vp[%u][%u]:", ie, ia);
                for (vector<Difference>::const_iterator i = vals_odim.begin(); i != vals_odim.end(); ++i)
                    printf(" %u:%d", i->idx, (int)i->val);
                printf("\n");
                const volume::LoadLog& llsp20 = vsp20.scan(ie).load_info().get_beam_load_log(ia);
                const volume::LoadLog& llodim = vodim.scan(ie).load_info().get_beam_load_log(ia);
                if (llsp20 != llodim)
                {
                    printf("sp20 vp[%u][%u] load log: ", ie, ia); llsp20.print(stdout);
                    printf("odim vp[%u][%u] load log: ", ie, ia); llodim.print(stdout);
                }
#endif
                ++failed_beams;
            }
        }
    }

    wassert(actual(failed_beams) == 0);
}

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("read");

void Tests::register_tests() {

add_method("read_sp20", []() {
    using namespace radarelab::volume;
    // Test loading of a radar volume via SP20
    const Site& gat = Site::get("GAT");
    SP20Loader loader;
    Scans<double> ssp20;
    loader.vol_z = &ssp20;
    loader.load("../testdata/DBP2_070120141530_GATTATICO");

    wassert(actual(ssp20.size()) == 6);
    wassert(actual(round(ssp20.scan(0).elevation * 10000)) ==  5274);
    wassert(actual(round(ssp20.scan(1).elevation * 10000)) == 14064);
    wassert(actual(round(ssp20.scan(2).elevation * 10000)) == 22854);
    wassert(actual(round(ssp20.scan(3).elevation * 10000)) == 31644);
    wassert(actual(round(ssp20.scan(4).elevation * 10000)) == 41313);
    wassert(actual(round(ssp20.scan(5).elevation * 10000)) == 50103);

    Volume<double> vsp20;
    algo::azimuthresample::MaxOfClosest<double> resampler;
    resampler.resample_volume(ssp20, vsp20, 1);
    // Check the contents of what we read
    wassert(test_0120141530gat_SP20(vsp20));
});

add_method("read_odim", []() {
    // Test loading of a radar volume via SP20
    static const char* fname = "../testdata/MSG1400715300U.101.h5";
    const Site& site(Site::get("GAT"));
    Volume<double> volume;
    CUM_BAC::read_odim_volume(volume, site, fname, false);
    // Check the contents of what we read
    wassert(test_0120141530gat_ODIM(volume));
});

add_method("read_sp20_1", []() {
    // Test loading of a radar volume via SP20
    static const char* fname = "../testdata/DBP2_060220140140_GATTATICO";

    const Site& site(Site::get("GAT"));
    Volume<double> volume;
    CUM_BAC::read_sp20_volume(volume, site, fname, 0, false);
    // TODO: Check the contents of what we read
    //wruntest(test_0120141530gat, volume);
});

add_method("read_sp20_2", []() {
    // Test loading of a radar volume via SP20
    static const char* fname = "../testdata/DBP2_020520141110_BOLOGNA";
    const Site& site(Site::get("SPC"));
    Volume<double> volume;
    CUM_BAC::read_sp20_volume(volume, site, fname, 0, false);
    // Make sure the elevations are what we want
    wassert(actual(volume.size()) == 11);
    wassert(actual(round(volume.scan(0).elevation * 10)) ==    5);
    wassert(actual(round(volume.scan(1).elevation * 10)) ==   14);
    wassert(actual(round(volume.scan(2).elevation * 10)) ==   23);
    wassert(actual(round(volume.scan(3).elevation * 10)) ==   32);
    wassert(actual(round(volume.scan(4).elevation * 10)) ==   41);
    wassert(actual(round(volume.scan(5).elevation * 10)) ==   50);
    wassert(actual(round(volume.scan(6).elevation * 10)) ==   70);
    wassert(actual(round(volume.scan(7).elevation * 10)) ==   95);
    wassert(actual(round(volume.scan(8).elevation * 10)) ==  130);
    wassert(actual(round(volume.scan(9).elevation * 10)) ==  180);
    wassert(actual(round(volume.scan(10).elevation * 10)) == 250);
    // TODO: Check the contents of what we read
    //wruntest(test_0120141530gat, cb->volume);
});

}

}
