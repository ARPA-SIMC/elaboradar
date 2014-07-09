#include <wibble/tests.h>
#include "cum_bac.h"
#include "config.h"
#include "volume/sp20.h"
#include "volume/odim.h"
#include "volume/resample.h"
#include "site.h"
#include "logging.h"
#include <stdio.h>
#include <vector>

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct read_shar {
    Logging logging;
};
TESTGRP(read);

namespace {

void test_0120141530gat_SP20(WIBBLE_TEST_LOCPRM, const volume::LoadInfo& li, const Volume<double>& v)
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
    wassert(actual(li.acq_date) == 1389108600);
    wassert(actual(v.scan(0).cell_size) == 250);

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats;
    v.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 139396);
    wassert(actual(stats.count_zeros[1]) == 181135);
    wassert(actual(stats.count_zeros[2]) == 192056);
    wassert(actual(stats.count_zeros[3]) == 196000);
    wassert(actual(stats.count_zeros[4]) == 195715);
    wassert(actual(stats.count_zeros[5]) == 195897);
    wassert(actual(stats.count_ones[0]) ==  97);
    wassert(actual(stats.count_ones[1]) == 217);
    wassert(actual(stats.count_ones[2]) == 119);
    wassert(actual(stats.count_ones[3]) ==  70);
    wassert(actual(stats.count_ones[4]) ==  66);
    wassert(actual(stats.count_ones[5]) ==  43);
    wassert(actual(stats.count_others[0]) == 58107);
    wassert(actual(stats.count_others[1]) == 16248);
    wassert(actual(stats.count_others[2]) ==  5425);
    wassert(actual(stats.count_others[3]) ==  1530);
    wassert(actual(stats.count_others[4]) ==  1819);
    wassert(actual(stats.count_others[5]) ==  1660);
    wassert(actual(stats.sum_others[0]) == 5475665);
    wassert(actual(stats.sum_others[1]) == 1130268);
    wassert(actual(stats.sum_others[2]) ==  332252);
    wassert(actual(stats.sum_others[3]) ==   52204);
    wassert(actual(stats.sum_others[4]) ==  100286);
    wassert(actual(stats.sum_others[5]) ==  101241);
}

void test_0120141530gat_ODIM(WIBBLE_TEST_LOCPRM, const volume::LoadInfo& li, const Volume<double>& v)
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
    wassert(actual(li.acq_date) == 1389108600);
    wassert(actual(v.scan(0).cell_size) == 250);

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats;
    v.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 139391);
    wassert(actual(stats.count_zeros[1]) == 181129);
    wassert(actual(stats.count_zeros[2]) == 192055);
    wassert(actual(stats.count_zeros[3]) == 196011);
    wassert(actual(stats.count_zeros[4]) == 195722);
    wassert(actual(stats.count_zeros[5]) == 195897);
    wassert(actual(stats.count_ones[0]) ==  98);
    wassert(actual(stats.count_ones[1]) == 217);
    wassert(actual(stats.count_ones[2]) == 120);
    wassert(actual(stats.count_ones[3]) ==  67);
    wassert(actual(stats.count_ones[4]) ==  62);
    wassert(actual(stats.count_ones[5]) ==  49);
    wassert(actual(stats.count_others[0]) == 58111);
    wassert(actual(stats.count_others[1]) == 16254);
    wassert(actual(stats.count_others[2]) ==  5425);
    wassert(actual(stats.count_others[3]) ==  1522);
    wassert(actual(stats.count_others[4]) ==  1816);
    wassert(actual(stats.count_others[5]) ==  1654);
    wassert(actual(stats.sum_others[0]) == 5454271);
    wassert(actual(stats.sum_others[1]) == 1125972);
    wassert(actual(stats.sum_others[2]) ==  330926);
    wassert(actual(stats.sum_others[3]) ==   51911);
    wassert(actual(stats.sum_others[4]) ==   99965);
    wassert(actual(stats.sum_others[5]) ==  100920);
}

namespace {
struct Difference
{
    unsigned idx;
    double val;
    Difference() : idx(0), val(0) {}
    Difference(unsigned idx, double val)
        : idx(idx), val(val) {}
};
}

void test_loadinfo_equal(WIBBLE_TEST_LOCPRM, const volume::LoadInfo& vsp20, const volume::LoadInfo& vodim)
{
    using namespace std;

    wassert(actual(vsp20.scans.size()) == vodim.scans.size());

    for (unsigned ie = 0; ie < vsp20.scans.size(); ++ie)
    {
        WIBBLE_TEST_INFO(testinfo);
        testinfo() << "elevation " << ie;

        wassert(actual(vsp20.scans[ie].beam_info.size()) == vodim.scans[ie].beam_info.size());
    }
}

void test_volumes_equal(WIBBLE_TEST_LOCPRM, const Volume<double>& vsp20, const Volume<double>& vodim)
{
    using namespace std;

    wassert(actual(vsp20.size()) == vodim.size());

    unsigned failed_beams = 0;
    for (unsigned ie = 0; ie < vsp20.size(); ++ie)
    {
        WIBBLE_TEST_INFO(testinfo);
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

}

template<> template<>
void to::test<1>()
{
    using namespace cumbac::volume;
    // Test loading of a radar volume via SP20
    LoadInfo load_info;
    const Site& gat = Site::get("GAT");
    SP20Loader loader(gat, false, false);
    Scans<double> ssp20;
    loader.vol_z = &ssp20;
    loader.load_info = &load_info;
    loader.load("testdata/DBP2_070120141530_GATTATICO");
    Volume<double> vsp20;
    volume_resample<double>(ssp20, loader.azimuth_maps, vsp20, merger_max_of_closest<double>);
    // Check the contents of what we read
    wruntest(test_0120141530gat_SP20, load_info, vsp20);
}

template<> template<>
void to::test<2>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/MSG1400715300U.101.h5";
    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
    bool res = cb->read_odim_volume(fname, 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // Check the contents of what we read
    wruntest(test_0120141530gat_ODIM, cb->load_info, cb->volume);
    delete cb;
}

template<> template<>
void to::test<3>()
{
    using namespace std;
    using namespace cumbac::volume;
    Scans<double> ssp20;
    Volume<double> vsp20;
    LoadInfo liSP20;
    Scans<double> sodim;
    Volume<double> vodim;
    LoadInfo liODIM;

    // FIXME: get rid of the static elev_array as soon as it is convenient to do so
    const Site& gat = Site::get("GAT");

    SP20Loader sp20(gat, false, false);
    sp20.load_info = &liSP20;
    sp20.vol_z = &ssp20;
    sp20.load("testdata/DBP2_070120141530_GATTATICO");

    ODIMLoader odim(gat, false, false);
    odim.load_info = &liODIM;
    odim.vol_z = &sodim;
    odim.load("testdata/MSG1400715300U.101.h5");

    volume_resample<double>(ssp20, sp20.azimuth_maps, vsp20, merger_max_of_closest<double>);
    volume_resample<double>(sodim, odim.azimuth_maps, vodim, merger_max_of_closest<double>);

    wruntest(test_volumes_equal, vsp20, vodim);
    wruntest(test_loadinfo_equal, liSP20, liODIM);
}

template<> template<>
void to::test<4>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
    bool res = cb->read_sp20_volume(fname, 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // TODO: Check the contents of what we read
    //wruntest(test_0120141530gat, cb->volume);
    delete cb;
}

template<> template<>
void to::test<5>()
{
    using namespace std;
    using namespace cumbac::volume;
    Scans<double> ssp20;
    Volume<double> vsp20;
    LoadInfo liSP20;
    Scans<double> s_mod;
    Volume<double> v_mod;
    LoadInfo li_mod;

    const Site& gat = Site::get("GAT");

    SP20Loader sp20(gat, false, true, 494);
    sp20.load_info = &liSP20;
    sp20.vol_z = &ssp20;
    sp20.load("testdata/DBP2_060220140140_GATTATICO");

    SP20Loader _mod(gat, false, false);
    _mod.load_info = &li_mod;
    _mod.vol_z = &s_mod;
    _mod.load("testdata/DBP2_060220140140_GATTATICO_mod");

    volume_resample<double>(ssp20, sp20.azimuth_maps, vsp20, merger_max_of_closest<double>);
    volume_resample<double>(s_mod, _mod.azimuth_maps, v_mod, merger_max_of_closest<double>);

    wruntest(test_volumes_equal, vsp20, v_mod);
    wruntest(test_loadinfo_equal, liSP20, li_mod);
}


template<> template<>
void to::test<6>()
{

    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/DBP2_020520141110_BOLOGNA";
    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "SPC");
    bool res = cb->read_sp20_volume(fname, 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // Make sure the elevations are what we want
    wassert(actual(cb->volume.size()) == 11);
    wassert(actual(round(cb->volume.scan(0).elevation * 10)) ==    5);
    wassert(actual(round(cb->volume.scan(1).elevation * 10)) ==   14);
    wassert(actual(round(cb->volume.scan(2).elevation * 10)) ==   23);
    wassert(actual(round(cb->volume.scan(3).elevation * 10)) ==   32);
    wassert(actual(round(cb->volume.scan(4).elevation * 10)) ==   41);
    wassert(actual(round(cb->volume.scan(5).elevation * 10)) ==   50);
    wassert(actual(round(cb->volume.scan(6).elevation * 10)) ==   70);
    wassert(actual(round(cb->volume.scan(7).elevation * 10)) ==   95);
    wassert(actual(round(cb->volume.scan(8).elevation * 10)) ==  130);
    wassert(actual(round(cb->volume.scan(9).elevation * 10)) ==  180);
    wassert(actual(round(cb->volume.scan(10).elevation * 10)) == 250);
    // TODO: Check the contents of what we read
    //wruntest(test_0120141530gat, cb->volume);
    delete cb;
}

}
