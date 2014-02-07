#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#include <stdio.h>
#include <vector>

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct read_sp20_shar {
    Logging logging;
};
TESTGRP(read_sp20);

namespace {

void test_0120141530gat(WIBBLE_TEST_LOCPRM, const Volume& v)
{
    // Ensure that nbeam_elev has been filled with the right values
    wassert(actual(v.nbeam_elev[0]) == 400);
    wassert(actual(v.nbeam_elev[1]) == 400);
    wassert(actual(v.nbeam_elev[2]) == 400);
    wassert(actual(v.nbeam_elev[3]) == 400);
    wassert(actual(v.nbeam_elev[4]) == 400);
    wassert(actual(v.nbeam_elev[5]) == 400);
    wassert(actual(v.nbeam_elev[6]) == 0);
    wassert(actual(v.nbeam_elev[7]) == 0);

    // Ensure that the beam sizes are what we expect
    wassert(actual(v.vol_pol[0][0].ray.size()) == 494);
    wassert(actual(v.vol_pol[1][0].ray.size()) == 494);
    wassert(actual(v.vol_pol[2][0].ray.size()) == 494);
    wassert(actual(v.vol_pol[3][0].ray.size()) == 494);
    wassert(actual(v.vol_pol[4][0].ray.size()) == 494);
    wassert(actual(v.vol_pol[5][0].ray.size()) == 494);

    // Ensure that the beam azimuth are what we expect
    wassert(actual(v.vol_pol[0][0].alfa) == 0);
    wassert(actual(v.vol_pol[0][1].alfa) == 10);
    wassert(actual(v.vol_pol[1][1].alfa) == 10);
    wassert(actual(v.vol_pol[2][1].alfa) == 10);
    wassert(actual(v.vol_pol[3][1].alfa) == 10);
    wassert(actual(v.vol_pol[4][1].alfa) == 10);
    wassert(actual(v.vol_pol[5][1].alfa) == 10);

    // Check other header fields
    wassert(actual(v.acq_date) == 1389108600);
    wassert(actual(v.size_cell) == 250);

    // for (int i = 0; i < 200; ++i)
    //     printf("%d ", (int)v.vol_pol[0][0].ray[i]);
    // printf("\n");

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats;
    v.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_zeros[5]) == 0);
    wassert(actual(stats.count_zeros[6]) == 0);
    wassert(actual(stats.count_ones[0]) == 146674);
    wassert(actual(stats.count_ones[1]) == 184613);
    wassert(actual(stats.count_ones[2]) == 193318);
    wassert(actual(stats.count_ones[3]) == 196292);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_ones[5]) == 196158);
    wassert(actual(stats.count_ones[6]) == 0);
    wassert(actual(stats.count_others[0]) == 50926);
    wassert(actual(stats.count_others[1]) == 12987);
    wassert(actual(stats.count_others[2]) ==  4282);
    wassert(actual(stats.count_others[3]) ==  1308);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.count_others[5]) ==  1442);
    wassert(actual(stats.count_others[6]) ==     0);
    wassert(actual(stats.sum_others[0]) == 4629202);
    wassert(actual(stats.sum_others[1]) == 890666);
    wassert(actual(stats.sum_others[2]) == 254745);
    wassert(actual(stats.sum_others[3]) == 45968);
    wassert(actual(stats.sum_others[4]) == 78321);
    wassert(actual(stats.sum_others[5]) == 88234);
    wassert(actual(stats.sum_others[6]) == 0);
}

}

template<> template<>
void to::test<1>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    CUM_BAC* cb = new CUM_BAC;
    bool res = cb->read_sp20_volume(fname, "GAT", 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // Check the contents of what we read
    wruntest(test_0120141530gat, cb->volume);
    delete cb;
}

template<> template<>
void to::test<2>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/MSG1400715300U.101.h5";
    CUM_BAC* cb = new CUM_BAC;
    bool res = cb->read_odim_volume(fname, "GAT", 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // Check the contents of what we read
    wruntest(test_0120141530gat, cb->volume);
    delete cb;
}

template<> template<>
void to::test<3>()
{
    using namespace std;
    Volume vsp20;
    Volume vodim;

    // FIXME: get rid of the static elev_array as soon as it is convenient to do so
    for (unsigned i = 0; i < NEL; ++i)
        elev_array[i] = elev_array_gat[i];

    vsp20.read_sp20("testdata/DBP2_070120141530_GATTATICO");
    vodim.read_odim("testdata/MSG1400715300U.101.h5");

    unsigned failed_beams = 0;
    for (unsigned ie = 0; ie < NEL; ++ie)
    {
        WIBBLE_TEST_INFO(testinfo);
        testinfo() << "elevation " << ie;

        wassert(actual(vsp20.nbeam_elev[ie]) == vodim.nbeam_elev[ie]);
        if (vsp20.nbeam_elev[ie] == 0) continue;

        for (unsigned ia = 0; ia < vsp20.nbeam_elev[ie]; ++ia)
        {
            testinfo() << "elevation " << ie << " angle " << ia;

            wassert(actual(vsp20.vol_pol[ie][ia].ray.size()) == vodim.vol_pol[ie][ia].ray.size());

            vector<unsigned char> vals_sp20;
            vector<unsigned char> vals_odim;
            for (unsigned ib = 0; ib < vsp20.vol_pol[ie][ia].ray.size(); ++ib)
            {
                if (vsp20.vol_pol[ie][ia].ray[ib] != vodim.vol_pol[ie][ia].ray[ib])
                {
                    vals_sp20.push_back(vsp20.vol_pol[ie][ia].ray[ib]);
                    vals_odim.push_back(vodim.vol_pol[ie][ia].ray[ib]);
                }
            }
            if (!vals_sp20.empty())
            {
                printf("sp20 vp[%u][%u]:", ie, ia);
                for (vector<unsigned char>::const_iterator i = vals_sp20.begin(); i != vals_sp20.end(); ++i)
                    printf(" %d", (int)*i);
                printf("\n");
                printf("odim vp[%u][%u]:", ie, ia);
                for (vector<unsigned char>::const_iterator i = vals_odim.begin(); i != vals_odim.end(); ++i)
                    printf(" %d", (int)*i);
                printf("\n");
                printf("sp20 vp[%u][%u] load log: ", ie, ia); vsp20.vol_pol[ie][ia].print_load_log(stdout);
                printf("odim vp[%u][%u] load log: ", ie, ia); vodim.vol_pol[ie][ia].print_load_log(stdout);
                ++failed_beams;
            }
        }
    }

    wassert(actual(failed_beams) == 0);
}

template<> template<>
void to::test<4>()
{
    // Test loading of a radar volume via SP20
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    CUM_BAC* cb = new CUM_BAC;
    bool res = cb->read_sp20_volume(fname, "GAT", 0);
    // Ensure that reading was successful
    wassert(actual(res).istrue());
    // Check the contents of what we read
    //wruntest(test_0120141530gat, cb->volume);
    delete cb;
}


}
