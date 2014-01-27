#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"

using namespace wibble::tests;

namespace tut {

struct read_sp20_shar {
    Logging logging;
};
TESTGRP(read_sp20);

namespace {

/// Compute some arbitrary stats on a volume
struct VolumeStats
{
    unsigned count_zeros[NEL];
    unsigned count_ones[NEL];
    unsigned count_others[NEL];
    unsigned sum_others[NEL];

    VolumeStats(const CUM_BAC* cb)
    {
        for (int iel = 0; iel < NEL; ++iel)
        {
            count_zeros[iel] = 0;
            count_ones[iel] = 0;
            count_others[iel] = 0;
            sum_others[iel] = 0;

            for (int ibeam = 0; ibeam < cb->nbeam_elev[iel]; ++ibeam)
            {
                for (int i = 0; i < MAX_DIM; ++i)
                {
                    int val = cb->vol_pol[iel][ibeam].ray[i];
                    switch (val)
                    {
                        case 0: count_zeros[iel]++; break;
                        case 1: count_ones[iel]++; break;
                        default:
                                count_others[iel]++;
                                sum_others[iel] += val;
                                break;
                    }
                }
            }
        }
    }
};

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
    // Ensure that nbeam_elev has been filled with the right values
    wassert(actual(cb->nbeam_elev[0]) == 400);
    wassert(actual(cb->nbeam_elev[1]) == 400);
    wassert(actual(cb->nbeam_elev[2]) == 400);
    wassert(actual(cb->nbeam_elev[3]) == 400);
    wassert(actual(cb->nbeam_elev[4]) == 400);
    // Check other header fields
    wassert(actual(cb->old_data_header.norm.maq.acq_date) == 1389108600);

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats(cb);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_ones[0]) == 146674);
    wassert(actual(stats.count_ones[1]) == 184613);
    wassert(actual(stats.count_ones[2]) == 193318);
    wassert(actual(stats.count_ones[3]) == 196292);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 50926);
    wassert(actual(stats.count_others[1]) == 12987);
    wassert(actual(stats.count_others[2]) ==  4282);
    wassert(actual(stats.count_others[3]) ==  1308);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) == 4629202);
    wassert(actual(stats.sum_others[1]) == 890666);
    wassert(actual(stats.sum_others[2]) == 254745);
    wassert(actual(stats.sum_others[3]) == 45968);
    wassert(actual(stats.sum_others[4]) == 78321);
    delete cb;
}

template<> template<>
void to::test<2>()
{
    // Test elabora_dato, con tutti i do_* a false
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats(cb);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_ones[0]) == 192288);
    wassert(actual(stats.count_ones[1]) == 193052);
    wassert(actual(stats.count_ones[2]) == 194509);
    wassert(actual(stats.count_ones[3]) == 196573);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 5312);
    wassert(actual(stats.count_others[1]) == 4548);
    wassert(actual(stats.count_others[2]) == 3091);
    wassert(actual(stats.count_others[3]) == 1027);
    wassert(actual(stats.count_others[4]) == 1440);
    wassert(actual(stats.sum_others[0]) == 493518);
    wassert(actual(stats.sum_others[1]) == 367575);
    wassert(actual(stats.sum_others[2]) == 222441);
    wassert(actual(stats.sum_others[3]) ==  41842);
    wassert(actual(stats.sum_others[4]) ==  78321);

    delete cb;
}

template<> template<>
void to::test<3>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_bloccorr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats(cb);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_ones[0]) == 177082);
    wassert(actual(stats.count_ones[1]) == 190387);
    wassert(actual(stats.count_ones[2]) == 193646);
    wassert(actual(stats.count_ones[3]) == 196318);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 20518);
    wassert(actual(stats.count_others[1]) ==  7213);
    wassert(actual(stats.count_others[2]) ==  3954);
    wassert(actual(stats.count_others[3]) ==  1282);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) == 1525814);
    wassert(actual(stats.sum_others[1]) ==  476555);
    wassert(actual(stats.sum_others[2]) ==  245904);
    wassert(actual(stats.sum_others[3]) ==   45719);
    wassert(actual(stats.sum_others[4]) ==   78321);

    delete cb;
}

}
