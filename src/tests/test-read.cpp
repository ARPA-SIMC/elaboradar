#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#include <stdio.h>

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

    VolumeStats(const Volume& v)
    {
        for (int iel = 0; iel < NEL; ++iel)
        {
            count_zeros[iel] = 0;
            count_ones[iel] = 0;
            count_others[iel] = 0;
            sum_others[iel] = 0;

            for (int ibeam = 0; ibeam < v.nbeam_elev[iel]; ++ibeam)
            {
                for (int i = 0; i < MAX_DIM; ++i)
                {
                    int val = v.vol_pol[iel][ibeam].ray[i];
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
    void PrintStats(){	
	printf("Nel   Zeros    Ones  Others     Sum\n"); 
	for (int iel =0; iel<NEL; ++iel){
	    printf("%4d%8d%8d%8d%8d\n",iel,count_zeros[iel],count_ones[iel],count_others[iel],sum_others[iel]);
	}
    }
};

template<typename T>
struct ArrayStats
{
    bool first;
    T min;
    T max;
    double avg;

    ArrayStats() : first(true) {}

    template<int A, int B>
    void fill2(const T (&arr)[A][B])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
            {
                if (first)
                {
                    min = arr[i][j];
                    max = arr[i][j];
                    avg = min / (A * B);
                    first = false;
                }
                else
                {
                    if (arr[i][j] < min)
                        min = arr[i][j];
                    if (arr[i][j] > max)
                        max = arr[i][j];
                    avg += (double)arr[i][j] / (A * B);
                }
            }
    }

    template<int A, int B, int C>
    void fill3(const T (&arr)[A][B][C])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                for (int k = 0; k < C; ++k)
                {
                    if (first)
                    {
                        min = arr[i][j][k];
                        max = arr[i][j][k];
                        avg = min / (A * B * C);
                        first = false;
                    }
                    else
                    {
                        if (arr[i][j][k] < min)
                            min = arr[i][j][k];
                        if (arr[i][j][k] > max)
                            max = arr[i][j][k];
                        avg += (double)arr[i][j][k] / (A * B * C);
                    }
                }
    }
};

void test_0120141530gat_odim(WIBBLE_TEST_LOCPRM, const Volume& v)
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
    wassert(actual(v.vol_pol[0][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[1][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[2][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[3][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[4][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[5][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[6][0].b_header.max_bin) == 0);

    // Ensure that the beam azimuth are what we expect
    wassert(actual(v.vol_pol[0][0].b_header.alfa) == 0);
    wassert(actual(v.vol_pol[0][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[1][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[2][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[3][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[4][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[5][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[6][1].b_header.alfa) == 0);

    // Check other header fields
    wassert(actual(v.acq_date) == 1389108600);
    wassert(actual(v.size_cell) == 250);

    // for (int i = 0; i < 200; ++i)
    //     printf("%d ", (int)v.vol_pol[0][0].ray[i]);
    // printf("\n");

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats(v);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_zeros[5]) == 7200);
    wassert(actual(stats.count_zeros[6]) == 0);
    wassert(actual(stats.count_ones[0]) == 145970);
    wassert(actual(stats.count_ones[1]) == 184986);
    wassert(actual(stats.count_ones[2]) == 193941);
    wassert(actual(stats.count_ones[3]) == 196293);
    wassert(actual(stats.count_ones[4]) == 196158);
    wassert(actual(stats.count_ones[5]) == 196158);
    wassert(actual(stats.count_ones[6]) == 0);
    wassert(actual(stats.count_others[0]) == 51630);
    wassert(actual(stats.count_others[1]) == 12614);
    wassert(actual(stats.count_others[2]) ==  3659);
    wassert(actual(stats.count_others[3]) ==  1307);
    wassert(actual(stats.count_others[4]) ==  1442);
    wassert(actual(stats.count_others[5]) ==  1442);
    wassert(actual(stats.count_others[6]) ==     0);
    wassert(actual(stats.sum_others[0]) == 4663175);
    wassert(actual(stats.sum_others[1]) == 858438);
    wassert(actual(stats.sum_others[2]) == 206212);
    wassert(actual(stats.sum_others[3]) == 45942);
    wassert(actual(stats.sum_others[4]) == 78314);
    wassert(actual(stats.sum_others[5]) == 88163);
    wassert(actual(stats.sum_others[6]) == 0);
}

void test_0120141530gat_sp20(WIBBLE_TEST_LOCPRM, const Volume& v)
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
    wassert(actual(v.vol_pol[0][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[1][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[2][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[3][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[4][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[5][0].b_header.max_bin) == 494);
    wassert(actual(v.vol_pol[6][0].b_header.max_bin) == 0);

    // Ensure that the beam azimuth are what we expect
    wassert(actual(v.vol_pol[0][0].b_header.alfa) == 0);
    wassert(actual(v.vol_pol[0][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[1][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[2][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[3][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[4][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[5][1].b_header.alfa) == 10);
    wassert(actual(v.vol_pol[6][1].b_header.alfa) == 0);

    // Check other header fields
    wassert(actual(v.acq_date) == 1389108600);
    wassert(actual(v.size_cell) == 250);

    // for (int i = 0; i < 200; ++i)
    //     printf("%d ", (int)v.vol_pol[0][0].ray[i]);
    // printf("\n");

    // Arbitrary stats on volume contents so we can check that we read data
    // that looks correct
    VolumeStats stats(v);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_zeros[5]) == 7200);
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
    wruntest(test_0120141530gat_sp20, cb->volume);
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
    wruntest(test_0120141530gat_odim, cb->volume);
    delete cb;
}

template<> template<>
void to::test<3>()
{
    // Test elabora_dato, con tutti i do_* a false
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    wassert(actual(cb->t_ground) == NODATAVPR);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats(cb->volume);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_ones[0]) == 192288);
    wassert(actual(stats.count_ones[1]) == 193052);
    wassert(actual(stats.count_ones[2]) == 194525);
    wassert(actual(stats.count_ones[3]) == 196576);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 5312);
    wassert(actual(stats.count_others[1]) == 4548);
    wassert(actual(stats.count_others[2]) == 3075);
    wassert(actual(stats.count_others[3]) == 1024);
    wassert(actual(stats.count_others[4]) == 1440);
    wassert(actual(stats.sum_others[0]) == 493518);
    wassert(actual(stats.sum_others[1]) == 358463);
    wassert(actual(stats.sum_others[2]) == 220768);
    wassert(actual(stats.sum_others[3]) ==  41677);
    wassert(actual(stats.sum_others[4]) ==  78321);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    qual_stats.fill3(cb->qual);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 0);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 4234);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill3(cb->flag_vpr);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 0);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 0);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill2(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 35);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 5);

    delete cb;
}

template<> template<>
void to::test<4>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->do_class = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    wassert(actual((int)(cb->t_ground * 100)) == 1010);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats(cb->volume);
    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_ones[0]) == 176003);
    wassert(actual(stats.count_ones[1]) == 190293);
    wassert(actual(stats.count_ones[2]) == 193588);
    wassert(actual(stats.count_ones[3]) == 196292);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 21597);
    wassert(actual(stats.count_others[1]) ==  7307);
    wassert(actual(stats.count_others[2]) ==  4012);
    wassert(actual(stats.count_others[3]) ==  1308);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) == 1538560);
    wassert(actual(stats.sum_others[1]) ==  478421);
    wassert(actual(stats.sum_others[2]) ==  246658);
    wassert(actual(stats.sum_others[3]) ==   45968);
    wassert(actual(stats.sum_others[4]) ==   78321);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    qual_stats.fill3(cb->qual);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 0);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 5955);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill3(cb->flag_vpr);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 1);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 92);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill2(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 35);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 10);

    cb->classifica_rain();

    ArrayStats<unsigned char> stratiform_stats;
    stratiform_stats.fill2(cb->stratiform);
    wassert(actual((unsigned)stratiform_stats.first).isfalse());
    wassert(actual((unsigned)stratiform_stats.min) == 0);
    wassert(actual((unsigned)stratiform_stats.max) == 0);
    wassert(actual((unsigned)(stratiform_stats.avg * 100)) == 0);

    delete cb;
}

template<> template<>
void to::test<5>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats(cb->volume);

    wassert(actual(stats.count_zeros[0]) == 7200);
    wassert(actual(stats.count_zeros[1]) == 7200);
    wassert(actual(stats.count_zeros[2]) == 7200);
    wassert(actual(stats.count_zeros[3]) == 7200);
    wassert(actual(stats.count_zeros[4]) == 7200);
    wassert(actual(stats.count_ones[0]) == 193092);
    wassert(actual(stats.count_ones[1]) == 193930);
    wassert(actual(stats.count_ones[2]) == 194527);
    wassert(actual(stats.count_ones[3]) == 196577);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) ==  4508);
    wassert(actual(stats.count_others[1]) ==  3670);
    wassert(actual(stats.count_others[2]) ==  3073);
    wassert(actual(stats.count_others[3]) ==  1023);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) ==  385721);
    wassert(actual(stats.sum_others[1]) ==  276535);
    wassert(actual(stats.sum_others[2]) ==  220710);
    wassert(actual(stats.sum_others[3]) ==   41674);
    wassert(actual(stats.sum_others[4]) ==   78321);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    qual_stats.fill3(cb->qual);
    wassert(actual((unsigned)qual_stats.min) == 0);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 5874);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill3(cb->flag_vpr);
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 1);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 91);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill2(cb->top);
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 35);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 3);

    delete cb;
}


}
