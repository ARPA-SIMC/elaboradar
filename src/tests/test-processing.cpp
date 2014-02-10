#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct process_shar {
    Logging logging;
};
TESTGRP(process);

namespace {

template<typename T>
struct ArrayStats
{
    bool first;
    T min;
    T max;
    double avg;
    unsigned count_zeros;
    unsigned count_ones;

    ArrayStats() : first(true), count_zeros(0), count_ones(0) {}

    void count_sample(const T& sample, unsigned item_count)
    {
        if (sample == 0)
            ++count_zeros;
        else if (sample == 1)
            ++count_ones;

        if (first)
        {
            min = sample;
            max = sample;
            avg = (double)min / item_count;
            first = false;
        }
        else
        {
            if (sample < min)
                min = sample;
            if (sample > max)
                max = sample;
            avg += (double)sample / item_count;
        }
    }

    template<int A, int B>
    void fill2(const T (&arr)[A][B])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                count_sample(arr[i][j], A * B);
    }

    template<int A, int B, int C>
    void fill3(const T (&arr)[A][B][C])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                for (int k = 0; k < C; ++k)
                    count_sample(arr[i][j][k], A * B * C);
    }
};

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

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    wassert(actual(cb->calcolo_vpr->t_ground) == NODATAVPR);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
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
    vpr_stats.fill3(cb->calcolo_vpr->flag_vpr);
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

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->do_class = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    wassert(actual((int)(cb->calcolo_vpr->t_ground * 100)) == 1010);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
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
    vpr_stats.fill3(cb->calcolo_vpr->flag_vpr);
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

    cb->calcolo_vpr->classifica_rain();

    ArrayStats<unsigned char> stratiform_stats;
    stratiform_stats.fill2(cb->calcolo_vpr->stratiform);
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

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);

    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 193091);
    wassert(actual(stats.count_ones[1]) == 193929);
    wassert(actual(stats.count_ones[2]) == 194526);
    wassert(actual(stats.count_ones[3]) == 196576);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) ==  4509);
    wassert(actual(stats.count_others[1]) ==  3671);
    wassert(actual(stats.count_others[2]) ==  3074);
    wassert(actual(stats.count_others[3]) ==  1024);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) ==  387010);
    wassert(actual(stats.sum_others[1]) ==  276931);
    wassert(actual(stats.sum_others[2]) ==  220729);
    wassert(actual(stats.sum_others[3]) ==   41677);
    wassert(actual(stats.sum_others[4]) ==   78321);

    delete cb;
}

template<> template<>
void to::test<6>()
{
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = false;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
    // in hray_inf
    // out hray
    // in dato_corrotto
    // in beam_blocking
    // in elev_fin
    // out qual
    // out flag_vpr
    // out top
    ArrayStats<unsigned char> stats_qual;
    stats_qual.fill3(cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 108000);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 0);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5907);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill3(cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 3072000);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 0);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 0);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill2(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 203210);
    wassert(actual((unsigned)stats_top.count_ones) == 4);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 35);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 10);

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 10);
    wassert(actual((unsigned)cb->cart.max()) == 226);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 0);
    wassert(actual((unsigned)cb->topxy.max()) == 35);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 0);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 16);
    wassert(actual((unsigned)cb->z_out.max()) == 226);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 29);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->top_1x1.max()) == 35);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    delete cb;
}

template<> template<>
void to::test<7>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
    ArrayStats<unsigned char> stats_qual;
    stats_qual.fill3(cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 108000);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 0);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5907);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill3(cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 251320);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 2820680);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 91);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill2(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 203210);
    wassert(actual((unsigned)stats_top.count_ones) == 4);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 35);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 10);

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 10);
    wassert(actual((unsigned)cb->cart.max()) == 226);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 0);
    wassert(actual((unsigned)cb->topxy.max()) == 35);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 1);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 16);
    wassert(actual((unsigned)cb->z_out.max()) == 226);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 29);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->top_1x1.max()) == 35);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 1);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 192);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    delete cb;
}

template<> template<>
void to::test<8>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
    ArrayStats<unsigned char> stats_qual;
    stats_qual.fill3(cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 108000);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 0);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5907);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill3(cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 251320);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 2820680);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 91);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill2(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 203210);
    wassert(actual((unsigned)stats_top.count_ones) == 4);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 35);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 10);

    cb->calcolo_vpr->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 1); // TODO: cosa deve dare?

    // TODO: cb->stampa_vpr()

    cb->class_conv_fixme_find_a_name();

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 10);
    wassert(actual((unsigned)cb->cart.max()) == 226);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 0);
    wassert(actual((unsigned)cb->topxy.max()) == 35);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 1);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 16);
    wassert(actual((unsigned)cb->z_out.max()) == 226);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 29);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->top_1x1.max()) == 35);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 1);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 192);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<9>()
{
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/ultimo_vpr");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = false;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    ArrayStats<unsigned char> stats_qual;
    stats_qual.fill3(cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 162505);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5888);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill3(cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 3072000);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 0);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 0);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill2(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 111644);
    wassert(actual((unsigned)stats_top.count_ones) == 423);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 75);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 565);

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 55);
    wassert(actual((unsigned)cb->cart.max()) == 255);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 3);
    wassert(actual((unsigned)cb->topxy.max()) == 75);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 26);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 14);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 0);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 66);
    wassert(actual((unsigned)cb->z_out.max()) == 255);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 25);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 97);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 15);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 3);
    wassert(actual((unsigned)cb->top_1x1.max()) == 75);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<10>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO_mod";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/ultimo_vpr");
    setenv("LAST_VPR", "testdata/last_vpr", 1);
    unlink("testdata/last_vpr");
    setenv("VPR_HMAX", "testdata/vpr_hmax", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    std::cout<<"INIZIO - "<<fname<<std::endl;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
std::cout<<"Dopo caratterizzo volume "<<std::endl;
    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    cb->calcolo_vpr->test_vpr=fopen("testdata/test_vpr","a+");
    ier = cb->calcolo_vpr->combina_profili();
std::cout<<"Dopo combina volume"<<std::endl;
    wassert(actual(ier) == 0);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 0);

    // TODO: cb->stampa_vpr()

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<11>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO_mod";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/ultimo_vpr");
    setenv("LAST_VPR", "testdata/last_vpr", 1);
    unlink("testdata/last_vpr");
    setenv("VPR_HMAX", "testdata/vpr_hmax", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();

    cb->calcolo_vpr->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    cb->calcolo_vpr->test_vpr=fopen("testdata/test_vpr","a+");
    ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 0);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 0);

    // TODO: cb->stampa_vpr()

    cb->class_conv_fixme_find_a_name();

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

#if 0
template<> template<>
void to::test<6>()
{
    // Test di tutto quello che può essere chiamato
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("LAST_VPR", "testdata/last_vpr_GAT", 1);
    setenv("VPR0_FILE", "testdata/vpr_SPC", 1);
    setenv("VPR_HEATING", "vpr_heat_GAT", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

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

    cb->classifica_rain();

    ier = cb->combina_profili();
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    delete cb;
}
#endif

}
